import os
import gc
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torch_geometric.data import Batch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score

from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies

def main(m_matrix_filepath, manifest_filepath, results_filepath):
    print("Starting main function")

    m_matrix_full = pd.read_parquet(m_matrix_filepath)
    manifest_df = pd.read_parquet(manifest_filepath)

    all_cols = m_matrix_full.columns.to_list()
    probe_cols = [c for c in all_cols if c.startswith('cg') or (c.startswith('ch') and not c.startswith('cha'))]
    probe_set = set(probe_cols)
    pheno_cols = [c for c in all_cols if c not in probe_set]
    pheno_df = m_matrix_full[pheno_cols].copy()

    cat_cols = pheno_df.select_dtypes(include=["category"]).columns
    for c in cat_cols:
        pheno_df[c] = pheno_df[c].astype("str")

    print(f"Writing pheno data to {results_filepath}")
    pheno_out = f"{results_filepath}/pheno_data.parquet"
    pheno_df.to_parquet(pheno_out)
    print(f"Wrote pheno_df to {pheno_out}")

    m_matrix_with_id = (
        m_matrix_full.copy()
        if 'Sample_Name' in m_matrix_full.columns
        else m_matrix_full.reset_index().rename(columns={m_matrix_full.index.name or 'index': 'Sample_Name'})
    )

    cols_to_keep = ['Sample_Name'] + [
        c for c in m_matrix_with_id.columns
        if c != 'Sample_Name' and c not in pheno_df.columns
    ]

    m_matrix_full_reduced = m_matrix_with_id[cols_to_keep]
    output_path = os.path.join(results_filepath, "m_matrix_full_reduced.parquet")
    m_matrix_full_reduced.to_parquet(output_path, index=False)
    print(f"Saved to: {output_path}")

    LABEL_MAP = {'Control': 0, 'PD': 1}

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Executing GAT Pipeline on: {device}")

    if pheno_df.index.name != "Sample_Name":
        pheno_df = pheno_df.set_index("Sample_Name")

    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']

    common_probes = list(set(m_matrix_full_reduced.columns).intersection(set(manifest_df['IlmnID'])))
    m_matrix_df = m_matrix_full_reduced[common_probes] # Discard unused columns to save RAM early
    print(f"Total overlapping autosome probes: {len(common_probes):,}")

    m_matrix_df.set_index("Sample_Name", inplace=True)

    print("Building chromosome 1D adjacency and genic graphs...")
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)

    y = pheno_df['Sample_Group'].map(LABEL_MAP)
    if y.isna().any():
        unknown_labels = sorted(pheno_df.loc[y.isna(), 'Sample_Group'].astype(str).unique().tolist())
        raise ValueError(f"Unexpected Sample_Group values: {unknown_labels}")
    y = y.to_numpy(dtype=np.int64)

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    roc_aucs, pr_aucs = [], []
    batch_size = 16  
    epochs = 40
    lr = 2e-4

    for fold, (train_idx, test_idx) in enumerate(skf.split(pheno_df, y)):
        print(f"\n--- Fold {fold + 1} ---")
        train_pheno = pheno_df.iloc[train_idx]
        test_pheno = pheno_df.iloc[test_idx]
        
        train_ds = WholeBloodMethylationDataset(m_matrix_df, train_pheno, cell_cols, chr_topologies)
        test_ds = WholeBloodMethylationDataset(m_matrix_df, test_pheno, cell_cols, chr_topologies)
        
        train_loader = DataLoader(train_ds, 
                                batch_size=batch_size, 
                                shuffle=True, 
                                collate_fn=chromosome_collate_fn,
                                num_workers=4,
                                pin_memory=True)
        test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False, collate_fn=chromosome_collate_fn)
        
        model = ChromosomeParallelGAT().to(device)
        
        # Explicit dtype declaration for pos_weight
        train_labels = train_pheno['Sample_Group'].map(LABEL_MAP)
        if train_labels.isna().any():
            unknown_labels = sorted(train_pheno.loc[train_labels.isna(), 'Sample_Group'].astype(str).unique().tolist())
            raise ValueError(f"Unexpected Sample_Group values in training fold: {unknown_labels}")
        train_labels = train_labels.to_numpy(dtype=np.float32)
        positives = float((train_labels == 1.0).sum())
        negatives = float((train_labels == 0.0).sum())
        if positives == 0:
            raise ValueError("Training fold has no positive samples, cannot compute pos_weight")
        pos_weight = torch.tensor([negatives / positives], device=device, dtype=torch.float32)
        criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-3)
        
        # Training
        model.train()
        for epoch in range(epochs):
            total_loss = 0.0
            for batched_chrs, u_cells, batch_y in train_loader:
                
                # Asynchronous memory transfer
                batched_chrs = [c.to(device, non_blocking=True) for c in batched_chrs]
                u_cells = u_cells.to(device, non_blocking=True)
                batch_y = batch_y.to(device, non_blocking=True)
                
                optimizer.zero_grad()
                logits = model(batched_chrs, u_cells).squeeze(1)
                loss = criterion(logits, batch_y)
                loss.backward()
                optimizer.step()
                
                total_loss += loss.item()
                
                # VRAM Fragmentation Protection
                del batched_chrs, u_cells, batch_y, logits, loss
                
        # Evaluation
        model.eval()
        preds, truths = [], []
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in test_loader:
                batched_chrs = [c.to(device, non_blocking=True) for c in batched_chrs]
                u_cells = u_cells.to(device, non_blocking=True)
                
                logits = model(batched_chrs, u_cells).squeeze(1)
                probs = torch.sigmoid(logits).cpu().numpy()
                
                preds.extend(probs)
                truths.extend(batch_y.numpy())
                
                del batched_chrs, u_cells, logits
                
        fold_roc = roc_auc_score(truths, preds)
        fold_pr = average_precision_score(truths, preds)
        roc_aucs.append(fold_roc)
        pr_aucs.append(fold_pr)
        print(f"Fold {fold + 1} | ROC AUC: {fold_roc:.4f} | PR AUC: {fold_pr:.4f}")

        # Save the model weights for downstream biological extraction
        torch.save(model.state_dict(), f"/workspace/results/gat_fold_{fold + 1}.pt")

        del model, optimizer, train_ds, test_ds
        torch.cuda.empty_cache()
        gc.collect()
        
    print("\n" + "=" * 40)
    print(f"Mean GAT ROC AUC: {np.mean(roc_aucs):.4f} ± {np.std(roc_aucs):.4f}")
    print(f"Mean GAT PR AUC:  {np.mean(pr_aucs):.4f} ± {np.std(pr_aucs):.4f}")



if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])
