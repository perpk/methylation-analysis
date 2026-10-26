import os
import gc
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import torch.multiprocessing

from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies

def evaluate_cohorts():
    torch.multiprocessing.set_sharing_strategy('file_system')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Executing GAT Evaluation Pipeline on: {device}")

    # =====================================================================
    # CONFIGURATION
    # =====================================================================
    # Adjust this path to wherever your global array manifest is stored
    manifest_filepath = "/workspace/data/infinium450k_manifest.parquet" 
    
    cohorts = {
        "PEG1": "/workspace/data/peg1",
        "SGPD": "/workspace/data/sgpd",
        "PPMI": "/workspace/data/ppmi"
    }

    cohort_m_matrix = {
        "PEG1": "GSE111629_data_test.parquet",
        "SGPD": "GSE145361_data_corrected.parquet",
        "PPMI": "ppmi_data_test.parquet"
    }

    output_dir = "/workspace/results/confusion_matrices/"
    os.makedirs(output_dir, exist_ok=True)

    LABEL_MAP = {'Control': 0, 'PD': 1}
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']
    
    # Load manifest once globally
    manifest_df = pd.read_parquet(manifest_filepath)

    # =====================================================================
    # COHORT LOOP
    # =====================================================================
    for cohort_name, cohort_path in cohorts.items():
        print(f"\n{'='*60}\nEvaluating Cohort: {cohort_name}\n{'='*60}")
        
        pheno_path = os.path.join(cohort_path, "pheno_data.parquet")
        m_matrix_path = os.path.join(cohort_path, "m_matrix_full_reduced.parquet")
        
        if not os.path.exists(pheno_path) or not os.path.exists(m_matrix_path):
            print(f"Skipping {cohort_name}: Processed data files not found in {cohort_path}")
            continue

        pheno_df = pd.read_parquet(pheno_path)
        m_matrix_df = pd.read_parquet(m_matrix_path)
        
        # FIX 1: Rename Neu to Gran BEFORE feeding to the dataset!
        if 'Neu' in pheno_df.columns:
            pheno_df = pheno_df.rename(columns={'Neu': 'Gran'})
            
        if pheno_df.index.name != "Sample_Name":
            pheno_df = pheno_df.set_index("Sample_Name")
        m_matrix_df.set_index('Sample_Name', inplace=True)
        
        # FIX 2: Enforce deterministic sorting to fix the `set()` randomization bug
        # By sorting alphabetically, we guarantee the column order matches perfectly 
        # with how PyTorch Geometric builds the node topologies.
        raw_intersection = set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID']))
        common_probes = sorted(list(raw_intersection)) # <-- The critical fix
        
        m_matrix_df = m_matrix_df[common_probes]
        print(f"Total overlapping autosome probes (Sorted): {len(common_probes):,}")
        
        print("Building chromosome topologies...")
        chr_topologies = build_chromosome_topologies(manifest_df, common_probes)
        
        # 3. Setup the Exact Cross-Validation Splitter
        y = pheno_df['Sample_Group'].map(LABEL_MAP)
        
        # Safety Check: Did the mapping fail?
        if y.isna().any():
            raise ValueError(f"Label Mapping failed! Check for NaNs. Unique values: {pheno_df['Sample_Group'].unique()}")
            
        y = y.to_numpy(dtype=np.int64)
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        all_oof_probs = []
        all_oof_targets = []
        
        # Instantiate model architecture and move to device
        model = ChromosomeParallelGAT().to(device)
        
        # 4. Iterate over the folds
        for fold, (train_idx, test_idx) in enumerate(skf.split(pheno_df, y)):
            print(f"  -> Processing Fold {fold + 1}...")
            
            checkpoint_path = f"{cohort_path}/gat_fold_{fold + 1}.pt"
            if not os.path.exists(checkpoint_path):
                print(f"     [!] Checkpoint not found: {checkpoint_path}")
                continue
                
            # Load weights
            checkpoint = torch.load(checkpoint_path, map_location=device)
            model.load_state_dict(checkpoint['model_state_dict'])
            model.eval() # Disable dropout/batchnorm for inference
            
            # Recreate test subset exactly as done in training
            test_pheno = pheno_df.iloc[test_idx]
            test_ds = WholeBloodMethylationDataset(m_matrix_df, test_pheno, cell_cols, chr_topologies)
            
            test_loader = DataLoader(
                test_ds, 
                batch_size=16, # Matching your training script
                shuffle=False, 
                collate_fn=chromosome_collate_fn,
                num_workers=6,
                pin_memory=True
            )

            # 5. Inference Loop
            with torch.no_grad():
                for batched_chrs, u_cells, batch_y in test_loader:
                    # Move to GPU via non_blocking to match your training style
                    batched_chrs = [c.to(device, non_blocking=True) for c in batched_chrs]
                    u_cells = u_cells.to(device, non_blocking=True)
                    
                    logits = model(batched_chrs, u_cells).squeeze(1)
                    probs = torch.sigmoid(logits).cpu().numpy()
                    
                    if probs.ndim == 0:
                        probs = np.expand_dims(probs, 0)
                        
                    all_oof_probs.extend(probs)
                    all_oof_targets.extend(batch_y.numpy())
                    
                    # Prevent VRAM spikes
                    del batched_chrs, u_cells, logits
                    
            # Clear memory between folds
            torch.cuda.empty_cache()
            gc.collect()

        if len(all_oof_targets) == 0:
            print(f"No validation data processed for {cohort_name}.")
            continue

        # =====================================================================
        # METRICS & PLOTTING 
        # =====================================================================
        y_true = np.array(all_oof_targets)
        y_probs = np.array(all_oof_probs)

        # Optimize Threshold via F1 Score
        precision, recall, thresholds = precision_recall_curve(y_true, y_probs)
        f1_scores = (2 * precision * recall) / (precision + recall + 1e-8)
        optimal_idx = np.argmax(f1_scores)
        optimal_threshold = thresholds[optimal_idx]

        y_preds = (y_probs >= optimal_threshold).astype(int)
        cm = confusion_matrix(y_true, y_preds)
        tn, fp, fn, tp = cm.ravel()
        
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        auc_score = roc_auc_score(y_true, y_probs)

        print(f"\n     --- {cohort_name} Aggregated Metrics ---")
        print(f"     Optimal Threshold : {optimal_threshold:.4f}")
        print(f"     Global ROC AUC    : {auc_score:.4f}")
        print(f"     Sensitivity       : {sensitivity:.1%}")
        print(f"     Specificity       : {specificity:.1%}\n")

        # Plotting
        plt.figure(figsize=(7, 6))
        cm_percentages = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        
        labels = np.asarray([
            f"True Negative\n{tn}\n({cm_percentages[0, 0]:.1%})",
            f"False Positive\n{fp}\n({cm_percentages[0, 1]:.1%})",
            f"False Negative\n{fn}\n({cm_percentages[1, 0]:.1%})",
            f"True Positive\n{tp}\n({cm_percentages[1, 1]:.1%})"
        ]).reshape(2, 2)

        ax = sns.heatmap(
            cm, annot=labels, fmt="", cmap='Blues', cbar=False, square=True,
            annot_kws={"size": 13, "weight": "bold"}, linewidths=2, linecolor='black'
        )

        ax.set_title(f'Aggregated Validation Performance\nCohort: {cohort_name}', fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel('Predicted Diagnosis (GATv2)', fontsize=14, fontweight='bold', labelpad=15)
        ax.set_ylabel('True Clinical Diagnosis', fontsize=14, fontweight='bold', labelpad=15)
        ax.set_xticklabels(['Control (0)', 'Parkinson\'s (1)'], fontsize=12)
        ax.set_yticklabels(['Control (0)', 'Parkinson\'s (1)'], fontsize=12, rotation=0, va='center')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"Fig_CM_{cohort_name}.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f"Fig_CM_{cohort_name}.pdf"), dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    evaluate_cohorts()