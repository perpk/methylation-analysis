import pandas as pd
from sklearn.model_selection import StratifiedKFold
import torch
from torch.utils.data import DataLoader
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies

cohorts = {
    "peg1": {
        "results_path": "/workspace/results/peg1",
        "pheno_data": "/workspace/data/peg1/GSE111629_pheno_data.parquet",
        "m_matrix": "/workspace/data/peg1/GSE111629_m_matrix_full_reduced.parquet"
    },
    "sgpd": {
        "results_path": "/workspace/results/sgpd",
        "pheno_data": "/workspace/data/sgpd/pheno_data.parquet",
        "m_matrix": "/workspace/data/sgpd/GSE145361_data_corrected.parquet"
    },
    "ppmi": {
        "results_path": "/workspace/results/ppmi",
        "pheno_data": "/workspace/data/ppmi/ppmi_pheno_data.parquet",
        "m_matrix": "/workspace/data/ppmi/m_matrix_full_reduced.parquet"
    }
}

manifest_df = pd.read_parquet("/workspace/data/infinium450k_manifest.parquet")

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for cohort, path in cohorts.items():
    results_path = path["results_path"]
    master_pheno_df = pd.read_parquet(path["pheno_data"])
    master_m_matrix = pd.read_parquet(path["m_matrix"])

    fold_splits = list(skf.split(master_pheno_df, master_pheno_df['Sample_Group']))

    fold_idx = 0 
    train_idx, test_idx = fold_splits[fold_idx]

    test_pheno = master_pheno_df.iloc[test_idx].copy()

    common_probes = list(set(master_m_matrix.columns).intersection(set(manifest_df['IlmnID'])))
    
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)
    val_ds = WholeBloodMethylationDataset(master_m_matrix, test_pheno, common_probes, chr_topologies)

    val_loader = DataLoader(
        val_ds, 
        batch_size=16, 
        shuffle=False,
        collate_fn=chromosome_collate_fn,
        num_workers=6,
        pin_memory=True,
        persistent_workers=True,
        prefetch_factor=1
    )

    checkpoint_path = f"{results_path}/fold_{fold_idx + 1}_checkpoint.pt"
    model = ChromosomeParallelGAT(num_node_classes=7, chr_embed_dim=16, cell_prop_dim=6)
    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    mlp_weights = checkpoint['model_state_dict']['classifier.0.weight']
    genome_weights = mlp_weights[:, :-6]  # Isolate the 352 genetic dimensions

    chr_mlp_importances = {}
    for c_idx in range(22):
        start_col = c_idx * 16
        end_col = (c_idx + 1) * 16
        chr_slice = genome_weights[:, start_col:end_col]
        
        # Calculate the mean absolute weight for this specific chromosome
        chr_mlp_importances[c_idx] = torch.mean(torch.abs(chr_slice)).item()

    gat_alpha_accum = {c_idx: [] for c_idx in range(22)}

    with torch.no_grad():
        for batched_chrs, _, _ in val_loader:
            
            for c_idx in range(22):
                batch = batched_chrs[c_idx]
                
                # Forward pass through GAT layers
                emb_func = model.func_embedding(batch.func_type)
                node_feat = torch.cat([batch.x, emb_func], dim=1)
                
                h = torch.relu(model.gat1(node_feat, batch.edge_index))
                h = torch.relu(model.gat2(h, batch.edge_index))
                
                # Extract raw pre-softmax logits from the pooling gate
                gate_logits = model.gate_nn(h).squeeze(-1)
                
                # Apply softmax dynamically based on patient batching
                alpha = torch.softmax(gate_logits, batch.batch)
                
                # Reshape to [patients_in_batch, nodes_in_chromosome]
                num_patients = int(batch.batch.max().item() + 1)
                nodes_per_patient = int(batch.num_nodes / num_patients)
                alpha_matrix = alpha.view(num_patients, nodes_per_patient)
                
                # Average across patients in this batch and store
                gat_alpha_accum[c_idx].append(alpha_matrix.mean(dim=0).cpu().numpy())

    global_probe_ranking = []

    for c_idx in range(22):
        # Average the GAT attention scores across all batches for this chromosome
        mean_gat_alphas = np.mean(np.array(gat_alpha_accum[c_idx]), axis=0)
        
        # Get the MLP weight for this chromosome
        chr_mlp_weight = chr_mlp_importances[c_idx]
        
        # CALCULATE COMPOUND IMPORTANCE
        compound_scores = mean_gat_alphas * chr_mlp_weight
        
        # Load the specific 450k manifest for this chromosome to get Illumina IDs and Gene names
        # Ensure this dataframe is in the exact same order as your graph nodes!
        chr_df = pd.read_parquet(f"processed_chrs/chr_{c_idx + 1}_probes.parquet")
        
        chr_df['chromosome'] = c_idx + 1
        chr_df['gat_alpha'] = mean_gat_alphas
        chr_df['mlp_chr_weight'] = chr_mlp_weight
        chr_df['compound_importance'] = compound_scores
        
        global_probe_ranking.append(chr_df)

    # 5. Concatenate all 22 chromosomes into a single master dataframe
    master_df = pd.concat(global_probe_ranking, ignore_index=True)

    # 6. Sort globally by Compound Importance
    master_df = master_df.sort_values(by='compound_importance', ascending=False).reset_index(drop=True)

    print("--- Top 20 Global Driver Probes for Parkinson's Classification ---")
    print(master_df[['probe_id', 'chromosome', 'gene_symbol', 'gat_alpha', 'mlp_chr_weight', 'compound_importance']].head(20))

    # Save for your thesis / pathway analysis
    master_df.to_csv("global_compound_importance_ranking.csv", index=False)




