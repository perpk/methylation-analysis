import gc
import os
import glob
import pandas as pd
import numpy as np
import torch
from torch.utils.data import DataLoader
import sys

from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies

def extract_fold_attention(fold_num, pheno_df, m_matrix_df, manifest_df, pt_filepath, chr_topologies):
    """Extracts and averages attention weights for ALL Parkinson's patients in a single fold."""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"  -> Extracting Fold {fold_num} on {device}...")
    
    # 1. Filter for Parkinson's cases only
    pheno_pd = pheno_df[pheno_df['Sample_Group'] == 'PD'].copy()
    
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']
    infer_ds = WholeBloodMethylationDataset(m_matrix_df, pheno_pd, cell_cols, chr_topologies)
    
    # batch_size=16 is safe here because we use the Modulo trick below
    infer_loader = DataLoader(infer_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn, num_workers=4)

    model = ChromosomeParallelGAT().to(device)
    checkpoint = torch.load(pt_filepath, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval() 

    all_batches_list = []

    # 2. Batch Loop (Process all PD patients)
    with torch.no_grad():
        for batched_chrs, u_cells, batch_y in infer_loader:
            batched_chrs = [c.to(device) for c in batched_chrs]
            u_cells = u_cells.to(device)
            
            logits, attention_dict = model(batched_chrs, u_cells, extract_attention=True)
            
            for chr_num in range(1, 23): 
                chr_key = f'chr{chr_num}'
                edges, weights = attention_dict[chr_key]['layer_1']
                
                edges_cpu = edges.cpu().numpy()
                weights_cpu = weights.mean(dim=1).cpu().numpy() 
                
                # The Modulo Trick: Revert stitched batch indices back to 0-to-N
                num_nodes = len(chr_topologies[chr_num]['probes']) # Ensure 'probes' matches your dict key
                src_idx = edges_cpu[0] % num_nodes
                tgt_idx = edges_cpu[1] % num_nodes
                
                batch_df = pd.DataFrame({
                    'Chromosome': chr_num,
                    'Source_Node_Idx': src_idx,
                    'Target_Node_Idx': tgt_idx,
                    'Attention_Weight': weights_cpu
                })
                collapsed_batch = batch_df.groupby(
                    ['Chromosome', 'Source_Node_Idx', 'Target_Node_Idx'], as_index=False
                )['Attention_Weight'].sum()
                
                all_batches_list.append(collapsed_batch)
                
    # 3. Aggregate all batches to get the FOLD's representative average
    df_all_batches = pd.concat(all_batches_list, ignore_index=True)
    fold_df = df_all_batches.groupby(
        ['Chromosome', 'Source_Node_Idx', 'Target_Node_Idx'], as_index=False
    )['Attention_Weight'].sum()
    fold_df['Attention_Weight'] = fold_df['Attention_Weight'] / len(pheno_pd)    
    
    # 4. Map Integer Indices to string CpG names once at the end (saves RAM)
    def map_cpg(row, col_name):
        # chr_key = f"chr{int(row['Chromosome'])}"
        chr_index = int(row['Chromosome'])
        return chr_topologies[chr_index]['probes'][int(row[col_name])]

    fold_df['Source_CpG'] = fold_df.apply(lambda r: map_cpg(r, 'Source_Node_Idx'), axis=1)
    fold_df['Target_CpG'] = fold_df.apply(lambda r: map_cpg(r, 'Target_Node_Idx'), axis=1)
    
    # Keep only the necessary columns
    return fold_df[['Chromosome', 'Source_CpG', 'Target_CpG', 'Attention_Weight']]


def process_cohort_consensus(cohort_name, m_matrix, pheno, pt_dir, manifest, results_filepath):
    """Calculates Mean, Std, and SNR across all 5 folds for a single cohort."""
    print(f"\n=== Processing Cohort: {cohort_name} ===")
    
    common_probes = list(set(m_matrix.columns).intersection(set(manifest['IlmnID'])))
    print(f"Building topologies for {cohort_name}...")
    chr_topologies = build_chromosome_topologies(manifest, common_probes)
    
    # REPLACED GLOB.GLOB: Deterministically generate the 5 expected paths
    pt_files = []
    for fold in range(1, 6):
        # Adjust the prefix if your files are named differently (e.g., 'gat_fold_{fold}.pt')
        expected_filename = f"gat_fold_{fold}.pt" 
        full_path = os.path.join(pt_dir, expected_filename)
        
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"CRITICAL ERROR: Cannot find model at {full_path}")
            
        pt_files.append(full_path)

    fold_dfs = []
    for i, pt_file in enumerate(pt_files):
        df_fold = extract_fold_attention(i + 1, pheno, m_matrix, manifest, pt_file, chr_topologies)
        fold_dfs.append(df_fold)
        
    # Combine all 5 folds
    df_all_folds = pd.concat(fold_dfs, ignore_index=True)
    
    # Calculate Mean, Standard Deviation, and SNR for each edge
    print(f"Calculating {cohort_name} Consensus (Mean, Std, SNR)...")
    cohort_consensus = df_all_folds.groupby(['Chromosome', 'Source_CpG', 'Target_CpG']).agg(
        Mean_Weight=('Attention_Weight', 'mean'),
        Std_Weight=('Attention_Weight', 'std')
    ).reset_index()
    
    # Handle division by zero for edges that have perfectly identical weights (rare, but causes NaNs)
    cohort_consensus['Std_Weight'] = cohort_consensus['Std_Weight'].fillna(1e-8).replace(0, 1e-8)
    cohort_consensus['SNR'] = cohort_consensus['Mean_Weight'] / cohort_consensus['Std_Weight']
    
    cohort_consensus = cohort_consensus.sort_values('SNR', ascending=False)
    
    # Save the Cohort-specific CSV
    out_path = os.path.join(results_filepath, f"{cohort_name}_edge_consensus.csv")
    cohort_consensus.to_csv(out_path, index=False)
    print(f"Saved {cohort_name} consensus to {out_path}")
    
    return cohort_consensus

def main(peg1_m_matrix_filepath, sgpd_m_matrix_filepath, ppmi_m_matrix_filepath,
         peg1_pt_dir, sgpd_pt_dir, ppmi_pt_dir,
         peg1_pheno_filepath, sgpd_pheno_filepath, ppmi_pheno_filepath,
         manifest_filepath, results_filepath):
    
    print("Reading data...")
    manifest = pd.read_parquet(manifest_filepath)
    
    cohort_paths = {
        # "PEG1": (peg1_m_matrix_filepath, peg1_pheno_filepath, peg1_pt_dir),
        "SGPD": (sgpd_m_matrix_filepath, sgpd_pheno_filepath, sgpd_pt_dir),
        "PPMI": (ppmi_m_matrix_filepath, ppmi_pheno_filepath, ppmi_pt_dir)
    }

    cohort_consensus_dfs = []
    
    # 1. Calculate Intra-Cohort Consensus
    for cohort_name, (m_path, pheno_path, pt_dir) in cohort_paths.items():        
        print(f"\n--- Loading {cohort_name} into RAM ---")
        m_matrix = pd.read_parquet(m_path)
        pheno = pd.read_parquet(pheno_path)
        
        df_consensus = process_cohort_consensus(cohort_name, m_matrix, pheno, pt_dir, manifest, results_filepath)
        df_consensus.rename(columns={
            'Mean_Weight': f'{cohort_name}_Mean',
            'Std_Weight': f'{cohort_name}_Std',
            'SNR': f'{cohort_name}_SNR'
        }, inplace=True)
        cohort_consensus_dfs.append(df_consensus)

    del m_matrix
    del pheno
    gc.collect()

    # 2. Calculate Master Inter-Cohort Consensus
    print("\n=== Calculating Master 3-Cohort Consensus ===")
    
    # Merge the three dataframes strictly on the edges they share
    master_df = cohort_consensus_dfs[0]
    for df in cohort_consensus_dfs[1:]:
        master_df = pd.merge(master_df, df, on=['Chromosome', 'Source_CpG', 'Target_CpG'], how='inner')
    
    # Calculate the Meta-Mean, Meta-Std, and Meta-SNR across the 3 Cohort Means
    master_df['Master_Mean'] = master_df[['PEG1_Mean', 'SGPD_Mean', 'PPMI_Mean']].mean(axis=1)
    master_df['Master_Std'] = master_df[['PEG1_Mean', 'SGPD_Mean', 'PPMI_Mean']].std(axis=1)
    
    master_df['Master_Std'] = master_df['Master_Std'].fillna(1e-8).replace(0, 1e-8)
    master_df['Master_SNR'] = master_df['Master_Mean'] / master_df['Master_Std']
    
    master_df = master_df.sort_values('Master_SNR', ascending=False)
    
    out_path = os.path.join(results_filepath, "MASTER_3Cohort_Edge_Consensus.csv")
    master_df.to_csv(out_path, index=False)
    print(f"Success! Master consensus saved to {out_path}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3],
         sys.argv[4], sys.argv[5], sys.argv[6],
         sys.argv[7], sys.argv[8], sys.argv[9],
         sys.argv[10], sys.argv[11])