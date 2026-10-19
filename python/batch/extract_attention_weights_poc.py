import torch
import pandas as pd
from torch.utils.data import DataLoader
# Import your updated class and dataset functions
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies

def extract_attention_weights(fold_num, m_matrix_filepath, manifest_filepath, results_filepath, pt_filepath):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Extracting attention from Fold {fold_num} on {device}")

    # 1. Load your data (Simplified for inference - you just need the test patients)
    m_matrix_df = pd.read_parquet(m_matrix_filepath)
    manifest_df = pd.read_parquet(manifest_filepath)
    pheno_df = pd.read_parquet(f"{results_filepath}/torchgat2/pheno_data.parquet")
    
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']
    common_probes = list(set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID'])))
    
    # Optional: You can filter pheno_df here to only extract attention from Parkinson's cases
    # pd_only_pheno = pheno_df[pheno_df['Sample_Group'] == 'PD']
    
    print("Building topologies...")
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)

    # 2. Instantiate Dataset and a DataLoader (batch_size=1 is easier to interpret)
    infer_ds = WholeBloodMethylationDataset(m_matrix_df, pheno_df, cell_cols, chr_topologies)
    infer_loader = DataLoader(infer_ds, batch_size=1, shuffle=False, collate_fn=chromosome_collate_fn)

    # 3. Load the pre-trained model from the .pt file
    model = ChromosomeParallelGAT().to(device)
    checkpoint = torch.load(pt_filepath, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval() # CRITICAL: Sets model to inference mode

    print(f"Model loaded. Validation ROC AUC for this fold was: {checkpoint['val_roc_auc']:.4f}")

    # 4. Extract Attention for a single batch
    with torch.no_grad():
        # Just grab the first batch (e.g., a single patient) for extraction
        batched_chrs, u_cells, batch_y = next(iter(infer_loader))
        
        batched_chrs = [c.to(device) for c in batched_chrs]
        u_cells = u_cells.to(device)
        
        # Trigger the modified forward pass
        logits, attention_dict = model(batched_chrs, u_cells, extract_attention=True)
        
    # 5. Parse the attention dictionary for Chromosome 1 (Example)
    # Extracting Layer 1 (Direct neighbor co-methylation)
    edges, weights = attention_dict['chr1']['layer_1']
    
    # edges is shape [2, num_edges]. weights is shape [num_edges, num_heads]
    edges_cpu = edges.cpu().numpy()
    weights_cpu = weights.mean(dim=1).cpu().numpy() # Average the attention heads
    
    # Create a DataFrame mapping the source node to the target node
    df_attention = pd.DataFrame({
        'Source_Node_Idx': edges_cpu[0],
        'Target_Node_Idx': edges_cpu[1],
        'Attention_Weight': weights_cpu
    })
    
    cpg_list = chr_topologies[1]['probes']
    
    # Map the integers to the strings
    df_attention['Source_CpG'] = df_attention['Source_Node_Idx'].map(lambda idx: cpg_list[idx])
    df_attention['Target_CpG'] = df_attention['Target_Node_Idx'].map(lambda idx: cpg_list[idx])
    
    # Reorder the dataframe so the CpG names are the first columns you see
    df_attention = df_attention[
        ['Source_CpG', 'Target_CpG', 'Attention_Weight', 'Source_Node_Idx', 'Target_Node_Idx']
    ]

    # Sort by the strongest co-methylation links
    df_attention = df_attention.sort_values(by='Attention_Weight', ascending=False)
    
    out_file = f"{results_filepath}/fold_{fold_num}_chr1_attention.csv"
    df_attention.to_csv(out_file, index=False)
    print(f"Saved top co-methylation weights to {out_file}")


def main():

    extract_attention_weights(
        fold_num=1,
        m_matrix_filepath="/workspace/data/GSE111629_data_test.parquet",
        manifest_filepath="/workspace/data/infinium450k_manifest.parquet",
        results_filepath="/workspace/results/peg1",
        pt_filepath="/workspace/results/peg1/torchgat2/gat_fold_1.pt"
    )

if __name__ == "__main__":
    main()
