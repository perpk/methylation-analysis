import torch
import pandas as pd
from torch.utils.data import DataLoader
# Import your updated class and dataset functions
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies
import sys

def main(m_matrix_filepath, manifest_filepath, results_filepath):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    fold_num = 1 # setting this for now to 1 - to be changed later

    print(f"Extracting attention from Fold {fold_num} on {device}")
    m_matrix_df = pd.read_parquet(m_matrix_filepath)
    manifest_df = pd.read_parquet(manifest_filepath)
    pheno_df = pd.read_parquet(f"{results_filepath}/pheno_data.parquet")
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']
    common_probes = list(set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID'])))
    print("Building topologies...")
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)
    infer_ds = WholeBloodMethylationDataset(m_matrix_df, pheno_df, cell_cols, chr_topologies)
    infer_loader = DataLoader(infer_ds, batch_size=1, shuffle=False, collate_fn=chromosome_collate_fn)
    model = ChromosomeParallelGAT().to(device)
    checkpoint = torch.load(pt_filepath, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval() # CRITICAL: Sets model to inference mode
    print(f"Model loaded. Validation ROC AUC for this fold was: {checkpoint['val_roc_auc']:.4f}")

    with torch.no_grad():
        batched_chrs, u_cells, batch_y = next(iter(infer_loader))
        
        batched_chrs = [c.to(device) for c in batched_chrs]
        u_cells = u_cells.to(device)
        
        logits, attention_dict = model(batched_chrs, u_cells, extract_attention=True)
        
    edges, weights = attention_dict['chr1']['layer_1']
    
    edges_cpu = edges.cpu().numpy()
    weights_cpu = weights.mean(dim=1).cpu().numpy()
    
    df_attention = pd.DataFrame({
        'Source_Node_Idx': edges_cpu[0],
        'Target_Node_Idx': edges_cpu[1],
        'Attention_Weight': weights_cpu
    })
    
    df_attention = df_attention.sort_values(by='Attention_Weight', ascending=False)
    
    out_file = f"{results_filepath}/fold_{fold_num}_chr1_attention.csv"
    df_attention.to_csv(out_file, index=False)
    print(f"Saved top co-methylation weights to {out_file}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
