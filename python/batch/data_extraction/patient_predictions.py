import pandas as pd
import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold
from torch.utils.data import DataLoader

from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn, build_chromosome_topologies


def export_patient_predictions(m_matrix_df, pheno_df, chr_topologies, cell_cols, results_filepath):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Extracting patient probabilities on: {device}")
    
    y = pheno_df['Sample_Group'].map({'Control': 0, 'PD': 1}).to_numpy(dtype=np.int64)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    patient_records = []
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(pheno_df, y)):
        print(f"Processing Fold {fold + 1}...")
        
        train_pheno = pheno_df.iloc[train_idx]
        test_pheno = pheno_df.iloc[test_idx]
        
        train_ds = WholeBloodMethylationDataset(m_matrix_df, train_pheno, cell_cols, chr_topologies)
        test_ds = WholeBloodMethylationDataset(m_matrix_df, test_pheno, cell_cols, chr_topologies)
        
        # shuffle=False guarantees indices perfectly align with outputs
        train_loader = DataLoader(train_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        test_loader = DataLoader(test_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        
        model = ChromosomeParallelGAT().to(device)
        fold_bin = torch.load(f"{results_filepath}/gat_fold_{fold + 1}.pt", map_location=device)
        model.load_state_dict(fold_bin['model_state_dict'])
        model.eval()
        
        # Extract Training Probabilities
        train_preds = []
        with torch.no_grad():
            for batched_chrs, u_cells, _ in train_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                train_preds.extend(torch.sigmoid(logits).cpu().numpy())
                
        for sn, truth, pred in zip(train_pheno.index, train_pheno['Sample_Group'], train_preds):
            patient_records.append({
                'Sample_Name': sn,
                'Fold': fold + 1,
                'Set': 'Train',
                'True_Label': truth,
                'PD_Probability': pred
            })

        # Extract Validation Probabilities
        val_preds = []
        with torch.no_grad():
            for batched_chrs, u_cells, _ in test_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                val_preds.extend(torch.sigmoid(logits).cpu().numpy())
                
        for sn, truth, pred in zip(test_pheno.index, test_pheno['Sample_Group'], val_preds):
            patient_records.append({
                'Sample_Name': sn,
                'Fold': fold + 1,
                'Set': 'Validation',
                'True_Label': truth,
                'PD_Probability': pred
            })
            
        del train_ds, test_ds, train_loader, test_loader, model
        torch.cuda.empty_cache()

    # Export to CSV
    predictions_df = pd.DataFrame(patient_records)
    export_path = f"{results_filepath}/patient_probability_distributions.csv"
    predictions_df.to_csv(export_path, index=False)
    print(f"Exported complete probability distributions to {export_path}")
    
    return predictions_df