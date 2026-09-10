import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from torch.utils.data import DataLoader
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn

def generate_evaluation_plots(m_matrix_df, pheno_df, chr_topologies, cell_cols, results_filepath):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # 1. Recreate the exact splits
    y = pheno_df['Sample_Group'].map({'Control': 0, 'PD': 1}).to_numpy(dtype=np.int64)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Trackers for Mean Curves
    mean_fpr = np.linspace(0, 1, 100)
    train_tprs, train_aucs, train_pr_aucs = [], [], []
    val_tprs, val_aucs, val_pr_aucs = [], [], []
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(pheno_df, y)):
        print(f"Evaluating Fold {fold + 1}...")
        
        # 2. Initialize Train and Validation Loaders
        train_pheno = pheno_df.iloc[train_idx]
        test_pheno = pheno_df.iloc[test_idx]
        
        train_ds = WholeBloodMethylationDataset(m_matrix_df, train_pheno, cell_cols, chr_topologies)
        test_ds = WholeBloodMethylationDataset(m_matrix_df, test_pheno, cell_cols, chr_topologies)
        
        train_loader = DataLoader(train_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        test_loader = DataLoader(test_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        
        # 3. Load the Saved Weights
        model = ChromosomeParallelGAT().to(device)
        fold_bin = torch.load(f"{results_filepath}/gat_fold_{fold + 1}.pt", map_location=device)
        model.load_state_dict(fold_bin['model_state_dict'])
        model.eval()
        
        # 4. Get Training Predictions & Plot Individual Lines
        train_preds, train_truths = [], []
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in train_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                train_preds.extend(torch.sigmoid(logits).cpu().numpy())
                train_truths.extend(batch_y.numpy())
                
        # Train ROC
        fpr, tpr, _ = roc_curve(train_truths, train_preds)
        train_aucs.append(auc(fpr, tpr))
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        train_tprs.append(interp_tpr)
        ax_roc.plot(fpr, tpr, alpha=0.15, color='#ffb347') # Pastel Orange
        
        # Train PR
        precision, recall, _ = precision_recall_curve(train_truths, train_preds)
        train_pr_aucs.append(average_precision_score(train_truths, train_preds))
        ax_pr.plot(recall, precision, alpha=0.15, color='#ffb347')

        # 5. Get Validation Predictions & Plot Individual Lines
        val_preds, val_truths = [], []
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in test_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                val_preds.extend(torch.sigmoid(logits).cpu().numpy())
                val_truths.extend(batch_y.numpy())
                
        # Validation ROC
        fpr, tpr, _ = roc_curve(val_truths, val_preds)
        val_aucs.append(auc(fpr, tpr))
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        val_tprs.append(interp_tpr)
        ax_roc.plot(fpr, tpr, alpha=0.15, color='#779ecb') # Pastel Blue
        
        # Validation PR
        precision, recall, _ = precision_recall_curve(val_truths, val_preds)
        val_pr_aucs.append(average_precision_score(val_truths, val_preds))
        ax_pr.plot(recall, precision, alpha=0.15, color='#779ecb')
        
        # VRAM Protection
        del train_ds, test_ds, train_loader, test_loader, model
        torch.cuda.empty_cache()

    # --- 6. Plot Formatting ---
    # ROC Plot
    mean_train_tpr = np.mean(train_tprs, axis=0)
    mean_train_tpr[-1] = 1.0
    ax_roc.plot(mean_fpr, mean_train_tpr, color='#d9822b', lw=2.5, label=f'Train Mean ROC (AUC = {np.mean(train_aucs):.2f} $\\pm$ {np.std(train_aucs):.2f})')
    
    mean_val_tpr = np.mean(val_tprs, axis=0)
    mean_val_tpr[-1] = 1.0
    ax_roc.plot(mean_fpr, mean_val_tpr, color='#3b719f', lw=2.5, label=f'Val Mean ROC (AUC = {np.mean(val_aucs):.2f} $\\pm$ {np.std(val_aucs):.2f})')
    
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=1.5, color='gray', label='Chance')
    ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Receiver Operating Characteristic (ROC)", xlabel="False Positive Rate", ylabel="True Positive Rate")
    ax_roc.legend(loc="lower right")
    
    # PR Plot (Using proxy artists for clean legends)
    baseline = np.sum(y) / len(y)
    ax_pr.axhline(y=baseline, color='gray', linestyle='--', label=f'Baseline ({baseline:.2f})')
    
    ax_pr.plot([], [], color='#d9822b', lw=2.5, label=f'Train Mean AP = {np.mean(train_pr_aucs):.2f} $\\pm$ {np.std(train_pr_aucs):.2f}')
    ax_pr.plot([], [], color='#3b719f', lw=2.5, label=f'Val Mean AP = {np.mean(val_pr_aucs):.2f} $\\pm$ {np.std(val_pr_aucs):.2f}')
    
    ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Precision-Recall Curve", xlabel="Recall", ylabel="Precision")
    ax_pr.legend(loc="lower left")
    
    plt.tight_layout()
    plt.savefig(f"{results_filepath}/roc_pr_auc_plots.png", dpi=300, bbox_inches="tight")
    plt.close()