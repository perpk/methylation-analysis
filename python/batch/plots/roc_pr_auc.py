import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, 
                             average_precision_score, accuracy_score, 
                             precision_score, recall_score, f1_score)
from torch.utils.data import DataLoader
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn

def generate_evaluation_plots(m_matrix_df, pheno_df, chr_topologies, cell_cols, results_filepath, panel_labels=('A', 'B'), file_suffix=""):
    """
    panel_labels: tuple of two strings for the ROC and PR panels respectively (e.g., ('A', 'B') or ('C', 'D'))
    file_suffix: string to append to saved filenames so multiple sets don't overwrite each other
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    y = pheno_df['Sample_Group'].map({'Control': 0, 'PD': 1}).to_numpy(dtype=np.int64)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    mean_fpr = np.linspace(0, 1, 100)
    mean_recall = np.linspace(0, 1, 100)
    
    train_tprs, train_aucs, train_pr_aucs, train_precs = [], [], [], []
    val_tprs, val_aucs, val_pr_aucs, val_precs = [], [], [], []
    metrics_records = []
    
    # Prefix formatting for files
    prefix = f"{file_suffix}_" if file_suffix else ""
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(pheno_df, y)):
        print(f"Evaluating Fold {fold + 1}...")
        
        train_pheno = pheno_df.iloc[train_idx]
        test_pheno = pheno_df.iloc[test_idx]
        
        train_ds = WholeBloodMethylationDataset(m_matrix_df, train_pheno, cell_cols, chr_topologies)
        test_ds = WholeBloodMethylationDataset(m_matrix_df, test_pheno, cell_cols, chr_topologies)
        
        train_loader = DataLoader(train_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        test_loader = DataLoader(test_ds, batch_size=16, shuffle=False, collate_fn=chromosome_collate_fn)
        
        model = ChromosomeParallelGAT().to(device)
        fold_bin = torch.load(f"/workspace/data/{file_suffix}/gat_fold_{fold + 1}.pt", map_location=device)
        model.load_state_dict(fold_bin['model_state_dict'])
        model.eval()
        
        # --- 1. Get Predictions ---
        train_preds, train_truths = [], []
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in train_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                train_preds.extend(torch.sigmoid(logits).cpu().numpy())
                train_truths.extend(batch_y.numpy())
                
        val_preds, val_truths = [], []
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in test_loader:
                batched_chrs = [c.to(device) for c in batched_chrs]
                logits = model(batched_chrs, u_cells.to(device)).squeeze(1)
                val_preds.extend(torch.sigmoid(logits).cpu().numpy())
                val_truths.extend(batch_y.numpy())

        # --- 2. Calculate Soft Metrics (AUCs) ---
        tr_fpr, tr_tpr, _ = roc_curve(train_truths, train_preds)
        tr_roc_auc = auc(tr_fpr, tr_tpr)
        tr_prec, tr_rec, _ = precision_recall_curve(train_truths, train_preds)
        tr_pr_auc = average_precision_score(train_truths, train_preds)
        
        train_aucs.append(tr_roc_auc)
        train_pr_aucs.append(tr_pr_auc)
        interp_tpr = np.interp(mean_fpr, tr_fpr, tr_tpr)
        interp_tpr[0] = 0.0
        train_tprs.append(interp_tpr)
        interp_prec = np.interp(mean_recall, tr_rec[::-1], tr_prec[::-1])
        train_precs.append(interp_prec)
        
        v_fpr, v_tpr, _ = roc_curve(val_truths, val_preds)
        v_roc_auc = auc(v_fpr, v_tpr)
        v_prec, v_rec, _ = precision_recall_curve(val_truths, val_preds)
        v_pr_auc = average_precision_score(val_truths, val_preds)
        
        val_aucs.append(v_roc_auc)
        val_pr_aucs.append(v_pr_auc)
        interp_tpr = np.interp(mean_fpr, v_fpr, v_tpr)
        interp_tpr[0] = 0.0
        val_tprs.append(interp_tpr)
        interp_prec = np.interp(mean_recall, v_rec[::-1], v_prec[::-1])
        val_precs.append(interp_prec)

        # --- 3. Calculate Hard Metrics ---
        tr_preds_bin = (np.array(train_preds) >= 0.5).astype(int)
        v_preds_bin = (np.array(val_preds) >= 0.5).astype(int)
        
        metrics_records.extend([{
            'Fold': fold + 1, 'Set': 'Train', 'ROC_AUC': tr_roc_auc, 'PR_AUC': tr_pr_auc,
            'Accuracy': accuracy_score(train_truths, tr_preds_bin),
            'Precision': precision_score(train_truths, tr_preds_bin, zero_division=0),
            'Recall': recall_score(train_truths, tr_preds_bin, zero_division=0),
            'F1': f1_score(train_truths, tr_preds_bin, zero_division=0)
        }, {
            'Fold': fold + 1, 'Set': 'Validation', 'ROC_AUC': v_roc_auc, 'PR_AUC': v_pr_auc,
            'Accuracy': accuracy_score(val_truths, v_preds_bin),
            'Precision': precision_score(val_truths, v_preds_bin, zero_division=0),
            'Recall': recall_score(val_truths, v_preds_bin, zero_division=0),
            'F1': f1_score(val_truths, v_preds_bin, zero_division=0)
        }])

        # --- 4. Generate Per-Fold Plot ---
        fig_fold, (ax_f_roc, ax_f_pr) = plt.subplots(1, 2, figsize=(16, 7))
        
        ax_f_roc.plot(tr_fpr, tr_tpr, color='#ffb347', lw=2.5, label=f'Train ROC (AUC = {tr_roc_auc:.4f})')
        ax_f_roc.plot(v_fpr, v_tpr, color='#779ecb', lw=2.5, label=f'Val ROC (AUC = {v_roc_auc:.4f})')
        ax_f_roc.plot([0, 1], [0, 1], linestyle='--', lw=2, color='gray', label='Chance')
        
        ax_f_roc.set_xlim([-0.05, 1.05])
        ax_f_roc.set_ylim([-0.05, 1.05])
        ax_f_roc.set_xlabel("False Positive Rate", fontsize=16)
        ax_f_roc.set_ylabel("True Positive Rate", fontsize=16)
        ax_f_roc.tick_params(axis='both', labelsize=14)
        ax_f_roc.set_title(panel_labels[0], loc="left", fontsize=22, fontweight="bold")
        ax_f_roc.set_title(f"Fold {fold + 1} - ROC", loc="center", fontsize=18)
        ax_f_roc.legend(loc="lower right", fontsize=14)
        ax_f_roc.grid(True, linestyle='--', alpha=0.5)
        
        baseline_fold = sum(val_truths) / len(val_truths)
        ax_f_pr.plot(tr_rec, tr_prec, color='#ffb347', lw=2.5, label=f'Train PR (AP = {tr_pr_auc:.4f})')
        ax_f_pr.plot(v_rec, v_prec, color='#779ecb', lw=2.5, label=f'Val PR (AP = {v_pr_auc:.4f})')
        ax_f_pr.axhline(y=baseline_fold, color='gray', linestyle='--', lw=2, label='Val Baseline')
        
        ax_f_pr.set_xlim([-0.05, 1.05])
        ax_f_pr.set_ylim([-0.05, 1.05])
        ax_f_pr.set_xlabel("Recall", fontsize=16)
        ax_f_pr.set_ylabel("Precision", fontsize=16)
        ax_f_pr.tick_params(axis='both', labelsize=14)
        ax_f_pr.set_title(panel_labels[1], loc="left", fontsize=22, fontweight="bold")
        ax_f_pr.set_title(f"Fold {fold + 1} - PR", loc="center", fontsize=18)
        ax_f_pr.legend(loc="lower left", fontsize=14)
        ax_f_pr.grid(True, linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        fig_fold.savefig(f"{results_filepath}/{prefix}fold_{fold + 1}_evaluation_curves.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig_fold)

        del train_ds, test_ds, train_loader, test_loader, model
        torch.cuda.empty_cache()

    # --- 5. Export Metrics to CSV ---
    metrics_df = pd.DataFrame(metrics_records)
    metrics_df.to_csv(f"{results_filepath}/{prefix}cv_evaluation_metrics.csv", index=False)
    print(f"\nSaved detailed evaluation metrics to {prefix}cv_evaluation_metrics.csv")

    # --- 6. Generate Global Mean Plot ---
    fig_mean, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(16, 7))
    
    mean_train_tpr = np.mean(train_tprs, axis=0)
    mean_train_tpr[-1] = 1.0
    ax_roc.plot(mean_fpr, mean_train_tpr, color='#d9822b', lw=3, label=f'Train Mean ROC (AUC = {np.mean(train_aucs):.2f} $\\pm$ {np.std(train_aucs):.2f})')
    
    mean_val_tpr = np.mean(val_tprs, axis=0)
    mean_val_tpr[-1] = 1.0
    ax_roc.plot(mean_fpr, mean_val_tpr, color='#3b719f', lw=3, label=f'Val Mean ROC (AUC = {np.mean(val_aucs):.2f} $\\pm$ {np.std(val_aucs):.2f})')
    
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=2, color='gray', label='Chance')
    ax_roc.set_xlim([-0.05, 1.05])
    ax_roc.set_ylim([-0.05, 1.05])
    ax_roc.set_xlabel("False Positive Rate", fontsize=16)
    ax_roc.set_ylabel("True Positive Rate", fontsize=16)
    ax_roc.tick_params(axis='both', labelsize=14)
    ax_roc.set_title(panel_labels[0], loc="left", fontsize=22, fontweight="bold")
    ax_roc.set_title(f"Mean ROC - {file_suffix.upper()}" if file_suffix else "Mean ROC", loc="center", fontsize=18)
    ax_roc.legend(loc="lower right", fontsize=14)
    ax_roc.grid(True, linestyle='--', alpha=0.5)
    
    baseline_global = np.sum(y) / len(y)
    ax_pr.axhline(y=baseline_global, color='gray', linestyle='--', lw=2, label=f'Overall Baseline ({baseline_global:.2f})')
    
    mean_train_prec = np.mean(train_precs, axis=0)
    ax_pr.plot(mean_recall, mean_train_prec, color='#d9822b', lw=3, label=f'Train Mean PR (AP = {np.mean(train_pr_aucs):.2f} $\\pm$ {np.std(train_pr_aucs):.2f})')
    
    mean_val_prec = np.mean(val_precs, axis=0)
    ax_pr.plot(mean_recall, mean_val_prec, color='#3b719f', lw=3, label=f'Val Mean PR (AP = {np.mean(val_pr_aucs):.2f} $\\pm$ {np.std(val_pr_aucs):.2f})')
    
    ax_pr.set_xlim([-0.05, 1.05])
    ax_pr.set_ylim([-0.05, 1.05])
    ax_pr.set_xlabel("Recall", fontsize=16)
    ax_pr.set_ylabel("Precision", fontsize=16)
    ax_pr.tick_params(axis='both', labelsize=14)
    ax_pr.set_title(panel_labels[1], loc="left", fontsize=22, fontweight="bold")
    ax_pr.set_title(f"Mean PR - {file_suffix.upper()}" if file_suffix else "Mean Precision-Recall", loc="center", fontsize=18)
    ax_pr.legend(loc="lower left", fontsize=14)
    ax_pr.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    fig_mean.savefig(f"{results_filepath}/{prefix}mean_evaluation_curves.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig_mean)