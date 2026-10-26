import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, 
                             average_precision_score, accuracy_score, 
                             precision_score, recall_score, f1_score, confusion_matrix)
from torch.utils.data import DataLoader
from gat import ChromosomeParallelGAT, WholeBloodMethylationDataset, chromosome_collate_fn

def generate_evaluation_plots(m_matrix_df, pheno_df, chr_topologies, cell_cols, results_filepath, panel_labels=('A', 'B', 'C'), file_suffix=""):
    """
    panel_labels: tuple of THREE strings for ROC, PR, and CM panels respectively (e.g., ('A', 'B', 'C'))
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
    
    # Global accumulators for the Aggregated Confusion Matrix
    global_val_preds = []
    global_val_truths = []
    global_val_preds_bin = []
    
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
                
        # Append to global lists for the confusion matrix later
        global_val_preds.extend(val_preds)
        global_val_truths.extend(val_truths)

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
        # Find the optimal threshold for this specific fold using Youden's J on the TRAIN set
        fpr_tr, tpr_tr, thresholds_tr = roc_curve(train_truths, train_preds)
        j_scores_tr = tpr_tr - fpr_tr
        opt_thresh_fold = thresholds_tr[np.argmax(j_scores_tr)]
        
        # Apply the fold-specific threshold instead of a static 0.5
        tr_preds_bin = (np.array(train_preds) >= opt_thresh_fold).astype(int)
        v_preds_bin = (np.array(val_preds) >= opt_thresh_fold).astype(int)
        
        global_val_preds_bin.extend(v_preds_bin)
        
        metrics_records.extend([{
            'Fold': fold + 1, 'Set': 'Train', 'ROC_AUC': tr_roc_auc, 'PR_AUC': tr_pr_auc,
            'Accuracy': accuracy_score(train_truths, tr_preds_bin),
            'Precision': precision_score(train_truths, tr_preds_bin, zero_division=0),
            'Recall': recall_score(train_truths, tr_preds_bin, zero_division=0),
            'F1': f1_score(train_truths, tr_preds_bin, zero_division=0),
            'Opt_Threshold': opt_thresh_fold
        }, {
            'Fold': fold + 1, 'Set': 'Validation', 'ROC_AUC': v_roc_auc, 'PR_AUC': v_pr_auc,
            'Accuracy': accuracy_score(val_truths, v_preds_bin),
            'Precision': precision_score(val_truths, v_preds_bin, zero_division=0),
            'Recall': recall_score(val_truths, v_preds_bin, zero_division=0),
            'F1': f1_score(val_truths, v_preds_bin, zero_division=0),
            'Opt_Threshold': opt_thresh_fold
        }])

        del train_ds, test_ds, train_loader, test_loader, model
        torch.cuda.empty_cache()

    # --- 4. Export Metrics to CSV ---
    metrics_df = pd.DataFrame(metrics_records)
    metrics_df.to_csv(f"{results_filepath}/{prefix}cv_evaluation_metrics.csv", index=False)
    
    # --- 5. Prepare Aggregated Data for Confusion Matrix ---
    print("\nGenerating Aggregated Confusion Matrix from fold-calibrated predictions...")
    y_true_all = np.array(global_val_truths)
    y_preds_all = np.array(global_val_preds_bin)
    
    cm = confusion_matrix(y_true_all, y_preds_all)
    tn, fp, fn, tp = cm.ravel()

    # --- 6. Generate Global 1x3 Mean Plot ---
    # Expanded figsize from 16 to 24 wide to fit all 3 panels perfectly
    fig_mean, (ax_roc, ax_pr, ax_cm) = plt.subplots(1, 3, figsize=(24, 7))
    
    # ---------------- PANEL A: ROC CURVE ----------------
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
    
    # ---------------- PANEL B: PR CURVE ----------------
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

    # ---------------- PANEL C: CONFUSION MATRIX ----------------
    cm_percentages = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    labels = np.asarray([
        f"True Negative\n{tn}\n({cm_percentages[0, 0]:.1%})",
        f"False Positive\n{fp}\n({cm_percentages[0, 1]:.1%})",
        f"False Negative\n{fn}\n({cm_percentages[1, 0]:.1%})",
        f"True Positive\n{tp}\n({cm_percentages[1, 1]:.1%})"
    ]).reshape(2, 2)

    # Note the ax=ax_cm argument here routes seaborn directly into matplotlib's grid
    sns.heatmap(
        cm, annot=labels, fmt="", cmap='Blues', cbar=False, square=True,
        annot_kws={"size": 15, "weight": "bold"}, linewidths=2, linecolor='black', ax=ax_cm
    )

    ax_cm.set_title(panel_labels[2], loc="left", fontsize=22, fontweight="bold")
    ax_cm.set_title("Aggregated CM (Fold-Calibrated)", loc="center", fontsize=18)
    ax_cm.set_xlabel('Predicted Diagnosis', fontsize=16, fontweight='bold', labelpad=15)
    ax_cm.set_ylabel('True Clinical Diagnosis', fontsize=16, fontweight='bold', labelpad=15)
    ax_cm.set_xticklabels(['Control (0)', 'PD (1)'], fontsize=14)
    ax_cm.set_yticklabels(['Control (0)', 'PD (1)'], fontsize=14, rotation=0, va='center')
    
    # ---------------- FINAL EXPORT ----------------
    plt.tight_layout()
    fig_mean.savefig(f"{results_filepath}/{prefix}evaluation_unified_panels.pdf", format="pdf", bbox_inches="tight")
    fig_mean.savefig(f"{results_filepath}/{prefix}evaluation_unified_panels.png", dpi=300, bbox_inches="tight")
    plt.close(fig_mean)
    print(f"Saved unified 3-panel figure to {prefix}evaluation_unified_panels.pdf")