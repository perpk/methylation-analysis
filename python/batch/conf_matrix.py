import torch
import numpy as np
from sklearn.metrics import precision_recall_curve, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import os

from gat import ChromosomeParallelGAT

# Instantiate the model architecture once
model = ChromosomeParallelGAT(num_node_classes=7, chr_embed_dim=16, cell_prop_dim=6)
model.eval() # Critical: disable dropout for inference

cohorts = {
    "PEG1": "/workspace/data/peg1",
    "SGPD": "/workspace/data/sgpd",
    "PPMI": "/workspace/data/ppmi"
}

output_dir = "/workspace/results/confusion_matrices/"
os.makedirs(output_dir, exist_ok=True)

# Loop through each cohort independently
for cohort_name, cohort_path in cohorts.items():
    print(f"\n{'='*50}\nEvaluating Cohort: {cohort_name}\n{'='*50}")
    
    # CRITICAL FIX 1: Reset accumulators FOR EACH COHORT
    all_oof_probs = []
    all_oof_targets = []
    
    for fold in range(1, 6):
        checkpoint_path = f"{cohort_path}/gat_fold_{fold}.pt"
        try:
            checkpoint = torch.load(checkpoint_path, map_location='cpu')
            model.load_state_dict(checkpoint['model_state_dict'])
        except FileNotFoundError:
            print(f"  -> Skipping fold {fold}: Checkpoint not found.")
            continue

        # ========================================================
        # TODO: Insert your validation loader logic here
        # val_loader = get_fold_val_loader(cohort_name, fold)
        # ========================================================
        
        with torch.no_grad():
            for batched_chrs, u_cells, batch_y in val_loader:
                # Forward pass
                logits = model(batched_chrs, u_cells)
                
                # Convert logits to probabilities
                probs = torch.sigmoid(logits).squeeze()
                
                # Ensure probs is iterable even if batch_size=1
                if probs.dim() == 0:
                    probs = probs.unsqueeze(0)
                    
                all_oof_probs.extend(probs.cpu().numpy())
                all_oof_targets.extend(batch_y.cpu().numpy())

    # Skip to next cohort if no data was processed
    if len(all_oof_targets) == 0:
        print(f"No validation data processed for {cohort_name}. Skipping...")
        continue

    # Convert pooled lists to numpy arrays
    y_true = np.array(all_oof_targets)
    y_probs = np.array(all_oof_probs)

    # Calculate global metrics
    auc_score = roc_auc_score(y_true, y_probs)
    precision, recall, thresholds = precision_recall_curve(y_true, y_probs)

    # Find the threshold that maximizes the F1 Score
    f1_scores = (2 * precision * recall) / (precision + recall + 1e-8)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx]
    best_f1 = f1_scores[optimal_idx]

    print(f"--- Threshold Optimization ---")
    print(f"Optimal Probability Threshold : {optimal_threshold:.4f}")
    print(f"Maximum F1 Score              : {best_f1:.4f}")
    print(f"Area Under Curve (AUROC)      : {auc_score:.4f}\n")

    # Apply the optimal threshold to get binary predictions
    y_preds = (y_probs >= optimal_threshold).astype(int)

    # Generate the Aggregated Confusion Matrix
    cm = confusion_matrix(y_true, y_preds)
    tn, fp, fn, tp = cm.ravel()
    
    # Calculate clinical metrics
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

    print(f"--- Clinical Metrics ({cohort_name}) ---")
    print(f"Sensitivity (Recall) : {sensitivity:.1%}")
    print(f"Specificity          : {specificity:.1%}\n")

    # ==========================================
    # Publication-Ready Plotting
    # ==========================================
    plt.figure(figsize=(7, 6))

    # Calculate row-wise percentages for the plot labels
    cm_percentages = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    # Create custom labels showing count AND percentage
    labels = np.asarray([
        f"True Negative\n{tn}\n({cm_percentages[0, 0]:.1%})",
        f"False Positive\n{fp}\n({cm_percentages[0, 1]:.1%})",
        f"False Negative\n{fn}\n({cm_percentages[1, 0]:.1%})",
        f"True Positive\n{tp}\n({cm_percentages[1, 1]:.1%})"
    ]).reshape(2, 2)

    # Create a clean, light-themed heatmap
    ax = sns.heatmap(
        cm, 
        annot=labels, 
        fmt="", # Set to empty string since our labels are pre-formatted strings
        cmap='Blues', 
        cbar=False, 
        square=True,
        annot_kws={"size": 13, "weight": "bold"},
        linewidths=2,
        linecolor='black' # Black borders look sharper in print than white
    )

    # Customize fonts and layout for a clean, academic look
    ax.set_title(f'Aggregated Cross-Validation Performance\nCohort: {cohort_name}', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Predicted Diagnosis (GATv2)', fontsize=14, fontweight='bold', labelpad=15)
    ax.set_ylabel('True Clinical Diagnosis', fontsize=14, fontweight='bold', labelpad=15)

    class_names = ['Control (0)', 'Parkinson\'s (1)']
    ax.set_xticklabels(class_names, fontsize=12)
    ax.set_yticklabels(class_names, fontsize=12, rotation=0, va='center')

    plt.tight_layout()
    
    # Export as both PNG and PDF
    png_path = os.path.join(output_dir, f"Fig_CM_{cohort_name}.png")
    pdf_path = os.path.join(output_dir, f"Fig_CM_{cohort_name}.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    
    print(f"Saved plots to:\n -> {png_path}\n -> {pdf_path}")
    
    # CRITICAL FIX 2: Close the figure to prevent overlapping plots
    plt.close()