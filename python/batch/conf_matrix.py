import torch
import numpy as np
from sklearn.metrics import precision_recall_curve, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

from gat import ChromosomeParallelGAT

all_oof_probs = []
all_oof_targets = []

# Instantiate the model architecture once
model = ChromosomeParallelGAT(num_node_classes=7, chr_embed_dim=16, cell_prop_dim=6)
model.eval() # Critical: disable dropout for inference

cohorts = {
    "peg1" : "/workspace/results/peg1",
    "sgpd" : "/workspace/results/sgpd",
    "ppmi" : "/workspace/results/ppmi"
}

# 2. Loop through all 5 folds
for cohort_name, cohort_path in cohorts.items():
    for fold in range(1, 6):
        checkpoint_path = f"{cohort_path}/gat_fold_{fold}.pt"
        try:
            checkpoint = torch.load(checkpoint_path, map_location='cpu')
            model.load_state_dict(checkpoint['model_state_dict'])
        except FileNotFoundError:
            print(f"Skipping fold {fold}: Checkpoint not found.")
            continue

        # ========================================================
        # Insert your logic here to get the validation loader for the current fold.
        # For example: val_loader = get_fold_val_loader(fold)
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

    # Convert pooled lists to numpy arrays
    y_true = np.array(all_oof_targets)
    y_probs = np.array(all_oof_probs)

    # 3. Calculate the Precision-Recall curve
    precision, recall, thresholds = precision_recall_curve(y_true, y_probs)

    # 4. Find the threshold that maximizes the F1 Score
    # F1 = 2 * (Precision * Recall) / (Precision + Recall)
    # We add a tiny epsilon to the denominator to prevent division by zero
    f1_scores = (2 * precision * recall) / (precision + recall + 1e-8)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx]
    best_f1 = f1_scores[optimal_idx]

    print(f"--- Global Threshold Optimization ---")
    print(f"Optimal Probability Threshold : {optimal_threshold:.4f}")
    print(f"Maximum F1 Score              : {best_f1:.4f}\n")

    # 5. Apply the optimal threshold to get binary predictions
    y_preds = (y_probs >= optimal_threshold).astype(int)

    # 6. Generate the Aggregated Confusion Matrix
    cm = confusion_matrix(y_true, y_preds)
    tn, fp, fn, tp = cm.ravel()

    print(f"--- Aggregated Cross-Validation Confusion Matrix for {cohort_name} ---")
    print(cm, "\n")
    print(f"True Negatives  (Correct Controls) : {tn}")
    print(f"False Positives (False Alarms)     : {fp}")
    print(f"False Negatives (Missed PD Cases)  : {fn}")
    print(f"True Positives  (Correct PD Cases) : {tp}")

    plt.figure(figsize=(8, 6))

    # Define labels for the axes
    class_names = ['Control', 'PD Case']

    # Create a clean, light-themed heatmap
    # cmap="Blues" provides a clean gradient, 'd' ensures integers are printed
    ax = sns.heatmap(
        cm, 
        annot=True, 
        fmt='d', 
        cmap='Blues', 
        cbar=False, 
        square=True,
        annot_kws={"size": 16},
        linewidths=1,
        linecolor='white'
    )

    # Customize fonts and layout for a clean, academic look
    ax.set_title(f'Aggregated Validation Confusion Matrix for {cohort_name}\n(Optimized Threshold)', fontsize=16, pad=20)
    ax.set_xlabel('Predicted Label', fontsize=14, labelpad=15)
    ax.set_ylabel('True Label', fontsize=14, labelpad=15)

    ax.set_xticklabels(class_names, fontsize=12)
    ax.set_yticklabels(class_names, fontsize=12, rotation=0)

    # Adjust layout and save the high-resolution image
    plt.tight_layout()
    plt.savefig(f"/workspace/results/{cohort_name}/aggregated_confusion_matrix_{cohort_name}.png", dpi=300, bbox_inches='tight')