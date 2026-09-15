import torch

# 1. Load the checkpoint (mapped to CPU for safety)
# Adjust the filename to whichever fold had the best Validation ROC
checkpoint_path = "/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd/gat_fold_1.pt"
checkpoint = torch.load(checkpoint_path, map_location='cpu')
state_dict = checkpoint['model_state_dict']

print(f"--- Fold {checkpoint.get('fold', 'Unknown')} Metrics ---")
print(f"Val ROC AUC: {checkpoint.get('val_roc_auc', 0):.4f}")
print(f"Val PR AUC:  {checkpoint.get('val_pr_auc', 0):.4f}\n")

# 2. Extract the weight matrix of the first MLP layer
# Shape is [64, fused_dim]
mlp_weights = state_dict['classifier.0.weight']

# 3. Slice the matrix into Genetics vs. Immunity
genome_weights = mlp_weights[:, :-6]  # Everything EXCEPT the last 6 columns
cell_weights = mlp_weights[:, -6:]    # ONLY the last 6 columns

# 4. Calculate the mean absolute weight for each modality
# We average across all neurons and all features in that block to get a single global score
mean_genome_weight = torch.mean(torch.abs(genome_weights))
mean_cell_weight = torch.mean(torch.abs(cell_weights))

# Calculate the individual cell weights for the breakdown
individual_cell_weights = torch.mean(torch.abs(cell_weights), dim=0)

# 5. Print the Head-to-Head Comparison
total_weight = mean_genome_weight + mean_cell_weight
genome_ratio = (mean_genome_weight / total_weight) * 100
cell_ratio = (mean_cell_weight / total_weight) * 100

print("--- Modality Importance Ratio ---")
print(f"Genetics (GAT Embeddings): {genome_ratio:.2f}%")
print(f"Immunity (Cell Proportions): {cell_ratio:.2f}%\n")

# 6. Print the individual cell breakdown
# Adjust this list to match the exact order of your u_cells columns!
cell_names = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"] 

print("--- Individual Cell Proportion Weights ---")
for name, score in zip(cell_names, individual_cell_weights):
    print(f"{name:10}: {score.item():.6f}")
