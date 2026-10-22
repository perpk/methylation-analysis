import pandas as pd

df = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/MASTER_3Cohort_Edge_Consensus.csv")

print(df["Master_SNR"].describe())

real_snr_values = df.loc[df['Master_SNR'] < 1e7, 'Master_SNR']

p95 = real_snr_values.quantile(0.95)
p99 = real_snr_values.quantile(0.99)

print(f"95th Percentile SNR: {p95:.2f}")
print(f"99th Percentile SNR: {p99:.2f}")

print(df[df["Master_SNR"] >= 51.05].shape)

print(df.columns)

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches

# # 1. Simulate the Training Dynamics (Loss Curves)
# epochs = np.arange(1, 101)
# # Training loss continuously drops as the model's backpropagation learns the data
# train_loss = 2.0 * np.exp(-0.06 * epochs) + 0.15
# # Validation loss drops initially, but begins to rise as the model overfits
# val_loss = 1.9 * np.exp(-0.06 * epochs) + 0.25 + 0.0003 * np.maximum(0, epochs - 45)**2

# # Identify the exact epoch where validation loss is minimized (Optimal early stopping)
# optimal_epoch = epochs[np.argmin(val_loss)]
# optimal_val_loss = np.min(val_loss)

# # 2. Set up the Figure Architecture (Dual Panel)
# fig, (ax_split, ax_loss) = plt.subplots(
#     nrows=2, ncols=1, 
#     figsize=(12, 9), 
#     gridspec_kw={'height_ratios': [1, 3]} # Make the loss curve larger
# )

# # ---------------------------------------------------------
# # PANEL A: The Structural 3-Way Data Split
# # ---------------------------------------------------------
# ax_split.barh("Total Dataset", 64, color='#1f77b4', edgecolor='white', label='Training Set (~64%)\nBackprop: ON')
# ax_split.barh("Total Dataset", 16, left=64, color='#ff7f0e', edgecolor='white', label='Validation Set (~16%)\nBackprop: OFF')
# ax_split.barh("Total Dataset", 20, left=80, color='#d62728', edgecolor='white', hatch='//', label='Test Set (20%)\nHidden / Unseen')

# # Formatting Panel A
# ax_split.set_xlim(0, 100)
# ax_split.set_xticks([0, 64, 80, 100])
# ax_split.set_xticklabels(['0%', '64%', '80%', '100%'], fontweight='bold')
# ax_split.set_title('Data Partitioning Within a Single Fold', fontweight='bold', pad=15)
# ax_split.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
# ax_split.spines['top'].set_visible(False)
# ax_split.spines['right'].set_visible(False)
# ax_split.spines['left'].set_visible(False)
# ax_split.spines['bottom'].set_visible(False)
# ax_split.tick_params(left=False)

# # ---------------------------------------------------------
# # PANEL B: The Learning Curves & Early Stopping
# # ---------------------------------------------------------
# ax_loss.plot(epochs, train_loss, color='#1f77b4', linewidth=2.5, label='Training Loss (Calculated on 64%)')
# ax_loss.plot(epochs, val_loss, color='#ff7f0e', linewidth=2.5, label='Validation Loss (Calculated on 16%)')

# # Mark the Early Stopping Point
# ax_loss.axvline(x=optimal_epoch, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
# ax_loss.scatter(optimal_epoch, optimal_val_loss, color='black', s=100, zorder=5)

# # Annotations explaining the network behavior
# ax_loss.annotate(
#     '1. Biological Learning\nBoth losses decrease.', 
#     xy=(20, 1.2), xycoords='data',
#     fontsize=10, ha='center', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
# )
# ax_loss.annotate(
#     '2. Early Stopping\nSave model weights here.', 
#     xy=(optimal_epoch, optimal_val_loss + 0.1), xytext=(optimal_epoch, optimal_val_loss + 0.4),
#     arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
#     fontsize=10, ha='center', fontweight='bold'
# )
# ax_loss.annotate(
#     '3. Overfitting / Memorization\nTrain loss drops, but Val loss rises.', 
#     xy=(80, 0.9), xycoords='data',
#     fontsize=10, ha='center', bbox=dict(boxstyle="round,pad=0.3", fc="#ffe6e6", ec="red", alpha=0.8)
# )

# # Formatting Panel B
# ax_loss.set_xlabel('Training Epochs (Time)', fontweight='bold', labelpad=10)
# ax_loss.set_ylabel('Model Loss (Error)', fontweight='bold', labelpad=10)
# ax_loss.set_title('Why the Validation Set Exists: Preventing Overfitting', fontweight='bold', pad=15)
# ax_loss.legend(loc='upper right', frameon=False)
# ax_loss.grid(True, linestyle='--', alpha=0.3)
# ax_loss.spines['top'].set_visible(False)
# ax_loss.spines['right'].set_visible(False)

# plt.tight_layout(pad=3.0)
# # plt.savefig("cross_validation_architecture.png", dpi=300, bbox_inches='tight')
# plt.show()

# # import glob
# # import re

# # import numpy as np
# # import torch
# # import matplotlib.pyplot as plt

# # checkpoint_dir_patterns = {
# #     "sgpd": "/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd/gat_fold_*.pt",
# #     "peg1": "/Users/kpax/Documents/study/phd/projects/methylation/results/peg1/torchgat2/gat_fold_*.pt",
# #     "ppmi": "/Users/kpax/Documents/study/phd/projects/methylation/results/ppmi/gat_fold_*.pt",
# # }

# # cell_names = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]

# # # Per-dataset, per-fold mean(|weight|) for each cell type: shape (n_folds, n_cells)
# # cell_weight_distributions = {}
# # # Per-dataset, raw |weight| for each cell type across all samples/folds: shape (n_total_samples, n_cells)
# # cell_weight_raw_distributions = {}
# # # Per-dataset, per-fold (genome_ratio, cell_ratio)
# # modality_ratio_distributions = {}
# # # Per-dataset, per-fold (genome_score, cell_score) actual attention scores
# # modality_score_distributions = {}

# # for dataset, pattern in checkpoint_dir_patterns.items():
# #     print(f"Processing dataset: {dataset}")
# #     checkpoint_paths = sorted(
# #         glob.glob(pattern),
# #         key=lambda p: int(re.search(r"fold_(\d+)", p).group(1)),
# #     )

# #     fold_individual_cell_weights = []
# #     fold_modality_ratios = []
# #     fold_modality_scores = []
# #     raw_cell_weights = []
# #     raw_genome_scores = []
# #     raw_cell_scores = []

# #     for checkpoint_path in checkpoint_paths:
# #         checkpoint = torch.load(checkpoint_path, map_location='cpu')
# #         state_dict = checkpoint['model_state_dict']

# #         print(f"--- Fold {checkpoint.get('fold', 'Unknown')} Metrics ---")
# #         print(f"Val ROC AUC: {checkpoint.get('val_roc_auc', 0):.4f}")
# #         print(f"Val PR AUC:  {checkpoint.get('val_pr_auc', 0):.4f}\n")

# #         mlp_weights = state_dict['classifier.0.weight']
# #         genome_weights = mlp_weights[:, :-6]
# #         cell_weights = mlp_weights[:, -6:]
# #         mean_genome_weight = torch.mean(torch.abs(genome_weights))
# #         mean_cell_weight = torch.mean(torch.abs(cell_weights))

# #         individual_cell_weights = torch.mean(torch.abs(cell_weights), dim=0)
# #         fold_individual_cell_weights.append(individual_cell_weights.detach().cpu().numpy())
# #         raw_cell_weights.append(torch.abs(cell_weights).detach().cpu().numpy())

# #         total_weight = mean_genome_weight + mean_cell_weight
# #         genome_ratio = (mean_genome_weight / total_weight) * 100
# #         cell_ratio = (mean_cell_weight / total_weight) * 100
# #         fold_modality_ratios.append((genome_ratio.item(), cell_ratio.item()))
# #         fold_modality_scores.append((mean_genome_weight.item(), mean_cell_weight.item()))
# #         raw_genome_scores.append(torch.abs(genome_weights).flatten().detach().cpu().numpy())
# #         raw_cell_scores.append(torch.abs(cell_weights).flatten().detach().cpu().numpy())

# #         print("--- Individual Cell Proportion Weights ---")
# #         for name, score in zip(cell_names, individual_cell_weights):
# #             print(f"{name:10}: {score.item():.6f}")

# #     fold_individual_cell_weights = np.stack(fold_individual_cell_weights, axis=0)
# #     fold_modality_ratios = np.array(fold_modality_ratios)
# #     fold_modality_scores = np.array(fold_modality_scores)

# #     cell_weight_distributions[dataset] = fold_individual_cell_weights
# #     cell_weight_raw_distributions[dataset] = np.concatenate(raw_cell_weights, axis=0)
# #     modality_ratio_distributions[dataset] = fold_modality_ratios
# #     modality_score_distributions[dataset] = (
# #         np.concatenate(raw_genome_scores, axis=0),
# #         np.concatenate(raw_cell_scores, axis=0),
# #     )

# #     print(f"--- {dataset.upper()} Aggregated Across {len(checkpoint_paths)} Folds ---")
# #     print("Individual Cell Proportion Weights (mean / median / std / p25 / p75):")
# #     for i, name in enumerate(cell_names):
# #         col = fold_individual_cell_weights[:, i]
# #         print(
# #             f"{name:10}: mean={col.mean():.6f}  median={np.median(col):.6f}  "
# #             f"std={col.std():.6f}  p25={np.percentile(col, 25):.6f}  "
# #             f"p75={np.percentile(col, 75):.6f}"
# #         )

# #     genome_col = fold_modality_ratios[:, 0]
# #     cell_col = fold_modality_ratios[:, 1]
# #     print("Modality Importance Ratio (mean / median / std / p25 / p75):")
# #     print(
# #         f"Genetics : mean={genome_col.mean():.2f}%  median={np.median(genome_col):.2f}%  "
# #         f"std={genome_col.std():.2f}%  p25={np.percentile(genome_col, 25):.2f}%  "
# #         f"p75={np.percentile(genome_col, 75):.2f}%"
# #     )
# #     print(
# #         f"Immunity : mean={cell_col.mean():.2f}%  median={np.median(cell_col):.2f}%  "
# #         f"std={cell_col.std():.2f}%  p25={np.percentile(cell_col, 25):.2f}%  "
# #         f"p75={np.percentile(cell_col, 75):.2f}%\n"
# #     )

# # fig, axes = plt.subplots(
# #     1,
# #     len(cell_weight_raw_distributions),
# #     figsize=(5 * len(cell_weight_raw_distributions), 5),
# #     sharey=True,
# #     squeeze=False,
# # )

# # for ax, (dataset, weights) in zip(axes[0], cell_weight_raw_distributions.items()):
# #     ax.boxplot(
# #         [weights[:, i] for i in range(len(cell_names))],
# #         tick_labels=cell_names,
# #     )
# #     n_folds = cell_weight_distributions[dataset].shape[0]
# #     ax.set_title(f"{dataset.upper()}")
# #     ax.set_ylabel("Absolute classifier weight")
# #     ax.tick_params(axis="x", rotation=45)
# #     ax.grid(axis="y", alpha=0.3)

# # fig.suptitle("Individual Cell Proportion Weight Distributions Across Folds")
# # fig.tight_layout()

# # fig2, axes2 = plt.subplots(
# #     1,
# #     len(modality_score_distributions),
# #     figsize=(5 * len(modality_score_distributions), 5),
# #     squeeze=False,
# # )

# # for ax, (dataset, scores) in zip(axes2[0], modality_score_distributions.items()):
# #     genome_scores, cell_scores = scores
# #     labels = ["Cell", "Genomic"]
# #     ax.boxplot(
# #         [cell_scores, genome_scores],
# #         tick_labels=labels,
# #     )
# #     ax.set_ylabel("|attention weight|")
# #     ax.set_title(f"{dataset.upper()}")
# #     ax.grid(axis="y", alpha=0.3)

# # fig2.suptitle("Modality Attention Scores Across Folds")
# # fig2.tight_layout()

# # fig3, ax3 = plt.subplots(figsize=(8, 0.8 + 1.2 * len(modality_ratio_distributions)))

# # datasets = list(modality_ratio_distributions.keys())
# # y_positions = np.arange(len(datasets))
# # bar_height = 0.35

# # genome_means = [modality_ratio_distributions[d][:, 0].mean() for d in datasets]
# # cell_means = [modality_ratio_distributions[d][:, 1].mean() for d in datasets]

# # ax3.barh(
# #     y_positions - bar_height / 2,
# #     cell_means,
# #     height=bar_height,
# #     label="Cell",
# # )
# # ax3.barh(
# #     y_positions + bar_height / 2,
# #     genome_means,
# #     height=bar_height,
# #     label="Genomic",
# # )
# # ax3.set_yticks(y_positions)
# # ax3.set_yticklabels([d.upper() for d in datasets])
# # ax3.set_xlim(0, 100)
# # ax3.set_xlabel("Percentage")
# # ax3.xaxis.set_major_formatter(lambda value, position: f"{value:.0f}%")
# # ax3.set_title("Mean Modality Importance Ratio by Dataset")
# # ax3.legend()
# # ax3.grid(axis="x", alpha=0.3)

# # fig3.tight_layout()

# # plt.show()
