import glob
import re

import numpy as np
import torch
import matplotlib.pyplot as plt

checkpoint_dir_patterns = {
    "sgpd": "/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd/gat_fold_*.pt",
    "peg1": "/Users/kpax/Documents/study/phd/projects/methylation/results/peg1/torchgat2/gat_fold_*.pt",
    "ppmi": "/Users/kpax/Documents/study/phd/projects/methylation/results/ppmi/gat_fold_*.pt",
}

cell_names = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]

# Trackers
cell_weight_distributions = {}
cell_weight_raw_distributions = {}
modality_ratio_distributions = {}

for dataset, pattern in checkpoint_dir_patterns.items():
    print(f"Processing dataset: {dataset}")
    checkpoint_paths = sorted(
        glob.glob(pattern),
        key=lambda p: int(re.search(r"fold_(\d+)", p).group(1)),
    )

    fold_individual_cell_weights = []
    fold_modality_ratios = []
    raw_cell_weights = []

    for checkpoint_path in checkpoint_paths:
        checkpoint = torch.load(checkpoint_path, map_location='cpu')
        state_dict = checkpoint['model_state_dict']

        print(f"--- Fold {checkpoint.get('fold', 'Unknown')} Metrics ---")
        print(f"Val ROC AUC: {checkpoint.get('val_roc_auc', 0):.4f}")
        print(f"Val PR AUC:  {checkpoint.get('val_pr_auc', 0):.4f}\n")

        mlp_weights = state_dict['classifier.0.weight']
        genome_weights = mlp_weights[:, :-6]
        cell_weights = mlp_weights[:, -6:]
        
        mean_genome_weight = torch.mean(torch.abs(genome_weights))
        mean_cell_weight = torch.mean(torch.abs(cell_weights))

        individual_cell_weights = torch.mean(torch.abs(cell_weights), dim=0)
        fold_individual_cell_weights.append(individual_cell_weights.detach().cpu().numpy())
        raw_cell_weights.append(torch.abs(cell_weights).detach().cpu().numpy())

        total_weight = mean_genome_weight + mean_cell_weight
        genome_ratio = (mean_genome_weight / total_weight) * 100
        cell_ratio = (mean_cell_weight / total_weight) * 100
        
        fold_modality_ratios.append((genome_ratio.item(), cell_ratio.item()))

    fold_individual_cell_weights = np.stack(fold_individual_cell_weights, axis=0)
    fold_modality_ratios = np.array(fold_modality_ratios)

    cell_weight_distributions[dataset] = fold_individual_cell_weights
    cell_weight_raw_distributions[dataset] = np.concatenate(raw_cell_weights, axis=0)
    modality_ratio_distributions[dataset] = fold_modality_ratios

# ---------------------------------------------------------
# Publication-Ready Plotting Setup (Grid Layout)
# ---------------------------------------------------------
datasets = list(modality_ratio_distributions.keys())
num_datasets = len(datasets)

# Create a master figure with a GridSpec layout
fig = plt.figure(figsize=(15, 12))
# Top row for Bar Chart, Bottom row for Boxplots (taller to fit x-labels)
gs = fig.add_gridspec(2, num_datasets, height_ratios=[1, 1.4], hspace=0.45)

color_cell = "#DD8452"
color_genomic = "#4C72B0"
color_immune_boxes = "#81C784" # Soft green to distinguish the cell breakdown

# --- Panel A: Mean Modality Ratio Bar Chart (Row 0, Spanning all columns) ---
ax_bar = fig.add_subplot(gs[0, :])

y_positions = np.arange(num_datasets)
bar_height = 0.35

genome_means = [modality_ratio_distributions[d][:, 0].mean() for d in datasets]
cell_means = [modality_ratio_distributions[d][:, 1].mean() for d in datasets]

bars_cell = ax_bar.barh(
    y_positions - bar_height / 2, cell_means,
    height=bar_height, label="Immunity (Cell)", color=color_cell, edgecolor="black"
)
bars_genome = ax_bar.barh(
    y_positions + bar_height / 2, genome_means,
    height=bar_height, label="Genetics (CpG)", color=color_genomic, edgecolor="black"
)

# Numeric percentage text
for i, (c_mean, g_mean) in enumerate(zip(cell_means, genome_means)):
    ax_bar.text(c_mean + 1, y_positions[i] - bar_height / 2, f"{c_mean:.1f}%", va='center', fontsize=14, color="black")
    ax_bar.text(g_mean + 1, y_positions[i] + bar_height / 2, f"{g_mean:.1f}%", va='center', fontsize=14, color="black")

ax_bar.set_yticks(y_positions)
ax_bar.set_yticklabels([d.upper() for d in datasets], fontsize=16)
ax_bar.set_xlim(0, 105)
ax_bar.set_xlabel("Modality Importance (%)", fontsize=16)
ax_bar.xaxis.set_major_formatter(lambda value, position: f"{value:.0f}%")
ax_bar.tick_params(axis="x", labelsize=14)
ax_bar.grid(axis="x", linestyle="--", alpha=0.6)

ax_bar.set_title("A", loc="left", fontsize=24, fontweight="bold", pad=20)
ax_bar.set_title("Mean Modality Importance Ratio Across Folds", loc="center", fontsize=18)
ax_bar.legend(fontsize=15, loc="lower right", framealpha=1)


# --- Panel B: Individual Cell Boxplots (Row 1, 1 column per dataset) ---
axes_box = []
for i, dataset in enumerate(datasets):
    # Share the y-axis across all boxplots
    if i == 0:
        ax = fig.add_subplot(gs[1, i])
        ax.set_ylabel("Absolute Classifier Weight", fontsize=16)
        ax.set_title("B", loc="left", fontsize=24, fontweight="bold", pad=20)
    else:
        ax = fig.add_subplot(gs[1, i], sharey=axes_box[0])
        plt.setp(ax.get_yticklabels(), visible=False)
    axes_box.append(ax)

    weights = cell_weight_raw_distributions[dataset]
    
    bp = ax.boxplot(
        [weights[:, j] for j in range(len(cell_names))],
        tick_labels=cell_names,
        patch_artist=True,
        widths=0.6
    )
    
    # Apply styling to boxplots
    for patch in bp['boxes']:
        patch.set_facecolor(color_immune_boxes)
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for median in bp['medians']:
        median.set_color("black")
        median.set_linewidth(2)

    ax.set_title(dataset.upper(), loc="center", fontsize=18, pad=15)
    # Rotate x labels so the 6 cell names fit cleanly side-by-side
    ax.tick_params(axis="x", rotation=45, labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.grid(axis="y", linestyle="--", alpha=0.6)


fig.tight_layout()

plt.savefig("/Users/kpax/Documents/study/phd/projects/methylation/results/combined_cell_weights_panel.pdf", format="pdf", bbox_inches="tight")
# plt.show()