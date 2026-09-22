import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

def main():
    # 1. Load the consensus data
    # Assuming your table is saved as 'consensus_probes.csv'
    df = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/peg1_ppmi_sgpd_common_annotated.csv")

    df = df[1:100]

    # 2. Clean up gene names (converting "PRDM15;PRDM15" to just "PRDM15")
    df['Clean_Gene'] = df['UCSC_RefGene_Name'].astype(str).apply(lambda x: x.split(';')[0])

    # 3. Setup the Matplotlib figure
    sns.set_theme(style="ticks", font_scale=1.2)
    fig, ax = plt.subplots(figsize=(14, 7))

    # Create a clean color palette for the distinct genes
    unique_genes = df['Clean_Gene'].unique()
    colors = sns.color_palette("Set1", len(unique_genes))
    gene_colors = dict(zip(unique_genes, colors))

    # 4. Plot each CpG probe exactly at its physical genomic coordinate
    for idx, row in df.iterrows():
        gene = row['Clean_Gene']
        color = gene_colors[gene]
        x_coord = row['MAPINFO']
        y_mean = row['importance_mean']
        y_std = row['importance_std']
        
        # Extract the raw scores for the three individual cohorts
        cohort_scores = [
            row['compound_importance_ppmi'], 
            row['compound_importance_peg1'], 
            row['compound_importance_sgpd']
        ]
        
        # A. Draw the physical error bars for the standard deviation
        ax.errorbar(x_coord, y_mean, yerr=y_std, fmt='none', ecolor='gray', 
                    elinewidth=1.5, capsize=4, zorder=1)
        
        # B. Plot the solid mean marker
        ax.scatter(x_coord, y_mean, color=color, s=80, edgecolor='black', 
                linewidth=0.8, zorder=3, label=gene if gene not in ax.get_legend_handles_labels()[1] else "")
        
        # C. Overlay the 3 individual cohorts as smaller, semi-transparent dots behind the mean
        ax.scatter([x_coord]*3, cohort_scores, color=color, alpha=0.4, s=30, zorder=2)

    # 5. Format the axes for genomic data
    # Format the X-axis to display in Megabases (Mb) for standard bioinformatics readability
    formatter = FuncFormatter(lambda x, pos: f"{x / 1_000_000:.1f} Mb")
    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlabel("Genomic Position on Chromosome 21 (hg19)", fontsize=14, labelpad=10, fontweight='bold')
    ax.set_ylabel("Cross-Cohort Compound Importance\n(Mean ± SD)", fontsize=14, labelpad=10, fontweight='bold')
    ax.set_title("Top Epigenetic Drivers of sPD Classification on Chromosome 21", fontsize=16, pad=15, fontweight='bold')

    # 6. Clean up the legend (removing duplicate gene entries)
    handles, labels = ax.get_legend_handles_labels()
    unique_labels_dict = dict(zip(labels, handles))
    ax.legend(unique_labels_dict.values(), unique_labels_dict.keys(), 
            title="Target Gene", bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.grid(axis='y', linestyle='--', alpha=0.6)
    sns.despine()
    plt.tight_layout()

    # Save as a high-resolution vector format suitable for LaTeX/PDF thesis rendering
    plt.savefig("chr21_consensus_locus_plot.pdf", format='pdf', bbox_inches='tight')
    plt.show()
    



if __name__ == "__main__":
    main()

