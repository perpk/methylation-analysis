import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import numpy as np


# 1. Load and clean the dataset
df = pd.read_csv('/Users/kpax/Documents/study/phd/projects/methylation/results/consensus_biomarkers_snr_ranked.csv')
df = df[df['Valid'] == True]  # Filter for valid probes

# Clean Gene Names and Regulatory Regions for plotting
df['Clean_Gene'] = df['UCSC_RefGene_Name'].astype(str).str.split(';').str[0]
df['Regulatory_Region'] = df['UCSC_RefGene_Group'].astype(str).str.split(';').str[0]

target_chromosomes = ['chr9', 'chr13', 'chr18', 'chr20', 'chr21', 'chr22']
y_metric = 'signal_to_noise' 

sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

for chrom in target_chromosomes:
    chrom_df = df[df['CHR'] == chrom].copy()
    if chrom_df.empty:
        continue
        
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Assign Hypo/Hyper colors
    chrom_df['Direction'] = np.where(chrom_df['delta_beta_mean'] > 0, 'Hypermethylated', 'Hypomethylated')
    color_map = {'Hypermethylated': '#d62728', 'Hypomethylated': '#1f77b4'}
    
    # Plot the individual CpG probes
    sns.scatterplot(
        data=chrom_df, 
        x='MAPINFO', 
        y=y_metric, 
        hue='Direction', 
        palette=color_map,
        alpha=0.7,
        s=35,
        edgecolor='black',
        linewidth=0.3,
        zorder=2,
        ax=ax
    )
    
    # 2. Visually Wrap the Top Genes
    # Group by gene and find the highest SNR within each gene cluster
    gene_scores = chrom_df.groupby('Clean_Gene')[y_metric].max().sort_values(ascending=False)
    
    # Select the top 6 most significant genes on this chromosome to highlight
    top_genes = gene_scores.head(6).index.tolist()
    
    y_max_overall = chrom_df[y_metric].max()
    texts = []
    
    for gene in top_genes:
        if gene == 'nan':
            continue
            
        gene_data = chrom_df[chrom_df['Clean_Gene'] == gene]
        
        # Calculate the physical span of the gene's probes
        min_pos = gene_data['MAPINFO'].min()
        max_pos = gene_data['MAPINFO'].max()
        
        # If the gene only has 1 probe, pad the span slightly so the highlight is visible
        if min_pos == max_pos:
            pad = 20000  # 20kb visual padding
            min_pos -= pad
            max_pos += pad
            
        # Draw a shaded vertical rectangle covering the gene's region
        ax.axvspan(min_pos, max_pos, color='gold', alpha=0.15, zorder=1)
        
        # Extract unique regulatory zones for this specific gene
        reg_zones = gene_data['Regulatory_Region'].dropna().unique()
        reg_str = ", ".join([r for r in reg_zones if r != 'nan'])
        
        # Check for distal distances from BEDTools
        dist_str = ""
        if 'distance_bp' in gene_data.columns:
            # Drop NaNs and zeros to find true distal elements
            dists = gene_data['distance_bp'].dropna()
            dists = dists[dists > 0]
            if not dists.empty:
                # If there are multiple distances, show the median or max
                dist_val = int(dists.median())
                dist_str = f"\n(+{dist_val} bp)"
        
        # Construct the final dynamic label
        label_text = f"{gene}\n({reg_str}){dist_str}"
        
        # Place the text at the center of the span, near the top of the plot
        mid_pos = (min_pos + max_pos) / 2
        local_y_max = gene_data[y_metric].max()
        
        t = ax.text(
            mid_pos, 
            local_y_max + (y_max_overall * 0.05), # Float slightly above the highest dot
            label_text,
            ha='center',
            va='bottom',
            fontsize=9,
            fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='grey', boxstyle='round,pad=0.3')
        )
        texts.append(t)
    
    # 3. Algorithmically adjust text labels in case two highlighted genes sit too close together
    adjust_text(
        texts, 
        ax=ax, 
        arrowprops=dict(arrowstyle='-', color='grey', lw=0.5, alpha=0.5)
    )
    
    # 4. Standardize Formatting
    def format_mb(x, pos):
        return f"{x * 1e-6:.1f} Mb"
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_mb))
    
    # Expand Y-axis slightly to make room for the text boxes
    ax.set_ylim(bottom=0, top=y_max_overall * 1.25)
    
    ax.set_title(f"Targeted Epigenetic Landscape: {chrom.capitalize()}", fontweight='bold', pad=15)
    ax.set_xlabel("Genomic Position (Mb)", fontweight='bold')
    ax.set_ylabel("Consensus Stability (Signal-to-Noise Ratio)", fontweight='bold')
    
    ax.legend(title="Methylation State", frameon=True, loc='upper right')
    sns.despine()
    
    plt.tight_layout()
    plt.savefig(f"/Users/kpax/Documents/study/phd/projects/methylation/results/annotated_locus_{chrom}.png", dpi=300)
    plt.close()
    print(f"Saved: annotated_locus_{chrom}.png")