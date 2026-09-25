import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Load and clean
df = pd.read_csv('/Users/kpax/Documents/study/phd/projects/methylation/results/consensus_biomarkers_snr_ranked.csv')
df = df[df['Valid'] == True]
df['Clean_Gene'] = df['UCSC_RefGene_Name'].astype(str).str.split(';').str[0]
df['Regulatory_Region'] = df['UCSC_RefGene_Group'].astype(str).str.split(';').str[0]
df['Regulatory_Region'] = df['Regulatory_Region'].replace('nan', 'Intergenic')

chromosomes = ['chr9', 'chr13', 'chr18', 'chr20', 'chr21', 'chr22']

for chrom in chromosomes:
    c_df = df[df['CHR'] == chrom].copy()
    if c_df.empty:
        continue
        
    gene_scores = c_df.groupby('Clean_Gene')['signal_to_noise'].max().sort_values(ascending=False)
    top_genes = gene_scores.head(20).index.tolist() 
    
    top_df = c_df[c_df['Clean_Gene'].isin(top_genes)].copy()
    top_df['Direction'] = np.where(top_df['delta_beta_mean'] > 0, '[+]', '[-]')
    
    # ---------------------------------------------------------
    # THE FIX: Single-line string to make the rotated label "thin"
    # ---------------------------------------------------------
    top_df['Gene_Region'] = top_df['Clean_Gene'] + ' (' + top_df['Regulatory_Region'] + ') ' + top_df['Direction']
    
    region_methylation = (
        top_df.groupby('Gene_Region')
        .agg(
            delta_beta_mean=('delta_beta_mean', 'mean'),
            delta_beta_std=('delta_beta_std', 'mean'),
        )
        .sort_values(by='delta_beta_mean', ascending=False)
    )
    
    colors = region_methylation['delta_beta_mean'].ge(0).map({True: 'red', False: 'blue'})
    
    # We can rely on a standard width now that labels are thin
    num_bars = len(region_methylation)
    dynamic_width = max(12.0, num_bars * 0.35)
    
    fig, ax = plt.subplots(figsize=(dynamic_width, 8))
    
    ax.bar(
        region_methylation.index,
        region_methylation['delta_beta_mean'],
        yerr=region_methylation['delta_beta_std'],
        color=colors,
        capsize=4,
        edgecolor='black',
        linewidth=0.8
    )
    
    ax.axhline(0, color='black', linewidth=1)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.1%}'))
    
    ax.set_title(f"Regional Methylation Direction for Top Consensus Genes: {chrom.capitalize()}", fontweight='bold', pad=15)
    ax.set_ylabel('Mean Methylation Change (Delta Beta)', fontweight='bold')
    
    # ---------------------------------------------------------
    # ALIGNMENT FIX: Anchor the right-side of the text to the tick
    # ---------------------------------------------------------
    ax.set_xticks(range(num_bars))
    ax.set_xticklabels(
        region_methylation.index, 
        rotation=90, 
        ha='right',      # Anchors the end of the string to the tick mark
        va='top',     # Centers it perfectly under the bar
        fontsize=10
    )
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    output_filename = f"/Users/kpax/Documents/study/phd/projects/methylation/results/waterfall_{chrom}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()