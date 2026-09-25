import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. PRE-PROCESSING
# Load your dataframe (replace 'consensus_probes.csv' with your actual file)
df = pd.read_csv('/Users/kpax/Documents/study/phd/projects/methylation/results/consensus_biomarkers_snr_ranked.csv')
df = df[df['Valid'] == True]  # Filter for valid probes

# Define Methylation State based on the cross-cohort mean effect size
# Hypo: delta_beta < 0 | Hyper: delta_beta > 0
df['Methylation_State'] = np.where(df['delta_beta_mean'] > 0, 'Hypermethylated', 'Hypomethylated')

# Clean the UCSC_RefGene_Group column to take only the primary annotation
# This converts "Body;5'UTR;5'UTR" into just "Body" for clean plotting
df['Regulatory_Region'] = df['UCSC_RefGene_Group'].astype(str).str.split(';').str[0]
df['Regulatory_Region'] = df['Regulatory_Region'].replace('nan', 'Intergenic/Unannotated')

# Force biological chromosome sorting (chr1, chr2 ... chr22) instead of alphabetical (chr1, chr10, chr2)
chr_order = [f'chr{i}' for i in [22, 21, 20, 18, 13, 9]]
df['CHR'] = pd.Categorical(df['CHR'], categories=chr_order, ordered=True)
df = df.dropna(subset=['CHR']) # Drop any weird contigs if they survived filtering

# Set universal publication aesthetics
sns.set_theme(style="ticks", context="paper", font_scale=1.2)
color_map = {'Hypermethylated': '#d62728', 'Hypomethylated': '#1f77b4'} # Red and Blue

# ==========================================
# PLOT 1: Hypo vs Hyper Distribution per Chromosome
# ==========================================
fig1, ax1 = plt.subplots(figsize=(12, 6))

# Create a cross-tabulation of chromosome vs methylation state
state_counts = pd.crosstab(df['CHR'], df['Methylation_State'])
state_counts.plot(kind='bar', stacked=True, color=[color_map['Hypermethylated'], color_map['Hypomethylated']], ax=ax1)

ax1.set_title("Distribution of Differentially Methylated Consensus Probes per Chromosome", fontweight='bold', pad=15)
ax1.set_xlabel("Chromosome", fontweight='bold')
ax1.set_ylabel("Number of Probes", fontweight='bold')
ax1.legend(title="Direction", frameon=False)
sns.despine()
plt.tight_layout()
fig1.savefig("/Users/kpax/Documents/study/phd/projects/methylation/results/fig1_chromosome_distribution.png", dpi=300)


# ==========================================
# PLOT 2: Genomic Locus Mapping (Horizontal Ideogram Style)
# ==========================================
fig2, ax2 = plt.subplots(figsize=(14, 8))

# Draw faint background lines representing the maximum spanned coordinate per chromosome
# This acts as the physical "ideogram" backbone
chr_max_lengths = df.groupby('CHR', observed=True)['MAPINFO'].max()
for y_pos, chromosome in enumerate(chr_order):
    if chromosome in chr_max_lengths:
        ax2.hlines(y=y_pos, xmin=0, xmax=chr_max_lengths[chromosome], color='lightgrey', linewidth=4, zorder=1)

# Overlay the probes at their exact genomic coordinates
sns.stripplot(
    data=df, 
    x='MAPINFO', 
    y='CHR', 
    hue='Methylation_State', 
    palette=color_map,
    jitter=0.2,    # Spreads overlapping probes slightly for density visualization
    size=3,        # Small dot size for resolution
    alpha=0.7, 
    zorder=2,
    ax=ax2
)

# Format the X-axis to display in Megabases (Mb)
def format_mb(x, pos):
    return f"{x * 1e-6:.0f} Mb"
ax2.xaxis.set_major_formatter(plt.FuncFormatter(format_mb))

ax2.set_title("Genomic Mapping of Consensus Epigenetic Drivers", fontweight='bold', pad=15)
ax2.set_xlabel("Genomic Position (Mb)", fontweight='bold')
ax2.set_ylabel("Chromosome", fontweight='bold')
ax2.legend(title="Direction", bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
sns.despine(left=True) # Remove left spine for cleaner ideogram look
plt.tight_layout()
fig2.savefig("/Users/kpax/Documents/study/phd/projects/methylation/results/fig2_genomic_ideograms.png", dpi=300)


# ==========================================
# PLOT 3: Distribution of Functional Groups per Chromosome
# ==========================================
fig3, ax3 = plt.subplots(figsize=(14, 6))

# Create a cross-tabulation of chromosome vs functional region (normalized to 100%)
region_counts = pd.crosstab(df['CHR'], df['Regulatory_Region'], normalize='index') * 100

# Use a qualitative colormap built for discrete categories
region_counts.plot(kind='bar', stacked=True, colormap='Set3', edgecolor='black', linewidth=0.5, ax=ax3)

ax3.set_title("Proportion of Functional Regulatory Regions per Chromosome", fontweight='bold', pad=15)
ax3.set_xlabel("Chromosome", fontweight='bold')
ax3.set_ylabel("Percentage of Probes (%)", fontweight='bold')

# Move the complex legend outside the plot area
ax3.legend(title="Regulatory Region", bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
sns.despine()
plt.tight_layout()
fig3.savefig("/Users/kpax/Documents/study/phd/projects/methylation/results/fig3_regulatory_regions.png", dpi=300)

plt.show()