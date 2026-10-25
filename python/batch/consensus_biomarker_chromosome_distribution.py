import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/workspace/data/consensus_biomarkers_snr_ranked.csv")
df = df[df['Valid'] == True]

plot_data = df[['CHR', 'signal_to_noise']].dropna()
chromosomes = sorted(plot_data['CHR'].unique(), key=lambda value: str(value))

upper_5 = df['signal_to_noise'].quantile(0.95)
print(f"upper 5% {upper_5:.6e}")
print(df[df['signal_to_noise'] >= upper_5].shape)
print(df[df['signal_to_noise'] >= upper_5]['CHR'].unique())
print(df[df['signal_to_noise'] >= upper_5]['CHR'].value_counts())

robust_probes = df[df['signal_to_noise'] >= upper_5]

def clean_illumina_string(val, default_label="Intergenic"):
        if pd.isna(val) or val == "":
            return default_label
        return list(dict.fromkeys(str(val).split(';')))[0]
        
robust_probes['Target_Gene'] = robust_probes['UCSC_RefGene_Name'].apply(clean_illumina_string)

# Drop NaNs and 'Intergenic' regions, keeping only actual gene names
elite_genes = robust_probes['Target_Gene'].dropna().unique().tolist()
elite_genes = [g for g in elite_genes if g != "Intergenic"]

print(f"Running pathway analysis on {len(elite_genes)} elite diagnostic genes.")

import gseapy as gp

databases = [
    'GO_Biological_Process_2026', 
    'Reactome_2022', 
    'KEGG_2021_Human',
    'Reactome_Pathways_2024',
    'GO_Molecular_Function_2026',
    'TRRUST_Transcription_Factors_2019',
    'TRANSFAC_and_JASPAR_PWMs',
    'ENCODE_TF_ChIP-seq_2015',
    'WikiPathways_2024_Human',
    'SynGO_2024',
    'WikiPathways_2024_Human',
    'Elsevier_Pathway_Collection'
]


enr = gp.enrichr(
    gene_list=elite_genes,
    gene_sets=databases,
    organism='human',
    outdir=None # Keeps it in Pandas rather than writing to disk immediately
)

# Filter for statistically significant pathways (p < 0.05)
df_pathways = enr.results
sig_pathways = df_pathways[df_pathways['P-value'] < 0.05].sort_values('P-value')

print(sig_pathways[['Term', 'Overlap', 'P-value', 'Genes']].head(5))
sig_pathways.to_csv("/workspace/results/sig_pathways.csv")

# fig, ax = plt.subplots(figsize=(12, 6))
# ax.boxplot(
# 	[plot_data.loc[plot_data['CHR'] == chromosome, 'signal_to_noise'] for chromosome in chromosomes],
# 	tick_labels=[str(chromosome) for chromosome in chromosomes], showfliers=False
# )
# ax.set_xlabel('Chromosome (CHR)')
# ax.set_ylabel('signal_to_noise')
# ax.set_title('signal_to_noise distribution by chromosome')
# ax.grid(axis='y', alpha=0.3)
# fig.tight_layout()
# plt.show()

