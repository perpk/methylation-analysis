import pandas as pd

df1 = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/peg1_ppmi_sgpd_common_annotated.csv")
df2 = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/peg1_ppmi_sgpd_consensus_signature.csv")

# 1. Drop the redundant or empty 'IlmnID' from df2 if it exists alongside 'cpg_id'
if 'IlmnID' in df2.columns and 'cpg_id' in df2.columns:
    df2 = df2.drop(columns=['IlmnID'])

# 2. Rename df2's unique columns to perfectly match df1's naming convention
rename_mapping = {
    'cpg_id': 'IlmnID',
    'chr': 'CHR',
    'start': 'MAPINFO',
    'target_gene': 'UCSC_RefGene_Name'  # Aligning the gene annotation column
}
df2_aligned = df2.rename(columns=rename_mapping)

# 3. Stack the rows vertically
# ignore_index=True ensures the final DataFrame has a clean, continuous index from 0 to N
df = pd.concat([df1, df2_aligned], axis=0, ignore_index=True)

# 4. Optional: If the two files have overlapping probes and you only want unique ones, 
# you can drop duplicates based on the Illumina ID:
# df = df.drop_duplicates(subset=['IlmnID'], keep='first').reset_index(drop=True)

print(df.head())

mean_threshold = df['importance_mean'].quantile(0.95)
print(f"Mean importance threshold (Top 5%): {mean_threshold:.6e}")

robust_probes = df[df['importance_mean'] >= mean_threshold].copy()

final_ranking = robust_probes.sort_values(by='signal_to_noise', ascending=False).reset_index(drop=True)

print(f"Total probes evaluated: {len(df)}")
print(f"Probes passing biological threshold (Top 5%): {len(final_ranking)}\n")
print("--- Top Consensus Biomarkers (Ranked by SNR) ---")
print(final_ranking.head(20))

# Incorporating delta-beta values per cohort and mean, median, std-deviation across cohorts alongside the SNR if one of the three cohorts has a delta-beta value is signed differently then an additional column named Valid shall contain FALSE otherwise TRUE

df_peg1 = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/peg1_beta_means.csv", index_col=0)
df_peg1['IlmnID'] = df_peg1.index
df_peg1 = df_peg1[['IlmnID', 'delta_beta']]
df_peg1.rename(columns={'delta_beta': 'delta_beta_peg1'}, inplace=True)

df_sgpd = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd_beta_means.csv", index_col=0)
df_sgpd['IlmnID'] = df_sgpd.index
df_sgpd = df_sgpd[['IlmnID', 'delta_beta']]
df_sgpd.rename(columns={'delta_beta': 'delta_beta_sgpd'}, inplace=True)

df_ppmi = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/ppmi_beta_means.csv", index_col=0)
df_ppmi['IlmnID'] = df_ppmi.index
df_ppmi = df_ppmi[['IlmnID', 'delta_beta']]
df_ppmi.rename(columns={'delta_beta': 'delta_beta_ppmi'}, inplace=True)

final_ranking = final_ranking.merge(df_peg1, on='IlmnID', how='left')
final_ranking = final_ranking.merge(df_sgpd, on='IlmnID', how='left')
final_ranking = final_ranking.merge(df_ppmi, on='IlmnID', how='left')

final_ranking['delta_beta_mean'] = final_ranking[['delta_beta_peg1', 'delta_beta_sgpd', 'delta_beta_ppmi']].mean(axis=1)
final_ranking['delta_beta_median'] = final_ranking[['delta_beta_peg1', 'delta_beta_sgpd', 'delta_beta_ppmi']].median(axis=1)
final_ranking['delta_beta_std'] = final_ranking[['delta_beta_peg1', 'delta_beta_sgpd', 'delta_beta_ppmi']].std(axis=1)
final_ranking['delta_beta_signal_to_noise'] = final_ranking['delta_beta_mean'] / final_ranking['delta_beta_std']

cols_to_check = ['delta_beta_peg1', 'delta_beta_sgpd', 'delta_beta_ppmi']

final_ranking['Valid'] = (final_ranking[cols_to_check] > 0).all(axis=1) | (final_ranking[cols_to_check] < 0).all(axis=1)

print(f"Probes with consistent delta-beta direction across cohorts: {final_ranking['Valid'].sum()}")

final_ranking.to_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/consensus_biomarkers_snr_ranked.csv", index=False)

