import pandas as pd
from mapping import map_distal_enhancers_to_genes, export_strict_test_bed, export_strict_background_bed

results_path = "/workspace/results"

# 1. Load the top 500 ranked probes from both datasets
peg1_df = pd.read_csv(f"{results_path}/peg1/global_compound_importance_ranking.csv")
ppmi_df = pd.read_csv(f"{results_path}/ppmi/global_compound_importance_ranking.csv")
sgpd_df = pd.read_csv(f"{results_path}/sgpd/global_compound_importance_ranking.csv")

manifest_df = pd.read_parquet("/workspace/data/infinium450k_manifest.parquet")

print(peg1_df.head())
print(ppmi_df.head())
print(sgpd_df.head())

# 2. Extract just the Illumina IDs as sets
peg1_probes = set(peg1_df['IlmnID'])
ppmi_probes = set(ppmi_df['IlmnID'])
sgpd_probes = set(sgpd_df['IlmnID'])

peg1_common_probes = list(set(peg1_probes).intersection(set(manifest_df['IlmnID'])))
ppmi_common_probes = list(set(ppmi_probes).intersection(set(manifest_df['IlmnID'])))
sgpd_common_probes = list(set(sgpd_probes).intersection(set(manifest_df['IlmnID'])))

peg1_annotated_probes = set(peg1_df[peg1_df['UCSC_RefGene_Name'].notna()]['IlmnID'])
ppmi_annotated_probes = set(ppmi_df[ppmi_df['UCSC_RefGene_Name'].notna()]['IlmnID'])
sgpd_annotated_probes = set(sgpd_df[sgpd_df['UCSC_RefGene_Name'].notna()]['IlmnID'])

common_annotated_probes = peg1_annotated_probes.intersection(ppmi_annotated_probes).intersection(sgpd_annotated_probes)
print(f"Common annotated probes between peg1, ppmi, and sgpd: {len(common_annotated_probes)}")

peg1_top_ranked_annotated = peg1_df[peg1_df['IlmnID'].isin(peg1_annotated_probes)].sort_values(by='compound_importance', ascending=False)
ppmi_top_ranked_annotated = ppmi_df[ppmi_df['IlmnID'].isin(ppmi_annotated_probes)].sort_values(by='compound_importance', ascending=False)
sgpd_top_ranked_annotated = sgpd_df[sgpd_df['IlmnID'].isin(sgpd_annotated_probes)].sort_values(by='compound_importance', ascending=False)

peg1_top_500_ranked_annotated = peg1_top_ranked_annotated[peg1_top_ranked_annotated['IlmnID'].isin(common_annotated_probes)].head(500)
ppmi_top_500_ranked_annotated = ppmi_top_ranked_annotated[ppmi_top_ranked_annotated['IlmnID'].isin(common_annotated_probes)].head(500)
sgpd_top_500_ranked_annotated = sgpd_top_ranked_annotated[sgpd_top_ranked_annotated['IlmnID'].isin(common_annotated_probes)].head(500)

peg1_top_500_ranked_annotated.to_csv(f"{results_path}/peg1/peg1_top_500_ranked_annotated.csv", index=False)
ppmi_top_500_ranked_annotated.to_csv(f"{results_path}/ppmi/ppmi_top_500_ranked_annotated.csv", index=False)
sgpd_top_500_ranked_annotated.to_csv(f"{results_path}/sgpd/sgpd_top_500_ranked_annotated.csv", index=False)

common_annotated_df = ppmi_top_500_ranked_annotated.copy()
common_annotated_df.rename(columns={'compound_importance': 'compound_importance_ppmi'}, inplace=True)
common_annotated_df = common_annotated_df.merge(peg1_top_500_ranked_annotated[['IlmnID', 'compound_importance']], left_on='IlmnID', right_on='IlmnID', how='left')
common_annotated_df.rename(columns={'compound_importance': 'compound_importance_peg1'}, inplace=True)
common_annotated_df = common_annotated_df.merge(sgpd_top_500_ranked_annotated[['IlmnID', 'compound_importance']], left_on='IlmnID', right_on='IlmnID', how='left')
common_annotated_df.rename(columns={'compound_importance': 'compound_importance_sgpd'}, inplace=True)

print(common_annotated_df.head())

common_annotated_df['importance_mean'] = common_annotated_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].mean(axis=1)
common_annotated_df['importance_median'] = common_annotated_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].median(axis=1)
common_annotated_df['importance_std'] = common_annotated_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].std(axis=1)

common_annotated_df.to_csv(f"{results_path}/peg1_ppmi_sgpd_common_annotated_top_500.csv", index=False)

export_strict_test_bed(
    drivers_csv_path=f"{results_path}/peg1/global_compound_importance_ranking.csv",
    manifest_df=manifest_df,
    top_n=500,
    output_path=f"{results_path}/top_500_peg1_strict.bed",
    sort_by="compound_importance"
)

export_strict_test_bed(
    drivers_csv_path=f"{results_path}/ppmi/global_compound_importance_ranking.csv",
    manifest_df=manifest_df,
    top_n=500,
    output_path=f"{results_path}/top_500_ppmi_strict.bed",
    sort_by="compound_importance"
)

export_strict_test_bed(
    drivers_csv_path=f"{results_path}/sgpd/global_compound_importance_ranking.csv",
    manifest_df=manifest_df,
    top_n=500,
    output_path=f"{results_path}/top_500_sgpd_strict.bed",
    sort_by="compound_importance"
)

export_strict_background_bed(
        manifest_df=manifest_df,
        common_probes=peg1_common_probes,
        output_path=f"{results_path}/peg1_consesus_background_universe_strict.bed"
    )

export_strict_background_bed(
        manifest_df=manifest_df,
        common_probes=ppmi_common_probes,
        output_path=f"{results_path}/ppmi_consensus_background_universe_strict.bed"
    )

export_strict_background_bed(
        manifest_df=manifest_df,
        common_probes=sgpd_common_probes,
        output_path=f"{results_path}/sgpd_consensus_background_universe_strict.bed"
    )

mapped_enhancers_peg1 = map_distal_enhancers_to_genes(f"{results_path}/top_500_peg1_strict.bed", f"{results_path}/peg1")
mapped_enhancers_ppmi = map_distal_enhancers_to_genes(f"{results_path}/top_500_ppmi_strict.bed", f"{results_path}/ppmi")
mapped_enhancers_sgpd = map_distal_enhancers_to_genes(f"{results_path}/top_500_sgpd_strict.bed", f"{results_path}/sgpd")

mapped_enhancers_peg1 = mapped_enhancers_peg1.merge(peg1_df[['IlmnID', 'compound_importance']], left_on='cpg_id', right_on='IlmnID', how='left', suffixes=('', '_importance'))
mapped_enhancers_ppmi = mapped_enhancers_ppmi.merge(ppmi_df[['IlmnID', 'compound_importance']], left_on='cpg_id', right_on='IlmnID', how='left', suffixes=('', '_importance'))
mapped_enhancers_sgpd = mapped_enhancers_sgpd.merge(sgpd_df[['IlmnID', 'compound_importance']], left_on='cpg_id', right_on='IlmnID', how='left', suffixes=('', '_importance'))

print(mapped_enhancers_peg1.head())
print(mapped_enhancers_ppmi.head())
print(mapped_enhancers_sgpd.head())

mapped_peg1 = set(mapped_enhancers_peg1['cpg_id'])
mapped_ppmi = set(mapped_enhancers_ppmi['cpg_id'])
mapped_sgpd = set(mapped_enhancers_sgpd['cpg_id'])

# 3. Calculate the exact intersection (probes that appear in the top 500 of BOTH cohorts)
consensus_probes = mapped_peg1.intersection(mapped_ppmi).intersection(mapped_sgpd)
print(f"There are {len(consensus_probes)} consensus probes between PEG1, PPMI, and SGPD in the top 500 ranked distal enhancers.")

print(f"Number of overlapping probes in Top 500: {len(consensus_probes)}")

# 4. Filter the PPMI dataframe to only contain the consensus probes
consensus_df = mapped_enhancers_ppmi[mapped_enhancers_ppmi['cpg_id'].isin(consensus_probes)].copy()

# Sort them by their importance in PPMI to see the strongest shared drivers
consensus_df = consensus_df.sort_values(by='compound_importance', ascending=False).reset_index(drop=True)

# 5. Check exactly how many of these consensus probes are 'ch.' (non-CpG)
ch_count = sum(consensus_df['cpg_id'].str.startswith('ch.'))
print(f"Number of 'ch.' (non-CpG) probes in the consensus: {ch_count}")

# 6. Save this high-confidence signature for pathway analysis

consensus_df_distal = consensus_df.copy()
consensus_df.rename(columns={'compound_importance': 'compound_importance_ppmi'}, inplace=True)
consensus_df.merge(mapped_enhancers_peg1[['compound_importance']], left_on='cpg_id', right_on='cpg_id', how='left', suffixes=('', '_peg1'), inplace=True)
consensus_df.merge(mapped_enhancers_sgpd[['compound_importance']], left_on='cpg_id', right_on='cpg_id', how='left', suffixes=('', '_sgpd'), inplace=True)

consensus_df['importance_mean'] = consensus_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].mean(axis=1)
consensus_df['importance_median'] = consensus_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].median(axis=1)
consensus_df['importance_std'] = consensus_df[['compound_importance_ppmi', 'compound_importance_peg1', 'compound_importance_sgpd']].std(axis=1)

consensus_df.to_csv('/workspace/results/peg1_ppmi_sgpd_consensus_signature.csv', index=False)

print("\n--- Top 20 Consensus Probes ---")
print(consensus_df.head(20))