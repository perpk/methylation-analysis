import pandas as pd
from mapping import map_distal_enhancers_to_genes, export_strict_test_bed, export_strict_background_bed

results_path = "/workspace/results"

# 1. Load the top 500 ranked probes from both datasets
peg1_df = pd.read_csv(f"{results_path}/peg1/global_compound_importance_ranking.csv")
ppmi_df = pd.read_csv(f"{results_path}/ppmi/global_compound_importance_ranking.csv")

manifest_df = pd.read_parquet("/workspace/data/infinium450k_manifest.parquet")

print(peg1_df.head())
print(ppmi_df.head())

# 2. Extract just the Illumina IDs as sets
peg1_probes = set(peg1_df['IlmnID'])
ppmi_probes = set(ppmi_df['IlmnID'])

peg1_common_probes = list(set(peg1_probes).intersection(set(manifest_df['IlmnID'])))
ppmi_common_probes = list(set(ppmi_probes).intersection(set(manifest_df['IlmnID'])))

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

mapped_enhancers_peg1 = map_distal_enhancers_to_genes(f"{results_path}/top_500_peg1_strict.bed", f"{results_path}/peg1")
mapped_enhancers_ppmi = map_distal_enhancers_to_genes(f"{results_path}/top_500_ppmi_strict.bed", f"{results_path}/ppmi")
print(mapped_enhancers_peg1.head())
print(mapped_enhancers_ppmi.head())

mapped_peg1 = set(mapped_enhancers_peg1['cpg_id'])
mapped_ppmi = set(mapped_enhancers_ppmi['cpg_id'])

# 3. Calculate the exact intersection (probes that appear in the top 500 of BOTH cohorts)
consensus_probes = mapped_peg1.intersection(mapped_ppmi)

print(f"Number of overlapping probes in Top 500: {len(consensus_probes)}")

# 4. Filter the PPMI dataframe to only contain the consensus probes
consensus_df = mapped_enhancers_ppmi[mapped_enhancers_ppmi['cpg_id'].isin(consensus_probes)].copy()

# Sort them by their importance in PPMI to see the strongest shared drivers
consensus_df = consensus_df.sort_values(by='compound_importance', ascending=False).reset_index(drop=True)

# 5. Check exactly how many of these consensus probes are 'ch.' (non-CpG)
ch_count = sum(consensus_df['cpg_id'].str.startswith('ch.'))
print(f"Number of 'ch.' (non-CpG) probes in the consensus: {ch_count}")

# 6. Save this high-confidence signature for pathway analysis
consensus_df.to_csv('/workspace/results/peg1_ppmi_consensus_signature.csv', index=False)

print("\n--- Top 20 Consensus Probes ---")
print(consensus_df.head(20))