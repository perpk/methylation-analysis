import pandas as pd

# 1. Load the top 500 ranked probes from both datasets
peg1_df = pd.read_csv('/workspace/results/peg1/global_compound_importance_ranking.csv').head(500)
ppmi_df = pd.read_csv('/workspace/results/ppmi/global_compound_importance_ranking.csv').head(500)

# 2. Extract just the Illumina IDs as sets
peg1_probes = set(peg1_df['probe_id'])
ppmi_probes = set(ppmi_df['probe_id'])

# 3. Calculate the exact intersection (probes that appear in the top 500 of BOTH cohorts)
consensus_probes = peg1_probes.intersection(ppmi_probes)

print(f"Number of overlapping probes in Top 500: {len(consensus_probes)}")

# 4. Filter the PPMI dataframe to only contain the consensus probes
consensus_df = ppmi_df[ppmi_df['probe_id'].isin(consensus_probes)].copy()

# Sort them by their importance in PPMI to see the strongest shared drivers
consensus_df = consensus_df.sort_values(by='compound_importance', ascending=False).reset_index(drop=True)

# 5. Check exactly how many of these consensus probes are 'ch.' (non-CpG)
ch_count = sum(consensus_df['probe_id'].str.startswith('ch.'))
print(f"Number of 'ch.' (non-CpG) probes in the consensus: {ch_count}")

# 6. Save this high-confidence signature for pathway analysis
consensus_df.to_csv('/workspace/results/peg1_ppmi_consensus_signature.csv', index=False)

print("\n--- Top 20 Consensus Probes ---")
print(consensus_df[['probe_id', 'chromosome', 'MAPINFO', 'gene_symbol']].head(20))