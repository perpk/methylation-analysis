import pandas as pd
import gseapy as gp

dmr_files = { 
    "peg1" : "/workspace/results/peg1/peg1_dmr_results.csv",
    "sgpd" : "/workspace/results/sgpd/GSE145361_DMR_results.csv",
    "ppmi" : "/workspace/results/ppmi/PPMI_DMR_results_df.csv" 
}

results_path = "/workspace/results"

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

for dmr_name, dmr_path in dmr_files.items():
    print(f"Processing DMR file for {dmr_name}: {dmr_path}")

    dmr_df = pd.read_csv(dmr_path)

    overlapping_genes = dmr_df['overlapping_genes']
    master_gene_list = (
        overlapping_genes.dropna()
        .str.split(',')
        .explode()
        .str.strip()
        .loc[lambda genes: genes.ne('')]
        .unique()
    )

    print(f"Running Enrichment against {len(databases)} databases...")
    enrichment = gp.enrichr(
        gene_list=master_gene_list,
        gene_sets=databases,
        organism='human',
        outdir=None 
    )

    results_df = enrichment.results

    # Filter for strict FDR significance
    fdr_sig = results_df[results_df['Adjusted P-value'] < 0.05]
    print(f"Found {len(fdr_sig)} strictly significant terms (FDR < 0.05)")
    print(f"Found {len(results_df[results_df['P-value'] < 0.05])} nominaly significant terms (pvalue < 0.05)")

    if len(fdr_sig) > 0:
        # Plot the top 15 strictly significant pathways
        gp.barplot(
            enrichment.results, 
            column="Adjusted P-value", 
            title="Consolidated PD Network Enrichment",
            top_term=15, 
            figsize=(8, 6),
            ofname=f"{results_path}/{dmr_name}/{dmr_name}_consolidated_pd_network_enrichment.png"
        )
    else:
        print("No terms survived FDR < 0.05. Check nominally significant terms:")
        print(results_df[results_df['P-value'] < 0.01][['Gene_set', 'Term', 'P-value', 'Overlap']].head(10))

    results_df.to_csv(f"{results_path}/{dmr_name}/{dmr_name}_consolidated_pd_network_enrichment_results.csv", index=False)