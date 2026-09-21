from asyncio import sleep

import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import sys

def split_gene_names(gene_string):
    """
    Splits a string of gene names separated by ';' into a list of individual gene names.
    """
    if pd.isna(gene_string):
        return []
    return [gene.strip() for gene in str(gene_string).split(';') if gene.strip()]

def main():
    N = 1_000
    results_path = "/workspace/results"

    # UCSC_RefGene_Name
    peg1_ppmi_sgpd_common_annotated = pd.read_csv(f"{results_path}/peg1_ppmi_sgpd_common_annotated.csv")

    # target_gene
    peg1_ppmi_sgpd_consensus_signature = pd.read_csv(f"{results_path}/peg1_ppmi_sgpd_consensus_signature.csv")

    # target_gene
    peg1_ppmi_consensus_signature = pd.read_csv(f"{results_path}/peg1_ppmi_consensus_signature.csv")

    peg1_ppmi_sgpd_common_annotated.sort_values(by='importance_mean', ascending=False, inplace=True)
    peg1_ppmi_sgpd_consensus_signature.sort_values(by='importance_mean', ascending=False, inplace=True)
    peg1_ppmi_consensus_signature.sort_values(by='importance_mean', ascending=False, inplace=True)

    peg1_ppmi_sgpd_common_annotated_top = peg1_ppmi_sgpd_common_annotated.head(N)
    peg1_ppmi_sgpd_consensus_signature_top = peg1_ppmi_sgpd_consensus_signature.head(N)
    peg1_ppmi_consensus_signature_top = peg1_ppmi_consensus_signature.head(N)

    peg1_ppmi_sgpd_common_genes = split_gene_names(';'.join(peg1_ppmi_sgpd_common_annotated_top['UCSC_RefGene_Name'].dropna().unique()))
    peg1_ppmi_sgpd_consensus_genes = split_gene_names(';'.join(peg1_ppmi_sgpd_consensus_signature_top['target_gene'].dropna().unique()))
    peg1_ppmi_consensus_genes = split_gene_names(';'.join(peg1_ppmi_consensus_signature_top['target_gene'].dropna().unique()))

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

    to_enrich = {
        "all_cohorts_annotated_probes": peg1_ppmi_sgpd_common_genes,
        "all_cohorts_distal": peg1_ppmi_sgpd_consensus_genes,
        "peg1_ppmi_distal_ch_probes": peg1_ppmi_consensus_genes
    }

    gsea_results = {}
    for name, gene_list in to_enrich.items():
        print(f"Performing GSEA for {name} with {len(gene_list)} genes...")
        enrichment = gp.enrichr(
            gene_list=gene_list, 
            gene_sets=databases, 
            organism='human', 
            outdir=None
        )
        sleep(10)  # To avoid hitting API rate limits
        results_df = enrichment.results
        fdr_sig = results_df[results_df['Adjusted P-value'] < 0.05]
        print(f"{name}: Found {len(fdr_sig)} strictly significant terms (FDR < 0.05)")
        print(f"{name}: Found {len(results_df[results_df['P-value'] < 0.05])} nominaly significant terms (pvalue < 0.05)")

        results_df.to_csv(f"{results_path}/{name}_gsea_results.csv", index=False)

if __name__ == "__main__":
    main()

