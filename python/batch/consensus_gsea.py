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

    genes_to_lookup = ["DYRK1A", "GNAS", "KCNQ2", "RNF160", "NCAM2", "PRODH", "PTGES", "SH2D3C", "AP1B1", "GPC6", "SEC11C", "LPIN2"]

    use_specific_genes = True

    final_ranking = pd.read_csv("/workspace/results/consensus_biomarkers_snr_ranked.csv")

    final_ranking_valid = final_ranking[final_ranking['Valid'] == True]

    final_ranking_genes = split_gene_names(';'.join(final_ranking_valid['UCSC_RefGene_Name'].dropna().unique()))

    final_ranking_gsea_results = "final_ranking_gsea_results"
    if use_specific_genes:
        final_ranking_genes = list(set([gene for gene in final_ranking_genes if gene in genes_to_lookup]))
        final_ranking_gsea_results = "final_ranking_gsea_results_specific_genes"

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

    print(f"Performing GSEA for with {len(final_ranking_genes)} genes...")
    enrichment = gp.enrichr(
        gene_list=final_ranking_genes, 
        gene_sets=databases, 
        organism='human', 
        outdir=None
    )
    results_df = enrichment.results
    fdr_sig = results_df[results_df['Adjusted P-value'] < 0.05]
    print(f"Found {len(fdr_sig)} strictly significant terms (FDR < 0.05)")
    print(f"Found {len(results_df[results_df['P-value'] < 0.05])} nominaly significant terms (pvalue < 0.05)")
    results_df = results_df[~results_df['Term'].str.contains('mouse', case=False, na=False)].copy()
    results_df.to_csv(f"/workspace/results/{final_ranking_gsea_results}.csv", index=False)

if __name__ == "__main__":
    main()

