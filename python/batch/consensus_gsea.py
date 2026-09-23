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

    final_ranking = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/consensus_biomarkers_snr_ranked.csv")

    final_ranking_valid = final_ranking[final_ranking['Valid'] == True]

    final_ranking_genes = split_gene_names(';'.join(final_ranking_valid['UCSC_RefGene_Name'].dropna().unique()))

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
    results_df.to_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/final_ranking_gsea_results.csv", index=False)

if __name__ == "__main__":
    main()

