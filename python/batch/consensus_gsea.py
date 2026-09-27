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

    run_preranked = True

    final_ranking = pd.read_csv("/workspace/results/consensus_biomarkers_snr_ranked.csv")

    final_ranking_valid = final_ranking[final_ranking['Valid'] == True]

    final_ranking_genes = split_gene_names(';'.join(final_ranking_valid['UCSC_RefGene_Name'].dropna().unique()))

    final_ranking_gsea_results = "final_ranking_gsea_results"
    if use_specific_genes:
        final_ranking_genes = list(set([gene for gene in final_ranking_genes if gene in genes_to_lookup]))
        final_ranking_gsea_results = "final_ranking_gsea_results_specific_genes"
    if run_preranked:
        final_ranking_gsea_results = "final_ranking_gsea_results_preranked"

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

    if run_preranked:
        df = final_ranking_valid.copy()
        df['Clean_Gene'] = df['UCSC_RefGene_Name'].astype(str).str.split(';').str[0]
        df = df[df['Clean_Gene'] != 'nan']
        gene_ranks = df.groupby('Clean_Gene')['signal_to_noise'].max().reset_index()
        gene_ranks = gene_ranks.sort_values('signal_to_noise', ascending=False).reset_index(drop=True)
        rnk_df = gene_ranks[['Clean_Gene', 'signal_to_noise']]
        print(f"Executing Preranked GSEA on {len(rnk_df)} unique mapped genes...")
        pre_res = gp.prerank(
            rnk=rnk_df, 
            gene_sets=databases, 
            threads=4, 
            permutation_num=1000, 
            outdir=None,
            format='png', 
            seed=42,
            max_size=1000
        )
        results = pre_res.res2d
        sig_results = results[results['FDR q-val'] < 0.05].copy()
        if sig_results.empty:
            print("No pathways reached FDR < 0.05. Consider looking at nominal p-values or adjusting the SNR metric.")
        else:
            print(f"Found {len(sig_results)} strictly significant terms (FDR < 0.05)")
            sig_results = sig_results.sort_values('NES', ascending=False)
            sig_results.to_csv(f"/workspace/results/{final_ranking_gsea_results}_preranked.csv", index=False)

    else:
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

