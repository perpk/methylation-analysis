import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import sys

def main(results_path):
    # ==========================================
    # 1. Extract the Top Annotated Genes
    # ==========================================
    print("Extracting top annotated GAT drivers...")
    drivers_df = pd.read_csv(f"{results_path}/pd_top_folds_epigenetic_drivers.csv")
    mapped_enhancers_df = pd.read_csv(f"{results_path}/top_500_pd_enhancers_strict.bed", sep="\t", header=None, names=['chrom', 'start', 'end', 'target_gene', 'distance_bp'])

    # Isolate probes that HAVE a known gene annotation
    annotated_df = drivers_df[
        drivers_df['UCSC_RefGene_Name'].notna() & 
        drivers_df['UCSC_RefGene_Name'].astype(str).str.lower() != 'nan'
    ]

    top_annotated = annotated_df.sort_values(by='Attention_Score', ascending=False).head(3000) #2000

    annotated_genes = []
    for gene_string in top_annotated['UCSC_RefGene_Name']:
        # Illumina often lists multiple transcripts/genes per probe separated by ';'
        genes = str(gene_string).split(';')
        annotated_genes.extend(genes)

    # ==========================================
    # 2. Consolidate with Mapped Distal Genes
    # ==========================================
    print("Consolidating with mapped distal enhancers...")
    # Apply the 50kb distance filter to the pybedtools output we made earlier
    strict_mapped_df = mapped_enhancers_df[mapped_enhancers_df['distance_bp'] <= 50000]
    distal_genes = strict_mapped_df['target_gene'].dropna().tolist()

    # Combine both lists and convert to a set to ensure unique gene names
    master_gene_list = list(set(annotated_genes + distal_genes))

    print(f" -> Annotated genes found: {len(set(annotated_genes))}")
    print(f" -> Distal enhancer genes found: {len(set(distal_genes))}")
    print(f" -> Total unique genes for GSEA: {len(master_gene_list)}\n")

    # ==========================================
    # 3. Execute Consolidated GSEA
    # ==========================================
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
            ofname=f"{results_path}/consolidated_pd_network_enrichment.png"
        )
    else:
        print("No terms survived FDR < 0.05. Check nominally significant terms:")
        print(results_df[results_df['P-value'] < 0.01][['Gene_set', 'Term', 'P-value', 'Overlap']].head(10))

    results_df.to_csv(f"{results_path}/consolidated_pd_network_enrichment_results.csv", index=False)



if __name__ == "__main__":
    main(results_path=sys.argv[1])