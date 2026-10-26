import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap
import re
import math
import string
import gseapy as gp
import os

def sort_gene_string(gene_str):
    """Sorts gene strings alphabetically to perfectly deduplicate shared pathways."""
    if pd.isna(gene_str):
        return gene_str
    genes = str(gene_str).split(';') 
    return ';'.join(sorted([g.strip() for g in genes]))

def perform_gsea_and_plot_panels(gene_list, databases, output_dir, file_prefix, title="Enrichment", size_multiplier=150):
    """
    Executes GSEA via gseapy, filters for FDR < 0.05, deduplicates, and generates 
    a multi-panel (A, B, C) publication-ready bubble plot.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Clean the gene list (remove NaNs, 'Intergenic')
    clean_genes = [g for g in gene_list if pd.notna(g) and g != "Intergenic"]
    print(f"Executing Enrichr on {len(clean_genes)} genes across {len(databases)} databases...")

    # 1. Execute gseapy Enrichr
    enr = gp.enrichr(
        gene_list=clean_genes,
        gene_sets=databases,
        organism='human',
        outdir=None # Keep in memory
    )
    
    df = enr.results
    
    # Save the raw master CSV for your supplementary data
    raw_csv_path = os.path.join(output_dir, f"{file_prefix}_raw_results.csv")
    df.to_csv(raw_csv_path, index=False)
    print(f"Raw results saved to {raw_csv_path}")

    # 2. Filter for databases that actually have FDR-significant hits
    valid_dbs = []
    for db in databases:
        db_df = df[(df['Gene_set'] == db) & (df['Adjusted P-value'] < 0.05)]
        if not db_df.empty:
            valid_dbs.append(db)
            
    n_panels = len(valid_dbs)
    if n_panels == 0:
        print("No databases had FDR < 0.05. Exiting plotting.")
        return

    print(f"Found {n_panels} databases with significant results. Building panels...")

    # 3. Dynamic Grid Layout (2 columns wide, as many rows as needed)
    cols = 2 if n_panels > 1 else 1
    rows = math.ceil(n_panels / cols)
    
    # Create the figure
    fig, axes = plt.subplots(rows, cols, figsize=(10 * cols, 8 * rows))
    
    # Flatten axes array for easy iteration, handling single-panel edge cases
    if n_panels > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
        
    panel_labels = string.ascii_uppercase # A, B, C, D...

    # 4. Process and Plot Each Valid Database
    for i, db in enumerate(valid_dbs):
        ax = axes[i]
        
        # Isolate and clean data
        fdr_sign = df[(df['Gene_set'] == db) & (df['Adjusted P-value'] < 0.05)].copy()
        fdr_sign['Term'] = fdr_sign['Term'].str.replace(r' \(GO:.*\)', '', regex=True)
        
        if 'Genes' in fdr_sign.columns:
            fdr_sign['Sorted_Genes'] = fdr_sign['Genes'].apply(sort_gene_string)
            fdr_sign = fdr_sign.sort_values('Adjusted P-value', ascending=True)
            fdr_sign = fdr_sign.drop_duplicates(subset=['Sorted_Genes'], keep='first')

        # Extract Top 15 and calculate metrics
        top_terms = fdr_sign.head(15).copy()
        top_terms['-log10(FDR)'] = -np.log10(top_terms['Adjusted P-value'])
        
        if 'Overlap' in top_terms.columns:
            top_terms['Gene_Count'] = top_terms['Overlap'].astype(str).str.split('/').str[0].astype(int)
        else:
            top_terms['Gene_Count'] = 1 

        top_terms['Term'] = top_terms['Term'].map(lambda term: textwrap.fill(str(term), width=30))
        top_terms = top_terms.sort_values('-log10(FDR)', ascending=True) # Sort for bottom-to-top plotting

        # Build the Scatter Plot
        sc = ax.scatter(
            top_terms['-log10(FDR)'],
            top_terms['Term'],
            s=top_terms['Gene_Count'] * size_multiplier, 
            c=top_terms['-log10(FDR)'],
            cmap='viridis',
            alpha=0.9,
            edgecolor='black',
            linewidth=0.6,
            zorder=3
        )
        
        ax.grid(axis='y', linestyle='--', alpha=0.6, zorder=1)
        ax.grid(axis='x', linestyle='--', alpha=0.3, zorder=1)
        ax.set_xlabel(r'$-\log_{10}(\text{Adjusted P-value})$', fontweight='bold', labelpad=10)
        
        # Format title and panel label
        clean_title = str(db).replace('_', ' ')
        ax.set_title(f'{clean_title}', fontweight='bold', pad=15)
        
        # Add the A, B, C label to the top left of the subplot
        ax.annotate(f"{panel_labels[i]}", xy=(-0.4, 1.05), xycoords='axes fraction',
                    fontsize=20, fontweight='bold', va='bottom', ha='right')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Attach Legends specific to this panel
        cbar = fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
        cbar.set_label(r'$-\log_{10}(\text{FDR})$', rotation=270, labelpad=20, fontweight='bold')

        handles, labels = sc.legend_elements(prop="sizes", alpha=0.6, num=4, func=lambda s: s/size_multiplier)
        ax.legend(handles, labels, title="Gene Count", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

    # Hide any unused subplots (if you have an odd number of databases in a 2-column grid)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # 5. Save the Multi-Panel Figure
    plt.tight_layout(rect=[0, 0, 1, 1])
    
    # Save as both High-Res PNG (for immediate viewing) and PDF (for journal submission)
    png_path = os.path.join(output_dir, f'{file_prefix}_MultiPanel.png')
    pdf_path = os.path.join(output_dir, f'{file_prefix}_MultiPanel.pdf')
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    plt.close() 
    
    print(f"Successfully generated multi-panel figures:\n  -> {png_path}\n  -> {pdf_path}")


if __name__ == "__main__":
    def clean_illumina_string(val, default_label="Intergenic"):
        if pd.isna(val) or val == "":
            return default_label
        return list(dict.fromkeys(str(val).split(';')))[0]    
    df = pd.read_csv("/workspace/data/consensus_biomarkers_snr_ranked.csv")
    df = df[df['Valid'] == True]
    df['Target_Gene'] = df['UCSC_RefGene_Name'].apply(clean_illumina_string)

    elite_genes = df['Target_Gene'].unique().tolist()
    
    output_directory = "/workspace/results"
    
    # Run 1: The Regulatory Networks Figure (Figure 3 in your paper)
    regulator_dbs = [
        'ENCODE_TF_ChIP-seq_2015',
        'TRRUST_Transcription_Factors_2019',
        'TRANSFAC_and_JASPAR_PWMs'
    ]
    perform_gsea_and_plot_panels(
        gene_list=elite_genes, 
        databases=regulator_dbs, 
        output_dir=output_directory, 
        file_prefix="Fig3_Regulators",
        size_multiplier=1
    )
    
    # Run 2: The Biological Pathways Figure (Figure 4 in your paper)
    pathway_dbs = [
        'GO_Biological_Process_2026', 
        'KEGG_2021_Human',
        'Reactome_2022',
        'SynGO_2024'
    ]
    perform_gsea_and_plot_panels(
        gene_list=elite_genes, 
        databases=pathway_dbs, 
        output_dir=output_directory, 
        file_prefix="Fig4_Pathways"
    )