import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap
import re

def sort_gene_string(gene_str):
    """Sorts gene strings alphabetically to perfectly deduplicate shared pathways."""
    if pd.isna(gene_str):
        return gene_str
    genes = str(gene_str).split(';') 
    return ';'.join(sorted([g.strip() for g in genes]))

def create_plots(cohort_name, file_path, filename="consolidated_pd_network_enrichment_results.csv"):
    df = pd.read_csv(f"{file_path}{filename}")
    
    # Dynamically handle the column name (gseapy defaults to 'Gene_set', but handles 'Database' if renamed)
    db_column = 'Gene_set' if 'Gene_set' in df.columns else 'Database'
    
    # Extract all unique databases from the CSV
    unique_databases = df[db_column].dropna().unique()
    
    print(f"Found {len(unique_databases)} databases. Processing...")

    for db in unique_databases:
        # 1. Isolate the data for the current database
        db_df = df[df[db_column] == db].copy()
        
        # 2. Filter strictly for FDR significance
        fdr_sign = db_df.loc[db_df['Adjusted P-value'] < 0.05].copy()
        
        # Skip this database if nothing survived the FDR threshold
        if fdr_sign.empty:
            print(f"  -> Skipping '{db}': No FDR-significant terms found.")
            continue
            
        # 3. Clean GO IDs and Deduplicate
        fdr_sign['Term'] = fdr_sign['Term'].str.replace(r' \(GO:.*\)', '', regex=True)
        
        if 'Genes' in fdr_sign.columns:
            fdr_sign['Sorted_Genes'] = fdr_sign['Genes'].apply(sort_gene_string)
            fdr_sign = fdr_sign.sort_values('Adjusted P-value', ascending=True)
            fdr_sign = fdr_sign.drop_duplicates(subset=['Sorted_Genes'], keep='first')

        # 4. Extract Top 15 and calculate metrics
        top_terms = fdr_sign.head(15).copy()
        top_terms['-log10(FDR)'] = -np.log10(top_terms['Adjusted P-value'])
        
        if 'Overlap' in top_terms.columns:
            top_terms['Gene_Count'] = top_terms['Overlap'].astype(str).str.split('/').str[0].astype(int)
        else:
            top_terms['Gene_Count'] = 1 

        top_terms['Term'] = top_terms['Term'].map(lambda term: textwrap.fill(str(term), width=20))
        top_terms = top_terms.sort_values('-log10(FDR)', ascending=True)

        # 5. Build the Plot
        fig, ax = plt.subplots(figsize=(10, 8))

        sc = ax.scatter(
            top_terms['-log10(FDR)'],
            top_terms['Term'],
            s=top_terms['Gene_Count'] * 150, 
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
        
        # Format the title for readability (e.g., "KEGG_2021_Human" -> "KEGG 2021 Human")
        clean_title = str(db).replace('_', ' ')
        ax.set_title(f'{clean_title} Enrichment', fontweight='bold', pad=20)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        cbar = fig.colorbar(sc, ax=ax, pad=0.03)
        cbar.set_label(r'$-\log_{10}(\text{FDR})$', rotation=270, labelpad=20, fontweight='bold')

        handles, labels = sc.legend_elements(prop="sizes", alpha=0.6, num=4, func=lambda s: s/150)
        ax.legend(handles, labels, title="Gene Count", bbox_to_anchor=(1.05, 1), loc='lower left', frameon=False)

        plt.tight_layout(rect=[0, 0, 0.95, 1])
        
        # 6. Save uniquely formatted filename and close figure to prevent RAM overload
        safe_db_name = re.sub(r'[^A-Za-z0-9_]', '_', str(db))  # Strips weird characters
        output_filename = f'{file_path}{cohort_name}_{safe_db_name}_FDR.png'
        
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.close() 
        
        print(f"  -> Saved: {cohort_name}_{safe_db_name}_FDR.png")

if __name__ == "__main__":
    create_plots("Consolidated_Cohort", "/Users/kpax/Documents/study/phd/projects/methylation/results/", "final_ranking_gsea_results_specific_genes.csv")