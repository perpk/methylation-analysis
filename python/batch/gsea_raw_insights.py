import pandas as pd

df_tff = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/Fig3_Regulators_raw_results.csv")

def sort_gene_string(gene_str):
    """Sorts gene strings alphabetically to perfectly deduplicate shared pathways."""
    if pd.isna(gene_str):
        return gene_str
    genes = str(gene_str).split(';') 
    return ';'.join(sorted([g.strip() for g in genes]))

genes = df_tff['Genes'].apply(sort_gene_string)