import pandas as pd

def prepare_cytoscape_subnetwork(master_nodes_csv, master_edges_csv, manifest_parquet, output_dir, top_n=500):
    """
    Extracts the Top N nodes, annotates them, and extracts the corresponding subgraph 
    from the Master Edge Consensus to generate two perfectly synced Cytoscape CSVs.
    """
    print(f"Extracting Top {top_n} sub-network...")
    
    # ==========================================
    # 1. BUILD THE ANNOTATED NODE TABLE
    # ==========================================
    df_nodes = pd.read_csv(master_nodes_csv)
    manifest = pd.read_parquet(manifest_parquet)
    
    # Isolate the top diagnostic drivers (Nodes)
    df_top_nodes = df_nodes.head(top_n).copy()
    top_cpg_list = df_top_nodes['IlmnID'].tolist()
    
    # Merge and clean the Illumina metadata
    manifest_sub = manifest[['IlmnID', 'UCSC_RefGene_Name', 'UCSC_RefGene_Group']]
    df_annotated_nodes = pd.merge(
        df_top_nodes, manifest_sub, left_on='IlmnID', right_on='IlmnID', how='left'
    ).drop(columns=['IlmnID'])
    
    def clean_illumina_string(val, default_label="Intergenic"):
        if pd.isna(val) or val == "":
            return default_label
        return list(dict.fromkeys(str(val).split(';')))[0]

    print(df_annotated_nodes.head())
        
    df_annotated_nodes['Target_Gene'] = df_annotated_nodes['UCSC_RefGene_Name'].apply(clean_illumina_string)
    df_annotated_nodes['Functional_Region'] = df_annotated_nodes['UCSC_RefGene_Group'].apply(clean_illumina_string)
    df_annotated_nodes = df_annotated_nodes.drop(columns=['UCSC_RefGene_Name', 'UCSC_RefGene_Group'])
    
    # ==========================================
    # 2. BUILD THE FILTERED EDGE TABLE
    # ==========================================
    # Now we load the massive Master Edge network
    df_edges = pd.read_csv(master_edges_csv)
    
    # Strict Subgraph Filter: Keep the edge ONLY if BOTH CpGs are in our Top N list.
    # This prevents the network from exploding into a massive hairball.
    df_subgraph_edges = df_edges[
        (df_edges['Source_CpG'].isin(top_cpg_list)) & 
        (df_edges['Target_CpG'].isin(top_cpg_list))
    ].copy()
    
    # ==========================================
    # 3. EXPORT FOR CYTOSCAPE
    # ==========================================
    nodes_out = f"{output_dir}/Cytoscape_Nodes_Top{top_n}.csv"
    edges_out = f"{output_dir}/Cytoscape_Edges_Top{top_n}.csv"
    
    df_annotated_nodes.to_csv(nodes_out, index=False)
    df_subgraph_edges.to_csv(edges_out, index=False)
    
    print(f"Success!")
    print(f"Nodes exported: {len(df_annotated_nodes)} -> {nodes_out}")
    print(f"Edges exported: {len(df_subgraph_edges)} -> {edges_out}")


def main():
    prepare_cytoscape_subnetwork(
        "/workspace/results/consensus_biomarkers_snr_ranked.csv", 
        "/workspace/results/MASTER_3Cohort_Edge_Consensus.csv", 
        "/workspace/data/infinium450k_manifest.parquet", 
        "/workspace/results", 
        top_n=100
    )

if __name__ == "__main__":
    main()