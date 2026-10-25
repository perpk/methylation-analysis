import pandas as pd
import networkx as nx

def generate_manuscript_metrics_table(edges_csv, nodes_csv, output_csv):
    """
    Calculates In-Degree, Out-Degree, Betweenness Centrality, and Louvain Communities,
    and formats the output as a clean, publication-ready CSV table.
    """
    print(f"Loading network data from {edges_csv}...")
    df_edges = pd.read_csv(edges_csv)
    df_nodes = pd.read_csv(nodes_csv)
    
    # 1. Build the Directed Graph
    G = nx.from_pandas_edgelist(
        df_edges, 
        source='Source_CpG', 
        target='Target_CpG', 
        edge_attr='Master_SNR',
        create_using=nx.DiGraph()
    )
    
    # 2. Calculate Degree Metrics
    in_degree_dict = dict(G.in_degree())
    out_degree_dict = dict(G.out_degree())
    
    # 3. Calculate Betweenness Centrality (Unweighted structure)
    print("Calculating Betweenness Centrality...")
    betweenness_dict = nx.betweenness_centrality(G, weight=None)
    
    # 4. Calculate Louvain Communities (requires undirected graph)
    print("Detecting Louvain Communities...")
    G_undirected = G.to_undirected()
    communities = nx.community.louvain_communities(G_undirected, weight='Master_SNR')
    
    community_dict = {}
    for comm_id, node_set in enumerate(communities):
        # We add 1 so communities are 1-indexed (Community 1, Community 2) for the paper
        for node in node_set:
            community_dict[node] = comm_id + 1
            
    # 5. Map metrics back to the Node table
    df_nodes['In_Degree'] = df_nodes['IlmnID'].map(in_degree_dict).fillna(0).astype(int)
    df_nodes['Out_Degree'] = df_nodes['IlmnID'].map(out_degree_dict).fillna(0).astype(int)
    df_nodes['Total_Degree'] = df_nodes['In_Degree'] + df_nodes['Out_Degree']
    df_nodes['Betweenness_Centrality'] = df_nodes['IlmnID'].map(betweenness_dict).fillna(0.0)
    df_nodes['Epigenetic_Module'] = df_nodes['IlmnID'].map(community_dict).fillna(-1).astype(int)
    
    # 6. Format for the Manuscript Table
    # Filter to only include nodes that actually survived the edge thresholding
    df_table = df_nodes[df_nodes['Total_Degree'] > 0].copy()
    
    # Sort by Out-Degree (Master Drivers) then Betweenness (Bottlenecks)
    df_table = df_table.sort_values(by=['Out_Degree', 'Betweenness_Centrality'], ascending=[False, False])
    
    # Round Betweenness to 4 decimal places for clean formatting
    df_table['Betweenness_Centrality'] = df_table['Betweenness_Centrality'].round(4)
    
    # Select and rename columns for the final publication table
    # Adjust 'Target_Gene' and 'Functional_Region' if your column names differ slightly
    cols_to_export = [
        'IlmnID', 'Target_Gene', 'Functional_Region', 
        'In_Degree', 'Out_Degree', 'Betweenness_Centrality', 'Epigenetic_Module'
    ]
    
    # Ensure we only try to export columns that actually exist
    available_cols = [c for c in cols_to_export if c in df_table.columns]
    df_table = df_table[available_cols]
    
    # Rename columns to publication standards
    rename_dict = {
        'IlmnID': 'CpG Site',
        'Target_Gene': 'Gene Symbol',
        'Functional_Region': 'Genomic Region',
        'In_Degree': 'In-Degree (Passenger)',
        'Out_Degree': 'Out-Degree (Driver)',
        'Betweenness_Centrality': 'Betweenness Centrality',
        'Epigenetic_Module': 'Community Module'
    }
    df_table = df_table.rename(columns=rename_dict)
    
    # Export to CSV
    df_table.to_csv(output_csv, index=False)
    print(f"Success! Formatted table saved to {output_csv}")
    
    return df_table

def main():
    df = generate_manuscript_metrics_table(
        edges_csv="/Users/kpax/Documents/study/phd/projects/methylation/results/chromosome_networks/Cytoscape_Edges_Top100.csv",
        nodes_csv="/Users/kpax/Documents/study/phd/projects/methylation/results/Cytoscape_Nodes_Top100.csv",
        output_csv="/Users/kpax/Documents/study/phd/projects/methylation/results/network_metrics_top_100.csv"
    )
    

if __name__ == "__main__":
    main()