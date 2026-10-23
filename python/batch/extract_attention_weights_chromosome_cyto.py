import pandas as pd
import os

def export_chromosome_networks(master_nodes_csv, master_edges_csv, target_chromosomes, output_dir):
    """
    Extracts Cytoscape-ready Node and Edge CSVs for specific chromosomes, 
    filtered by a strict Edge SNR threshold to eliminate noise.
    """
    print("Loading Master Dataframes into RAM...")
    df_edges = pd.read_csv(master_edges_csv)
    df_nodes = pd.read_csv(master_nodes_csv)
    df_nodes = df_nodes[df_nodes['Valid'] == True]

    real_snr_values = df_edges.loc[df_edges['Master_SNR'] < 1e7, 'Master_SNR']

    p95 = real_snr_values.quantile(0.95)
    p99 = real_snr_values.quantile(0.99)

    print(f"95th Percentile SNR: {p95:.2f}")
    print(f"99th Percentile SNR: {p99:.2f}")

    # Then use the 99th percentile as your absolute biological threshold
    snr_threshold = p95
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the Illumina cleaning function for the nodes
    def clean_illumina_string(val, default_label="Intergenic"):
        if pd.isna(val) or val == "":
            return default_label
        return list(dict.fromkeys(str(val).split(';')))[0]

    for chr_num in target_chromosomes:
        print(f"\n=== Processing Chromosome {chr_num} ===")
        
        # 1. Isolate the Chromosome and apply the SNR Threshold
        # (This cuts the hundreds of thousands of edges down to the biological core)
        chr_edges = df_edges[
            (df_edges['Chromosome'] == chr_num) & 
            (df_edges['Master_SNR'] >= snr_threshold)
        ].copy()
        
        edge_count = len(chr_edges)
        if edge_count == 0:
            print(f"Skipping Chr {chr_num}: No edges survived an SNR threshold of {snr_threshold}.")
            continue
            
        print(f"Extracted {edge_count} high-confidence edges.")
        
        # 2. Identify the unique CpGs that are part of these surviving edges
        # We need both the Source and Target lists to build the Node skeleton
        surviving_cpgs = set(chr_edges['Source_CpG']).union(set(chr_edges['Target_CpG']))
        
        # 3. Filter the Master Node table to ONLY include these surviving CpGs
        chr_nodes = df_nodes[df_nodes['IlmnID'].isin(surviving_cpgs)].copy()
        
        # 4. Clean the Illumina manifest strings for Cytoscape grouping
        chr_nodes['Target_Gene'] = chr_nodes['UCSC_RefGene_Name'].apply(clean_illumina_string)
        chr_nodes['Functional_Region'] = chr_nodes['UCSC_RefGene_Group'].apply(clean_illumina_string)
        chr_nodes = chr_nodes.drop(columns=['UCSC_RefGene_Name', 'UCSC_RefGene_Group'])
        
        # 5. Export the paired CSVs
        edges_out = os.path.join(output_dir, f"Chr{chr_num}_Edges_SNR_{snr_threshold}.csv")
        nodes_out = os.path.join(output_dir, f"Chr{chr_num}_Nodes_SNR_{snr_threshold}.csv")
        
        chr_edges.to_csv(edges_out, index=False)
        chr_nodes.to_csv(nodes_out, index=False)
        
        print(f"Exported {len(chr_nodes)} nodes to {nodes_out}")
        print(f"Exported {edge_count} edges to {edges_out}")

# ==========================================
# Execution Block
# ==========================================
# Your target chromosomes based on your top 6 genes
target_chromosomes = [9, 13, 18, 20, 21]

def main():
    # Set your SNR threshold. 
    # A threshold of 15-20 is usually a strong starting point to strip background noise.
    export_chromosome_networks(
        master_nodes_csv="/workspace/data/consensus_biomarkers_snr_ranked.csv",
        # master_nodes_csv="/workspace/data/consensus_biomarkers_snr_unfiltered.csv",
        master_edges_csv="/workspace/results/MASTER_3Cohort_Edge_Consensus.csv",
        target_chromosomes=target_chromosomes,
        output_dir="/workspace/results"
    )

if __name__ == "__main__":
    main()