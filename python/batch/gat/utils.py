def build_chromosome_topologies(manifest_df, common_probes, max_linear_dist_bp=1000):
    manifest = manifest_df[manifest_df['IlmnID'].isin(common_probes)].copy()
    manifest['CHR'] = manifest['CHR'].astype(str).str.replace('chr', '')
    manifest = manifest[manifest['CHR'].isin([str(i) for i in range(1, 23)])]
    manifest['MAPINFO'] = manifest['MAPINFO'].astype(int)
    
    chr_topologies = {}
    
    for c in range(1, 23):
        chr_str = str(c)
        sub = manifest[manifest['CHR'] == chr_str].sort_values(by='MAPINFO').reset_index(drop=True)
        probe_list = sub['IlmnID'].tolist()
        probe_to_idx = {p: i for i, p in enumerate(probe_list)}
        positions = sub['MAPINFO'].values
        n_nodes = len(probe_list)
        
        # A. Encode primary functional annotation
        func_labels = np.full(n_nodes, FUNCTIONAL_MAP['Other'], dtype=np.int64)
        for idx, row in sub.iterrows():
            grps = str(row.get('UCSC_RefGene_Group', '')).split(';')
            if grps and grps[0] in FUNCTIONAL_MAP:
                func_labels[idx] = FUNCTIONAL_MAP[grps[0]]
                
        # B. EXPLICIT SELF LOOPS (Fixes Empty Edge Crash & Mathematically sound for GAT)
        edges = [[i, i] for i in range(n_nodes)]
        
        # C. 1D Linear Adjacency Edges
        for i in range(n_nodes - 1):
            if 0 < (positions[i+1] - positions[i]) <= max_linear_dist_bp:
                edges.append([i, i + 1])
                edges.append([i + 1, i])
                
        # D. Shared Genic Annotation Edges
        gene_groups = {}
        for idx, row in sub.iterrows():
            genes = str(row.get('UCSC_RefGene_Name', '')).split(';')
            if genes and genes[0] not in ('', 'nan'):
                gene = genes[0]
                gene_groups.setdefault(gene, []).append(idx)
                
        for members in gene_groups.values():
            if 1 < len(members) <= 50:
                for i in range(len(members)):
                    for j in range(i + 1, len(members)):
                        edges.append([members[i], members[j]])
                        edges.append([members[j], members[i]])
                        
        # Because we initialized with self-loops, 'edges' is guaranteed non-empty
        edge_arr = np.unique(np.array(edges), axis=0).T
        edge_index = torch.tensor(edge_arr, dtype=torch.long)
            
        chr_topologies[c] = {
            'probes': probe_list,
            'edge_index': edge_index,
            'func_type': torch.tensor(func_labels, dtype=torch.long),
            'n_nodes': n_nodes
        }
        
    return chr_topologies

def chromosome_collate_fn(batch):
    batched_chromosomes = []
    graphs_per_sample = [item[0] for item in batch]
    u_tensor = torch.stack([item[1] for item in batch])
    y_tensor = torch.stack([item[2] for item in batch])
    
    for chr_idx in range(22):
        chr_list = [graphs[chr_idx] for graphs in graphs_per_sample]
        batched_chromosomes.append(Batch.from_data_list(chr_list))
        
    return batched_chromosomes, u_tensor, y_tensor