import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from torch_geometric.nn import GATv2Conv, GlobalAttention
from torch_geometric.utils import softmax
from torch_geometric.data import Batch, Data
import torch
import torch.nn as nn
import pandas as pd

# =====================================================================
# 1. THE EXPLAINER ARCHITECTURE
# =====================================================================
# This is identical to your training architecture, but the forward pass 
# is modified to explicitly return the biological attention weights.
class GATExplainer(nn.Module):
    def __init__(self, num_node_classes=7, chr_embed_dim=16, cell_prop_dim=6):
        super(GATExplainer, self).__init__()
        self.func_embedding = nn.Embedding(num_embeddings=num_node_classes, embedding_dim=8)
        self.gat1 = GATv2Conv(in_channels=9, out_channels=8, heads=2, concat=True, add_self_loops=False)
        self.gat2 = GATv2Conv(in_channels=16, out_channels=chr_embed_dim, heads=1, concat=True, add_self_loops=False)
        
        self.gate_nn = nn.Sequential(
            nn.Linear(chr_embed_dim, 8),
            nn.ReLU(),
            nn.Linear(8, 1)
        )
        self.pool = GlobalAttention(gate_nn=self.gate_nn)
        
        fused_dim = (22 * chr_embed_dim) + cell_prop_dim
        self.classifier = nn.Sequential(
            nn.Linear(fused_dim, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)
        )

    def forward(self, batched_chrs, u_cells):
        genome_attention = {}
        
        for c_idx in range(22):
            batch = batched_chrs[c_idx]
            
            # Forward pass through GAT layers
            emb_func = self.func_embedding(batch.func_type)
            node_feat = torch.cat([batch.x, emb_func], dim=1)
            h = torch.relu(self.gat1(node_feat, batch.edge_index))
            h = torch.relu(self.gat2(h, batch.edge_index))
            
            # --- EXTRACT BIOLOGICAL IMPORTANCE ---
            # Calculate the raw importance score for every single CpG site
            raw_gate_scores = self.gate_nn(h).view(-1)
            # Normalize scores so they sum to 1.0 across the chromosome
            cpg_attention = softmax(raw_gate_scores, batch.batch)
            
            # Store the attention weights for this chromosome
            # Shape: [num_nodes_in_batch]
            genome_attention[c_idx + 1] = cpg_attention.detach().cpu().numpy()
            
        return genome_attention

# =====================================================================
# 2. EXTRACTION LOGIC
# =====================================================================
def extract_biological_drivers(model_path, m_matrix_df, pheno_df, chr_topologies, cell_cols):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Extracting on: {device}")
    
    # 1. Initialize and Load Weights
    model = GATExplainer().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()
    
    # 2. Isolate Parkinson's Patients (We want to see what drives the PD signature)
    pd_samples = pheno_df[pheno_df['Sample_Group'] == 'PD'].index.tolist()
    print(f"Extracting drivers across {len(pd_samples)} PD patients...")
    
    # 3. Process Patients (Batch size 1 for clean extraction)
    # We will accumulate the attention weights to find the consensus drivers
    consensus_attention = {c: np.zeros(chr_topologies[c]['n_nodes']) for c in range(1, 23)}
    
    for sample in pd_samples:
        chr_graphs = []
        for c in range(1, 23):
            probes = chr_topologies[c]['probes']
            raw_x = torch.from_numpy(m_matrix_df.loc[sample, probes].values.astype(np.float32)).unsqueeze(1)
            data = Data(x=raw_x, edge_index=chr_topologies[c]['edge_index'], func_type=chr_topologies[c]['func_type'])
            chr_graphs.append(data)
            
        # Convert to PyG Batch format (batch size 1)
        batched_chrs = [Batch.from_data_list([g]).to(device) for g in chr_graphs]
        u = torch.tensor(pheno_df.loc[sample, cell_cols].values.astype(np.float32)).unsqueeze(0).to(device)
        
        # 4. Extract Attention
        with torch.no_grad():
            attention_dict = model(batched_chrs, u)
            
        for c in range(1, 23):
            consensus_attention[c] += attention_dict[c]

    # 5. Average the attention across all PD patients
    for c in range(1, 23):
        consensus_attention[c] /= len(pd_samples)
        
    return consensus_attention

# =====================================================================
# 3. MAP TO BIOLOGY & EXPORT
# =====================================================================
def map_and_export_drivers(consensus_attention, chr_topologies, manifest_df, export_path="pd_epigenetic_drivers.csv"):
    records = []
    
    for c in range(1, 23):
        probes = chr_topologies[c]['probes']
        weights = consensus_attention[c]
        
        for idx, probe in enumerate(probes):
            records.append({
                'IlmnID': probe,
                'Chromosome': c,
                'Attention_Score': weights[idx]
            })
            
    df_weights = pd.DataFrame(records)
    
    # Merge with Manifest to get Genes and Regions
    final_df = df_weights.merge(
        manifest_df[['IlmnID', 'UCSC_RefGene_Name', 'UCSC_RefGene_Group']], 
        on='IlmnID', 
        how='left'
    )
    
    # Sort by highest attention score
    final_df = final_df.sort_values(by='Attention_Score', ascending=False).reset_index(drop=True)
    final_df.to_csv(export_path, index=False)
    print(f"Top drivers exported to {export_path}")
    
    return final_df

