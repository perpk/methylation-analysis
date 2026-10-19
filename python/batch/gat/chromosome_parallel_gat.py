import torch
import torch.nn as nn
from torch_geometric.nn import GATv2Conv, AttentionalAggregation
import torch.nn.functional as F

class ChromosomeParallelGAT(nn.Module):
    def __init__(self, num_node_classes=7, chr_embed_dim=16, cell_prop_dim=6):
        super(ChromosomeParallelGAT, self).__init__()
        
        self.func_embedding = nn.Embedding(num_embeddings=num_node_classes, embedding_dim=8)
        
        # add_self_loops=False because we explicitly defined them in the topology builder
        # set dropout from nothing/default to 0.6 for regularization
        self.gat1 = GATv2Conv(in_channels=9, out_channels=8, heads=2, concat=True, add_self_loops=False)
        self.gat2 = GATv2Conv(in_channels=16, out_channels=chr_embed_dim, heads=1, concat=True, add_self_loops=False)
        
        self.gate_nn = nn.Sequential(
            nn.Linear(chr_embed_dim, 8),
            nn.ReLU(),
            nn.Linear(8, 1)
        )
        self.pool = AttentionalAggregation(gate_nn=self.gate_nn)
        
        fused_dim = (22 * chr_embed_dim) + cell_prop_dim
        self.classifier = nn.Sequential(
            nn.Linear(fused_dim, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)
        )

    # used with PEG1
    def forward(self, batched_chrs, u_cells, extract_attention=False):
        chr_embeddings = []
        
        # We will store the co-methylation weights here if requested
        attention_dict = {} 

        for c_idx in range(22):
            batch = batched_chrs[c_idx]
            
            emb_func = self.func_embedding(batch.func_type)
            node_feat = torch.cat([batch.x, emb_func], dim=1)
            
            if extract_attention:
                # PyTorch Geometric changes the output to a tuple when return_attention_weights=True
                # It returns: new_node_features, (edge_index, attention_coefficients)
                h, (edges_1, alpha_1) = self.gat1(node_feat, batch.edge_index, return_attention_weights=True)
                h = torch.relu(h)
                
                h, (edges_2, alpha_2) = self.gat2(h, batch.edge_index, return_attention_weights=True)
                h = torch.relu(h)
                
                # Store the edges and their weights, labeled by chromosome (1-22)
                attention_dict[f'chr{c_idx+1}'] = {
                    'layer_1': (edges_1, alpha_1),
                    'layer_2': (edges_2, alpha_2)
                }
            else:
                # Standard forward pass for your 40-epoch training loop (saves VRAM)
                h = torch.relu(self.gat1(node_feat, batch.edge_index))
                h = torch.relu(self.gat2(h, batch.edge_index))
            
            chr_emb = self.pool(h, batch.batch)
            chr_embeddings.append(chr_emb)
            
        genome_vector = torch.cat(chr_embeddings, dim=1)
        fused = torch.cat([genome_vector, u_cells], dim=1)
        logits = self.classifier(fused)
        
        if extract_attention:
            return logits, attention_dict
            
        return logits

    # def forward(self, batched_chrs, u_cells):
    #     chr_embeddings = []
    #     for c_idx in range(22):
    #         batch = batched_chrs[c_idx]
            
    #         emb_func = self.func_embedding(batch.func_type)
    #         node_feat = torch.cat([batch.x, emb_func], dim=1)
            
    #         # 1. Drop 60% of the raw CpG + functional embeddings
    #         node_feat = F.dropout(node_feat, p=0.1, training=self.training)
            
    #         h = F.elu(self.gat1(node_feat, batch.edge_index))
            
    #         # 2. Drop 60% of the hidden embeddings before the second GAT layer
    #         h = F.dropout(h, p=0.1, training=self.training)
            
    #         h = F.elu(self.gat2(h, batch.edge_index))
            
    #         chr_emb = self.pool(h, batch.batch)
    #         chr_embeddings.append(chr_emb)
            
    #     genome_vector = torch.cat(chr_embeddings, dim=1)
    #     fused = torch.cat([genome_vector, u_cells], dim=1)
        
    #     # 3. Drop 60% of the final fused vector before classification
    #     fused = F.dropout(fused, p=0.1, training=self.training)
        
    #     return self.classifier(fused)