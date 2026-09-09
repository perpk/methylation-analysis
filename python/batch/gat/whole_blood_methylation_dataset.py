import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data


LABEL_MAP = {'Control': 0.0, 'PD': 1.0}

class WholeBloodMethylationDataset(Dataset):
    def __init__(self, m_matrix_df, pheno_df, cell_cols, chr_topologies):
        self.sample_ids = pheno_df.index.tolist()
        
        # Explicit label encoding for BCEWithLogitsLoss
        labels = pheno_df['Sample_Group'].map(LABEL_MAP)
        if labels.isna().any():
            unknown_labels = sorted(pheno_df.loc[labels.isna(), 'Sample_Group'].astype(str).unique().tolist())
            raise ValueError(f"Unexpected Sample_Group values: {unknown_labels}")
        self.labels = labels.to_numpy(dtype=np.float32)
        
        self.cell_props = pheno_df[cell_cols].values.astype(np.float32)
        self.chr_topologies = chr_topologies
        
        self.chr_m_values = {}
        for c in range(1, 23):
            probes = chr_topologies[c]['probes']
            # Safe slice: guarantees no missing probes because df was pre-subsetted
            self.chr_m_values[c] = m_matrix_df.loc[self.sample_ids, probes].values.astype(np.float32)
            
    def __len__(self):
        return len(self.sample_ids)
        
    def __getitem__(self, idx):
        chr_graphs = []
        for c in range(1, 23):
            raw_x = torch.from_numpy(self.chr_m_values[c][idx]).unsqueeze(1)
            topo = self.chr_topologies[c]
            
            data = Data(
                x=raw_x,
                edge_index=topo['edge_index'],
                func_type=topo['func_type']
            )
            chr_graphs.append(data)
            
        u = torch.from_numpy(self.cell_props[idx])
        y = torch.tensor(self.labels[idx], dtype=torch.float32)
        return chr_graphs, u, y
