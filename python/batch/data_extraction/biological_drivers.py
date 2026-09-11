import torch
import numpy as np
from torch_geometric.data import Data, Batch
from models import GATExplainer
from utils import map_and_export_drivers

def extract_biological_drivers(model_path, m_matrix_df, pheno_df, chr_topologies, cell_cols):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # 1. Unpack the new Dictionary Checkpoint
    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    model = GATExplainer().to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    print(f"Extraction Loaded | Fold {checkpoint['fold']} | Val ROC: {checkpoint['val_roc_auc']:.4f}")
    
    # 2. Isolate PD Patients
    pd_samples = pheno_df[pheno_df['Sample_Group'] == 'PD'].index.tolist()
    
    # 3. Process Patients and Accumulate Weights
    consensus_attention = {c: np.zeros(chr_topologies[c]['n_nodes']) for c in range(1, 23)}
    
    for sample in pd_samples:
        chr_graphs = []
        for c in range(1, 23):
            probes = chr_topologies[c]['probes']
            raw_x = torch.from_numpy(m_matrix_df.loc[sample, probes].values.astype(np.float32)).unsqueeze(1)
            data = Data(x=raw_x, edge_index=chr_topologies[c]['edge_index'], func_type=chr_topologies[c]['func_type'])
            chr_graphs.append(data)
            
        batched_chrs = [Batch.from_data_list([g]).to(device) for g in chr_graphs]
        u = torch.tensor(pheno_df.loc[sample, cell_cols].values.astype(np.float32)).unsqueeze(0).to(device)
        
        with torch.no_grad():
            attention_dict = model(batched_chrs, u)
            
        for c in range(1, 23):
            consensus_attention[c] += attention_dict[c]

    # Average across PD patients
    for c in range(1, 23):
        consensus_attention[c] /= len(pd_samples)
        
    return consensus_attention


