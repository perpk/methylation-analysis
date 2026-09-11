import sys
import pandas as pd
import numpy as np
from gat import build_chromosome_topologies
from data_extraction import extract_biological_drivers, map_and_export_drivers

def main(m_matrix_path, pheno_path, manifest_path, results_filepath):
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']

    m_matrix_df = pd.read_parquet(m_matrix_path)
    pheno_df = pd.read_parquet(pheno_path)
    manifest_df = pd.read_parquet(manifest_path)

    common_probes = list(set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID'])))
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)

    top_folds = [1, 5]
    fold_paths = [f"{results_filepath}/gat_fold_{i}.pt" for i in top_folds]

    master_consensus = {c: np.zeros(chr_topologies[c]['n_nodes']) for c in range(1, 23)}
    print("Starting Targeted Extraction on Folds 1 and 5...")

    for model_path in fold_paths:
        fold_consensus = extract_biological_drivers(model_path, m_matrix_df, pheno_df, chr_topologies, cell_cols)
        
        for c in range(1, 23):
            master_consensus[c] += fold_consensus[c]

    for c in range(1, 23):
        master_consensus[c] /= len(top_folds)

    ensemble_drivers_df = map_and_export_drivers(
        master_consensus, 
        chr_topologies, 
        manifest_df, 
        export_path=f"{results_filepath}/pd_top_folds_epigenetic_drivers.csv"
    )

    print("\nTop 20 Consistently Weighted CpG Sites from Generalizable Folds:")
    print(ensemble_drivers_df.head(20))
    

if __name__ == "__main__":
    main(
        m_matrix_path=sys.argv[1],
        pheno_path=sys.argv[2],
        manifest_path=sys.argv[3],
        results_filepath=sys.argv[4]
    )