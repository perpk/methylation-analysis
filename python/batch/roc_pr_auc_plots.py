import pandas as pd
from plots import generate_evaluation_plots
import sys

from python.batch.gat.utils import build_chromosome_topologies

def main(m_matrix_path, pheno_path, manifest_path, results_filepath):
    cell_cols = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Gran']

    m_matrix_df = pd.read_parquet(m_matrix_path)
    pheno_df = pd.read_parquet(pheno_path)
    manifest_df = pd.read_parquet(manifest_path)

    if m_matrix_df.index.name != "Sample_Name":
        m_matrix_df = m_matrix_df.set_index("Sample_Name")

    if 'Neu' in pheno_df.columns:
        pheno_df = pheno_df.rename(columns={'Neu': 'Gran'})

    common_probes = list(set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID'])))
    chr_topologies = build_chromosome_topologies(manifest_df, common_probes)

    generate_evaluation_plots(m_matrix_df, pheno_df, chr_topologies, cell_cols, results_filepath)

    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

