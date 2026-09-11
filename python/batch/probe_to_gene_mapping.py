import pandas as pd
import sys
from mapping import map_distal_enhancers_to_genes, export_strict_test_bed, export_strict_background_bed

def main(m_matrix_path, manifest_path, results_path):
    mapped_enhancers_df = map_distal_enhancers_to_genes(f"{results_path}/top_500_pd_enhancers_strict.bed")
    print(mapped_enhancers_df.head())

    manifest_df = pd.read_parquet(manifest_path)
    m_matrix_df = pd.read_parquet(m_matrix_path)

    common_probes = list(set(m_matrix_df.columns).intersection(set(manifest_df['IlmnID'])))

    test_bed = export_strict_test_bed(
        drivers_csv_path=f"{results_path}/pd_ensemble_epigenetic_drivers.csv",
        manifest_df=manifest_df,
        top_n=500,
        output_path=f"{results_path}/top_500_pd_enhancers_strict.bed"
    )

    bg_bed = export_strict_background_bed(
        manifest_df=manifest_df,
        common_probes=common_probes,
        output_path=f"{results_path}/gat_background_universe_strict.bed"
    )

if __name__ == "__main__":
    main(
        m_matrix_path=sys.argv[1],
        manifest_path=sys.argv[2],
        results_path=sys.argv[3]
    )