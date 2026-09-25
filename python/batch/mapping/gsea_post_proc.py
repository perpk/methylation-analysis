
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap

def create_plots(cohort_name, file_path, filename="consolidated_pd_network_enrichment_results.csv", adjust_cohort_name=False):
    df = pd.read_csv(f"{file_path}{filename}")

    plot_title_cohort = cohort_name.upper()
    if adjust_cohort_name:
        plot_title_cohort = cohort_name.replace("_", " ").capitalize()

    corrected_p_sign = df.loc[df['Adjusted P-value'] < 0.05].copy()
    corrected_p_sign.sort_values('Adjusted P-value', inplace=True, ascending=True)
    corrected_p_sign['-log10(P-value)'] = -np.log10(corrected_p_sign['Adjusted P-value'])
    top_corrected = corrected_p_sign.head(15).copy()
    top_corrected['Term'] = top_corrected['Term'].map(
        lambda term: textwrap.fill(str(term), width=35)
    )

    nominal_p_sign = df.loc[(df['P-value'] < 0.05) & (df['Adjusted P-value'] > 0.05)].copy()
    nominal_p_sign.sort_values('P-value', inplace=True, ascending=True)
    nominal_p_sign['-log10(P-value)'] = -np.log10(nominal_p_sign['P-value'])
    top_terms = nominal_p_sign.head(15).copy()
    top_terms['Term'] = top_terms['Term'].map(
        lambda term: textwrap.fill(str(term), width=35)
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 9))

    sc1 = axes[0].scatter(
        top_corrected['-log10(P-value)'],
        top_corrected['Term'],
        c=top_corrected['-log10(P-value)'],
        cmap='viridis',
        alpha=0.8
    )
    axes[0].grid(True)
    axes[0].set_xlabel('-log10(Adjusted P-value)')
    axes[0].set_title(f'{plot_title_cohort} corrected significant terms')
    fig.colorbar(sc1, ax=axes[0])

    sc2 = axes[1].scatter(
        top_terms['-log10(P-value)'],
        top_terms['Term'],
        c=top_terms['-log10(P-value)'],
        cmap='viridis',
        alpha=0.8
    )
    axes[1].grid(True)
    axes[1].set_xlabel('-log10(P-value)')
    axes[1].set_title(f'{plot_title_cohort} nominally significant terms')
    fig.colorbar(sc2, ax=axes[1])

    plt.tight_layout()
    plt.savefig(f'{file_path}{cohort_name}_significant_terms.png')

def main():    
    create_plots("Consolidated Cohort GSEA Results", "/Users/kpax/Documents/study/phd/projects/methylation/results/", "final_ranking_gsea_results.csv", True)

if __name__ == "__main__":
    main()