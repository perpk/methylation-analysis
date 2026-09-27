import pandas as pd

import matplotlib.pyplot as plt

pheno_df_files = {
    "peg1": "/Users/kpax/Documents/study/phd/projects/methylation/results/peg1/torchgat2/pheno_data.parquet",
    "sgpd": "/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd/pheno_data.parquet",
    "ppmi": "/Users/kpax/Documents/study/phd/projects/methylation/results/ppmi/pheno_data.parquet"
}

# for cohort_name, pheno_file in pheno_df_files.items():
    # pheno_df = pd.read_parquet(pheno_file)
    # print(f"--- Case Distribution for {cohort_name} ---")
    # print(pheno_df.head())
    # print("\n")

    # # Example visualization: distribution of a column (replace 'column_name' with an actual column)
    # plt.figure(figsize=(8, 6))
    # plt.bar(pheno_df['Sample_Group'].value_counts().index, pheno_df['Sample_Group'].value_counts().values)
    # plt.title(f'PD vs Control for {cohort_name}')
    # plt.xlabel('Sample_Group')
    # plt.ylabel('Frequency')
    # plt.tight_layout()
    # plt.savefig(f"/Users/kpax/Documents/study/phd/projects/methylation/results/{cohort_name}/pd_vs_control_{cohort_name}.png", dpi=300, bbox_inches='tight')

for cohort_name, pheno_file in pheno_df_files.items():
    pheno_df = pd.read_parquet(pheno_file)
    if (cohort_name == "ppmi"):
       pheno_df['Age'] = pheno_df['ENROLL_AGE.y']
    print(f"--- Demographics for {cohort_name} ---")
    if 'Age'in pheno_df.columns:
        pheno_df['Age'].dropna(inplace=True)
        pheno_df['Age'] = pheno_df['Age'].astype(float)
        print(f"Mean Age: {pheno_df['Age'].mean():.2f}, Std Dev: {pheno_df['Age'].std():.2f}")
        print(f"Age Distribution:\n{pheno_df['Age'].describe()}\n")
        print(f"Age Distribution PD vs Control:\n{pheno_df.groupby('Sample_Group')['Age'].describe()}\n")
    if 'Sex' in pheno_df.columns:
        print(f"Sex Distribution:\n{pheno_df['Sex'].value_counts()}")
        print(f"Sex Distribution PD vs Control:\n{pheno_df.groupby('Sample_Group')['Sex'].value_counts()}\n")
    print(f"Total Samples: {len(pheno_df)}\n")
    print(f"--- Sample Group Distribution for {cohort_name} ---")
    print(pheno_df['Sample_Group'].value_counts())
