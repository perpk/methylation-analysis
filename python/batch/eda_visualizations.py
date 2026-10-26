import pandas as pd
import string
import matplotlib.pyplot as plt

pheno_df_files = {
    "peg1": "/Users/kpax/Documents/study/phd/projects/methylation/results/peg1/torchgat2/pheno_data.parquet",
    "sgpd": "/Users/kpax/Documents/study/phd/projects/methylation/results/sgpd/pheno_data.parquet",
    "ppmi": "/Users/kpax/Documents/study/phd/projects/methylation/results/ppmi/pheno_data.parquet"
}


# Dynamically determine the number of columns based on the dictionary length
num_cohorts = len(pheno_df_files)

# Create a 1-row grid with 'num_cohorts' columns. 
# Adjust figsize dynamically (e.g., 6 inches wide per plot)
fig, axes = plt.subplots(1, num_cohorts, figsize=(6 * num_cohorts, 6))

# Ensure 'axes' is always iterable even if there is only 1 cohort
if num_cohorts == 1:
    axes = [axes]

# Zip the dictionary items and the axes array together to iterate over both
for i, ((cohort_name, pheno_file), ax) in enumerate(zip(pheno_df_files.items(), axes)):
    pheno_df = pd.read_parquet(pheno_file)
    
    print(f"--- Case Distribution for {cohort_name} ---")
    print(pheno_df.head(), "\n")

    counts = pheno_df['Sample_Group'].value_counts()
    
    # Plot directly onto the specific axis (ax) instead of plt
    # Adding colors and edge black lines for a cleaner publication look
    bars = ax.bar(counts.index, counts.values, color=["#4C72B0", "#DD8452"], edgecolor="black", width=0.6)
    
    # Add publication panel labels (A, B, C) dynamically using the string module
    panel_label = string.ascii_uppercase[i]
    ax.set_title(panel_label, loc="left", fontsize=20, fontweight="bold", pad=15)
    
    # Add the specific cohort title centered
    ax.set_title(f'PD vs Control: {cohort_name.upper()}', loc="center", fontsize=18, pad=15)
    
    ax.set_xlabel('Sample Group', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Optional: Add exact count numbers on top of the bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height - (height * 0.3),
                f'{int(height)}', ha='center', va='bottom', fontsize=14, weight='bold')

plt.tight_layout()
# Save the combined figure (commented out by default)
plt.savefig("/Users/kpax/Documents/study/phd/projects/methylation/results/pd_vs_control_combined_pub.png", dpi=300, bbox_inches='tight')
# plt.show()

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
