import pandas as pd
import numpy as np 

df_tff = pd.read_csv("/Users/kpax/Documents/study/phd/projects/methylation/results/Fig3_Regulators_raw_results.csv")

genes = (
	df_tff['Genes']
	.dropna()
	.str.split(';')
	.explode()
	.str.strip()
	.loc[lambda s: s.ne('')]
	.unique()
	.tolist()
)

print(f"Gene unique count: ", len(genes))

print(genes)