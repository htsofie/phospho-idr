#%%
import pandas as pd

df = pd.read_excel("../data/raw/rat/14rat_data.xls")

df.sample(n=20, random_state=42).to_csv("../data/processed/rat/test_data.csv", index=False)

df = pd.read_csv("../data/processed/rat/test_data.csv")

pd.set_option('display.max_columns', None)
print(df.columns)

#%%
keep_columns = ['Uniprot', 'ENSEMBL', 'Position', 'Localization Prob', 'Number of Phospho (STY)', 'Amino Acid', 'Sequence Window', 'Modified Sequence', 'Total Brain #', 'Cortex #', 'Brainstem #', 'Cerebellum #', 'Testicle #', 'Pancreas #', 'Stomach #', 'Liver #', 'Fat #', 'Intestine #', 'Kidney #', 'Spleen #', 'Thymus #', 'Lung #', 'Muscle #', 'Heart #', 'Blood #', 'TOTAL MS SUM.1']
df = df[keep_columns]
df.to_csv("../data/processed/rat/test_data.csv", index=False)

print("shape:", df.shape)