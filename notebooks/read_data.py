#%%
import pandas as pd

#%%
df = pd.read_parquet("../data/processed/cleaned_phosphomouse_complete.parquet")
print(df.head())
print(df.columns)
print(df.shape)
