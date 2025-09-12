# %%
import pandas as pd
import os


# %% [markdown]
## Need to edit this to work for multiple data sets through config file
# Convert Excel to Parquet
df = pd.read_excel("data/raw/mouse/Phosphomouse_phosphorylation_sites.xlsb", engine='pyxlsb')
df.to_parquet("data/processed/mouse/Phosphomouse_phosphorylation_sites.parquet", index=False)
print(f"Converted: {df.shape[0]} rows × {df.shape[1]} columns")
print("Column headings:")
print(df.columns.tolist())