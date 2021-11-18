# %%
import pandas as pd
import numpy as np
from pathlib import Path

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, OrdinalEncoder

# %% import data (except methylation)
path = Path(r"/Users/anthonystephenson/Downloads/OneDrive_1_11-16-2021")
file_list = path.glob("**/*")

for file in file_list:
    if file.stem == "methylation":
        continue


clinical = pd.read_csv(path.joinpath("clinical.txt"), sep="\t", index_col=0)
cnv = pd.read_csv(path.joinpath("cnv.txt"), sep="\t", index_col=0)
protein = pd.read_csv(path.joinpath("protein.txt"), sep="\t", index_col=0)
mrna = pd.read_csv(path.joinpath("mrna.txt"), sep="\t", index_col=0)
mirna = pd.read_csv(path.joinpath("mirna.txt"), sep="\t", index_col=0)
gene_annotation = pd.read_csv(path.joinpath(
    "gene-annotation.txt"), sep="\t", index_col=0)
mirna_targets = pd.read_csv(path.joinpath(
    "mirna-targets.txt"), sep="\t", index_col=0)
protein_annotation = pd.read_csv(path.joinpath(
    "protein-annotation.txt"), sep="\t", index_col=0)
methylation_annotation = pd.read_csv(
    path.joinpath("methylation-annotation.txt"), sep="\t", index_col=0)

methylation = pd.read_csv(path.joinpath(
    "methylation.txt"), sep="\t", index_col=0, nrows=1000)

# %% miRNA
# we think miRNA and methylation are the most important, so focus on them
# how do we join miRNA targets (score) to sample/participant data?
# mirna columns are participant labels
# join to mirna targets w "precursor" col
# effect on gene expression something like product of score and mirna expression?
mirna_targets = mirna_targets.drop("mirna", axis=1, inplace=True)
mirna_targets.rename({"precursor": "mirna"}, inplace=True, axis=1)
mirna.reset_index().join(mirna_targets, on="mirna")

# %% Methylation

# %% look into data types a bit
clinical_values = {s: clinical[s].unique()
                   for s in clinical if clinical[s].dtype not in ['int64', 'float64']}

# %% process data
# 1. transform non-numeric data
# 2. pivot/transform data for joining

le = LabelEncoder()
oe = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=np.nan)
ohe = OneHotEncoder()

# %% combine to create dataset
protein.T.index.set_names(clinical.index.name)
data = clinical.join(protein.T, on="participant")

# %% set up training/validation data splits


# %% standardise as required

# %% dimensionality reduction
# pick number of components (or similar) by CV?


# %% fit models
