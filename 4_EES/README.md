- [Software paper](https://www.nature.com/articles/s41587-020-00803-5)
- [Vignette](https://nbviewer.org/github/KrishnaswamyLab/MELD/blob/main/notebooks/Wagner2018_Chordin_Cas9_Mutagenesis.ipynb)
- [Krishnaswamy Lab projects](https://krishnaswamylab.org/projects)

# Installation

**PHATE only works on Intel x86 Mac** :(

## Installing x86 Python
Follow instructions [at this link](https://towardsdatascience.com/how-to-use-manage-multiple-python-versions-on-an-apple-silicon-m1-mac-d69ee6ed0250) to set up a copy of your Terminal app to run with Rosetta and install x86 Python.

- NOTE 1: Must install a different major version of Python than currently installed on your Mac...
  - If you have v.3.11.4 for example, installing x86 3.11.5 won't work because pip still installs packages in same 3.11 dir, so instead get x86 3.10.
  - pyenv-alias plugin should solve this problem, but I've had no luck getting this to work.
- NOTE 2: MacOS Ventura won't let you duplicate the Terminal app. Instead, you can download & use iTerm.
- Further instructions here for [managing pyenv environments](https://realpython.com/intro-to-pyenv/)

## Installing packages in x86 Python
In your Rosetta terminal, set the global python version to your x86 Python install (change 3.xx.xx below to your version):

```
pyenv global 3.xx.xx
```

Now quit & reopen the Rosetta terminal. This step is important... without this step, my OS would not recognize the correct Python version I had just set.

Now you can install the needed packages. Even if you have previously installed these packages in a different Python version (including through Anaconda), you need to do it again here so that the x86 version is installed.

```
pip install pandas
pip install meld
pip install --user scprep
pip install --user magic-impute
# graphtools is a dependency of one of the above so it should automatically be installed. But if it isn't:
pip install graphtools
```

## Anaconda
If you're working in a standard Anaconda install (not x86 Python), you can't run PHATE. There is probably a way to install an x86 version of Anaconda Python, but I haven't tried that. If you want to run some of this code in Anaconda but omit the PHATE parts, follow these steps to set up a new conda environment:

```
# Create new conda env (replace <newenv> with whatever name you want)
conda create -n newenv
conda activate newenv
# Install pandas
conda install pandas
# Install pip so you can use pip to install the other packages
conda install pip
```

Then use pip to install meld, scprep, graphtools, magic as above.

# Input files

Create the input files (counts matrix and metadata) in R and export as CSV files.

The following is R code. For the **counts matrix**:

```
# Here, my Seurat object is named allRats
# Subset the Seurat object to only contain genes expressed in at least X cells (set to 10 cells here)
subset_gene_list <- rownames(allRats)[Matrix::rowSums(allRats) > 10]
subset_data <- subset(allRats, features = subset_gene_list)
# Convert to matrix and pull counts from DATA slot
subset_mtx <- as.matrix(GetAssayData(subset_data, slot = "data"))
# Write matrix to CSV file
write.table(subset_mtx, "../MELD/allRats_counts_new.csv", sep = ",", row.names = T, col.names = T, quote = F)
```

For **metadata**:

```
subset_metadata <- subset_data@meta.data
write.csv(subset_metadata, "../MELD/allRats_metadata.csv", quote = FALSE)
```
