- [Software paper](https://www.nature.com/articles/s41587-020-00803-5)
- [Vignette](https://nbviewer.org/github/KrishnaswamyLab/MELD/blob/main/notebooks/Wagner2018_Chordin_Cas9_Mutagenesis.ipynb)
- [Krishnaswamy Lab projects](https://krishnaswamylab.org/projects)

# Installation

**PHATE only works on Intel x86 Mac**

## Installing x86 Python
These are the steps to set up a copy of your Terminal app to run with Rosetta and install x86 Python. The following is modified from [this Medium story](https://towardsdatascience.com/how-to-use-manage-multiple-python-versions-on-an-apple-silicon-m1-mac-d69ee6ed0250), rewritten here (a) in case this link is/becomes inaccessible and (b) to add my own notes.

**1. Install Homebrew**

If you already have Homebrew installed, you're all set. Nothing new here.

Otherwise, see [instructions here](https://brew.sh/) to install (it's literally just 1 command).

**2. Install pyenv with Homebrew**

[pyenv](https://github.com/pyenv/pyenv) allows you to install different versions of Python and switch among them. This is how you will be able to run both arm64 Python in your regular terminal and an x86 version of Python in your Rosetta terminal (even at the same time if you need to!).

Installation is easiest using Homebrew:

```
brew install pyenv
```

After installing, configure pyenv by adding the following code to either ~/.zprofile and ~/.zshrc (if you use zsh shell) or to ~/.bash_profile and ~/.bashrc if you, like me, are a bash holdout. This is easiest by using `nano` to directly edit those files in terminal (e.g., `nano ~/.bashrc`) and add the following lines to the files:

```
##### Add to ~/.bash_profile or ~/.zprofile #####
eval "$(pyenv init --path)"

##### Add to ~/.bashrc or ~/.zshrc #####
if command -v pyenv 1>/dev/null 2>&1; then
    eval "$(pyenv init -)"
fi
```

Now close terminal and reopen it. You are now ready to install and manage Python versions! See [this site](https://realpython.com/intro-to-pyenv/) for more detailed information about using pyenv to manage Python versions. We will use pyenv to install x86 Python below (not yet).

**3. Install Rosetta (maybe)**

The blog linked above says to install Rosetta (which allows ARM Macs to run apps designed on Intel) with the following command in Terminal:

```
softwareupdate --install-rosetta
```

This gave me some errors (Apple M2 Max running macOS 13.4 Ventura). Some websites say Rosetta *does not* come already installed with macOS, but when I tried opening a copy of Terminal (see below) using Rosetta, I had no issues. So you might just try starting at step 2 below and return to step 1 if Rosetta does not work for you.

**4. Create a Rosetta terminal**

If your OS is older than Ventura, you can simply duplicate the Terminal app (right-click, duplicate).

If you're on macOS Ventura and can't duplicate Terminal, download and install a third-party terminal app such as [iterm2](https://iterm2.com/), which you can then duplicate.

Rename your duplicated terminal app to e.g. "iTerm Rosetta" to denote it as the copy that will open with Rosetta by default.

Now right-click the app -> "Get Info" -> check *Open using Rosetta*

This copy of the terminal app will now open with Rosetta by default, meaning the app will be running as if on a machine with x86/Intel architecture. You can test this by opening your new Rosetta terminal and typing `arch`

```
## i386

# If it says arm64 instead, your terminal app is not running Rosetta.
```

**5. Install Homebrew for x86**

This step isn't required for running x86 Python, but it's necessary if you think you might need to install other things specifically in the x86 environment down the line. Since Homebrew is so easy to install (it's literally the same command - see Step 1 above, just **install it in your Rosetta terminal**), you may as well go ahead and do it. By running the install command in the Rosetta terminal, Homebrew will automatically install the x86 version. The main difference is where Homebrew will be installed and where it will put packages you install using `brew`

NOTE: You *don't* need to install pyenv again in the Rosetta terminal.

**6. Configure x86 pyenv and Homebrew**

Now add the following code block to your ~/.bashrc or ~/.zshrc file:

```
# rosetta terminal setup
if [ $(arch) = "i386" ]; then
    alias brew86="/usr/local/bin/brew"
    alias pyenv86="arch -x86_64 pyenv"
fi
```

This code is an `if:then` statement: If your machine architecture is i386, then the following two aliases are set:
- The command `brew86` calls Homebrew installed in the x86 location (instead of typing `/usr/local/bin/brew`)
- The command `pyenv86` calls pyenv under x86 architecture (instead of typing `arch -x86_64 pyenv`)

This just makes it easier for you to call x86 brew (`brew86 install ...`) and x86 pyenv (`pyenv86 install ...`) in your Rosetta terminal.

**7. Install x86 Python!**

ISSUE: In my (limited) experience, you will need to install a different major version of Python than is currently installed on your Mac.

For example, if you have Python 3.11.4 installed (arm64), installing x86 Python 3.11.5 through pyenv in the Rosetta terminal won't work because `pip` will still install packages in the same 3.11 directory and thus *will not* re-install the x86 version of any packages you already have installed in that directory. You could probably remove all packages installed in the 3.11 directory, for example, and then see if x86 installing will work. But there are likely some packages you need in your normal (arm64) Python and x86 Python, such as matplotlib, so the easiest way around this is to just install a different major version, such as x86 3.10. That way you know all the packages and dependencies will be clean x86 installs and you can keep the arm64 3.11 packages separate.

(NOTE: I had no luck with the `pyenv-alias` plugin mentioned in step 5 on the blog, but you can try it.)

Now let's install an x86 version of Python **in the Rosetta terminal**. First, you can show a list of all versions available to you:

```
pyenv86 install --list
```

Since I already have Python 3.11.4 on my system, I chose the most recent version of 3.10 (3.10.13). For our purposes here (MELD, PHATE, etc.), exact version doesn't seem to matter, only that it's x86. Python 3.10.13 worked for me.

```
pyenv86 install -v 3.10.13
```

After the install finishes, run `pyenv versions` and it should list all versions of Python you have installed - arm64 and x86.

To switch to your new x86 Python version, set the global python version to your x86 Python install (change 3.10.13 below to your version, if different):

```
pyenv global 3.10.13
```

Now **quit & reopen the Rosetta terminal**. This step is important... without this step, my OS would not recognize the correct Python version I had just set.

## Installing packages in x86 Python

Now you can install the necessary packages. Even if you have previously installed these packages in a different Python version (including through Anaconda), you need to do it again here so that the x86 version is installed. **Be sure you are in your Rosetta terminal.**

First, install `pip` with x86 Homebrew: `brew86 install pip`

Now install the other packages:

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
