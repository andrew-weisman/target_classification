# Environment information

Tested environment at least to run on already-generated data: `/data/BIDS-HPC/public/software/conda/envs/r_env`

Packages installed using `conda` are in `/data/BIDS-HPC/public/software/conda/envs/r_env/conda_r_env.txt`.

Packages installed using `pip` are in `/data/BIDS-HPC/public/software/conda/envs/r_env/pip_r_env.txt`.

To generate `conda` requirements list from working environment:

```bash
conda list --explicit > conda_r_env.txt
```

To copy the `conda` environment into a brand new one:

```bash
conda create --name r_env --file conda_r_env.txt
```

To generate `pip` requirements list from working environment:

```bash
pip freeze > pip_r_env.txt
```

To copy the `pip` environment into the new one:

```bash
pip install -r pip_r_env.txt
```

`python` libraries we probably want to ensure we have (will be installed using the above requirements files and instructions):

```
import collections, glob, json, matplotlib, numpy, os, pandas, random, seaborn, sklearn, subprocess, sys
```

`R` libraries we probably want to ensure we have (will be installed using the above requirements files and instructions):

```
library("DESeq2")
library("ggplot2")
library("Rtsne")
```

Note the above process will install `pandas` version `1.1.0`, whereas Biowulf's default `python` module has `pandas` version `0.24.2`, which is insufficient. So, use the `conda` environment above.
