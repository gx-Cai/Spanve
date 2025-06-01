# Spanve: A Statistical Method for Detecting Downstream-friendly Spatially Variable Genes in Large-scale Spatial Transcriptomics Data

## Citation
```
@preprint{b.23.SpanveStatistical,
  title = {Spanve: An {{Statistical Method}} to {{Detect Clustering-friendly Spatially Variable Genes}} in {{Large-scale Spatial Transcriptomics Data}}},
  author = {{Guoxin Cai} and {Yichang Chen} and {Shuqing Chen} and {Xun Gu} and {Zhan Zhou}},
  date = {2023-01-01},
  journaltitle = {bioRxiv},
  pages = {2023.02.08.527623},
  doi = {10.1101/2023.02.08.527623},
  url = {http://biorxiv.org/content/early/2023/03/08/2023.02.08.527623.abstract},
}
```

Analysis code for the paper is available at `Evaluate` directory.

## Installation

- Install by pip (recommend):

```bash
pip install Spanve
```

- Install by source

if there are some problems with the pip installation, you can install it from source code.

```bash
git clone https://github.com/gx-Cai/Spanve.git # or download the zip file and unzip
cd Spanve
pip install -e .  # install in editable mode
```


- no install usage:

```bash
# install required packages
cd Evaluate/Softs
pip install spanve_requirements.txt
cp Spanve.py $your_path
```

## Usage

### cli usage

```bash
spanve --help
```

```
Usage: Spanve [OPTIONS]

Options:
  -i, --input_file PATH       input anndata file(h5ad file.)
  -r, --running_mode TEXT     running mode, default is f(c:cluster;
                              i:imputation; f:fitting)
  -s, --save_path PATH        save path
  -v, --verbose BOOLEAN       verbose
  -n, --njobs INTEGER
  -p, --preprocessed BOOLEAN  int preprocessed or not.
  --help                      Show this message and exit.****
```
command line usage can only run in a standard h5ad file, where there is a `anndata.AnnData.obsm` key named `'spatial'`.

### python usage

#### Quick Start

```python
from Spanve import Spanve
adata = sc.read_h5ad('data.h5ad')
spanve = Spanve(adata)
# -- fitting for spatial genes
spanve.fit()
spanve.save('result.csv')
# -- spatial imputation
X = adata_preprocess(adata)
X_ = spanve.impute_from_graph(X[:,spanve.rejects])
```

#### Details

see [tutorial notebook](tutorial.ipynb) or [html page](https://github.com/gx-Cai/Spanve/wiki/Tutorial).


## Known Issues

- When numpy version is 1.24 or higher, you may enrounter the following error:

```
AttributeError: module 'numpy' has no attribute 'int0'
```

solution: downgrade numpy, or **install by source code**. \[will fix in future release\]