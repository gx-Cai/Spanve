# Spanve: An ultra fast tool for detecting Spatial dependent Expressed Genes.

## Installation

- Install by pip (recommend):

```bash
pip install Spanve
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
command line usage can only run in a standard h5ad file, where there is a `anndata.AnnData.obsm` key named `'spatail'`.

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

see [tutorial notebook](tutorial.ipynb) or [html page](http://htmlpreview.github.io/?https://github.com/gx-Cai/Spanve/blob/main/tutorial.html).