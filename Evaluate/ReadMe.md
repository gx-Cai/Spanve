# Evaluation pipeline.

## Data Sources

The data used in this study includes the following:
- Datasets from the [10X Genomics](https://www.10xgenomics.com/resources/datasets).
- Datasets sorted by [QuKunLab's SpatialBenchmarking repository](https://github.com/QuKunLab/SpatialBenchmarking).
- DLPFC datasets downloaded by the [script](./Data/fetch_libd.R) based on BayesSpace. And converted to the format of `h5ad` by the `sceasy` package.
- Simulations generated from [our simulation pipeline](../data_Simulation_generate.ipynb).

NOTE: The data used in this study is not included in this repository. Please download the data from the above sources and place them in the corresponding folders.


## Build the conda environment

```bash
# using the conda_config.yml file
conda env create -f conda_config.yml
```

install additional R packages

```R
library(devtools)
library(remotes)
lib_loc = '/home/anaconda3/envs/SpaBench/lib/R/library' # path to the conda environment
remotes::install_github("RubD/Giotto",lib=lib_loc)
devtools::install_github('xzhoulab/SPARK',lib=lib_loc)
remotes::install_github('JEFworks-Lab/MERINGUE', build_vignettes = TRUE,lib=lib_loc)
devtools::install_github("Shufeyangyi2015310117/SC.MEB",lib=lib_loc)
```

## Run the evaluation pipeline [SE genes]

```bash
# activate the conda environment
conda activate SpaBench
cd Bench
# run the evaluation pipeline
# modify the python script to change the parameters
python run_bench.py
```
## Run the evaluation pipeline [Imputation cluster]

see jupyter notebooks under this folder.