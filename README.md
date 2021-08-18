# SafePathsRNAPC
This repository contains:
- The implementation of the computation of Safe Paths under different models of Path Cover in a Directed Acyclic Graph (DAG). C++ code.
- Experimental evaluation of the algorithm by computing contigs in a RNA Transcript Assembly problem. Jupyter Notebook.

A comprehensive explanation of both parts can be found in [Insert reference to publication here].

## C++ algorithm
First clone the repo:
```
 git clone https://github.com/algobio/SafePathsRNAPC.git
 ```
 
This project is a CMake project. To build this project with some runnables you should do

```
cd ../..
mkdir build
cd build
cmake ..
cmake .. # Issue: second cmake necessary to compile external library
make
```
> This C++ project downloads the [LEMON graph library](https://lemon.cs.elte.hu/trac/lemon), which is stored in a [Mercurial](https://www.mercurial-scm.org/) repository. As such, the installation requires Mercurial.


## Jupyter Notebooks
After compiling the C++ code you can replicate our experiments by running the Jupyter Notebooks in the folder `data`. These notebooks create intermediate files in the different subfolders of `data`. The notebooks are self-contained and must be run in the following order (indicated in the notebooks too):

- `data_manipulation/graph_creation.ipynb`
- `experiments/run_experiments.ipynb`
- `evaluation/compute_metrics.ipynb`
- `evaluation/compute_tables.ipynb`

These notebooks correspond to the experiments for `Homo sapiens`. The experiments for other species can be found (following the same structure) in the folders `data/mouse` (Mus musculus), `data/triticum_aestivum`, `data/hordeum_vulgare`, `data/fruit_fly` (Drosophila melanogaster) and `data/magnaporthe_oryzae`.

Once all these notebooks have be run, you can run the notebook `compute_summary_tables.ipynb`.

> The package requirements to run our jupyter notebooks can be found at `data/requirements.txt`

> The datasets used are two `BED` files (per dataset) created from `GFF` files in the [Enssembl project](https://www.ensembl.org/index.html). The script used to transform the GFF into BED files can be found at `data/scripts/gtf2bed.py`

 # Contact
 Any error, improvement or suggestion you can write me to `elarielcl` at Gmail.

