# SafePathsRNAPC
This repository contains:
- The implementation of the computation of Safe Paths under different models of Path Cover in a Directed Acyclic Graph (DAG). C++ code.
- Experimental evaluation of the algorithm by computing contigs in a RNA Transcript Assembly problem. Jupyter Notebook.

A comprehensive explanation of both parts can be found in [Insert reference to publication here].

## C++ algorithm
First clone the repo:
```
 git clone https://github.com/elarielcl/SafePathsRNAPC.git
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

## Jupyter Notebooks
After compiling the C++ code you can replicate our experiments by running the Jupyter Notebooks in the folder `data`. These notebooks create intermediate files in the different subfolders of `data`. The notebooks are self-contained and must be run in the following order (indicated in the notebooks too):

- `data_manipulation/graph_creation.ipynb`
- `experiments/run_experiments.ipynb`
- `evaluation/compute_metrics.ipynb`
- `evaluation/compute_tables.ipynb`


 # Contact
 Any error, improvement or suggestion you can write me to `elarielcl` at Gmail.
