# scRiboQuant
Workflow for single-cell ribo-seq analysis


**Author: @vivekbhr**

[![DOI](https://zenodo.org/badge/248052517.svg)](https://zenodo.org/badge/latestdoi/248052517) Bhardwaj V. (2021) scRiBoQuant: Workflow for processing single-cell Ribo-seq data  


## How to run

### 1. Install via conda

Create a new environment and install the workflow in the environment. (recommended)

```
conda create -n riboseq -c bioconda -c vivekbhr scriboquant
```

Alternative: Install the workflow directly in your base environment.

```
conda install -c bioconda -c vivekbhr scriboquant
```

### 1b. (Alternative) Install via pip from github

Create a new environment.

```
conda create -n riboseq python
```
Install the workflow from the master branch

```
conda activate riboseq
pip install git+https://github.com/vivekbhr/scRiboQuant.git
```


### 2. configure the config.yaml

config.yaml file is needed to provide files and arguments on the workflow that you don't want to repeat everytime you run the workflow, such as path to fasta and gtf files. Copy the `config.yaml` file from [here](./scriboquant/config.yaml) and modify as per your own requirements.

Files needed (provide in the config):

1) Path to the (indexed) genome fasta file
2) (optional) Path to bed file to mask certain regions in the genome and fasta files to append certain regions.
3) Cell barcodes (.txt/tsv file, first column must be cell-barcodes, no header)

### 3. Test the workflow

in order to test the workflow we can do a run on the real dataset, but with `--downsample <n_reads> ` command.

For example, to run the workflow on 1000 downsamples reads (provided all files in `config.yaml` are accessible)

```
conda activate riboseq

scRiboQuant -i <testdata_folder> --downsample 1000 -o . -c <your_config.yaml> -j <jobs> -cl
```

### 4. Submission parameters

#### Running on HPC Cluster
  - In the workflow command above, **j** is the number of parallel jobs you want to run, **-cl** means submit to cluster (default is to run locally). Therefore if you wish to run the workflow on a cluster, simply use the workflow with the -cl command on the submission node.

  - cluster configuration, such as memory and cluster submission command are placed in [cluster_config.yaml](./scriboquant/cluster_config.yaml), and can be modified to suite the users internal infrastructure.

#### Dry-run
In order to just test what the workflow would do, use the command `-s ' -np' `


### Other technical Notes

  - After running the pipeline, **LOG** file are stored in the **<output>/log/** directory and the workflow top-level log is in scRiboQuant.log file.

  - Currently the -o option is not very flexible and and pipeline works only when it's executed in the output directory.

  - Use the -t argument to specify a local directory for temp files. Default is to use the /tmp/ folder, which might have low space on cluster (unless tmpspace is specified in cluster_config.yaml)

  - **Manual interruption of the workflow**: Simple Ctrl+C is enough to cancel/inturrupt the workflow. However, in some cases re-running the workflow after inturruption might fail with message "Locked working directory". In that case, please run the workflow with `-s ' --unlock'` once.
  
  
### DAG (Directed Acyclic Graph) of the Workflow

<img src=./dag.png width=500 height=500>

**Note: This workflow is provided as is and feature/bug issues would not be entertained as I am not a member of the scRiboSeq team. Please feel free to copy/modify it for your project with attribution**
