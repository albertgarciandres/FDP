# _TEmerics_

_TEmerics_ is a [Snakemake][snakemake] workflow for the identification and
analysis of chimeric gene-TEs transcripts.

## Table of contents

1. [Installation](#installation)
    - [Cloning the repository](#cloning-the-repository)
    - [Dependencies](#dependencies)
    - [Setting up the virtual enviroment](#setting-up-the-virtual-enviroment)
2. [Usage](#usage)
    - [Preparing inputs](#preparing-inputs)
    - [Running the workflow](#running-the-workflow)
    - [Expected output files](#expected-output-files)
3. [Workflow description](#workflow-description)


## Installation 

The workflow lives inside this repository and will be available for you to run
after the following installation instructions displayed in this section.

### Cloning the repository

Navigate to the desired path on your file system, then clone the repository 
with the following commands:

```bash
git clone /TEmerics_path/
cd TEmerics
```

### Dependencies

For improved reproducibility and reusability of the workflow, as well as an
easy means to run it on a high performance computing (HPC) cluster managed,
e.g., by [Slurm][slurm], all steps of the workflow run inside isolated
[Conda][conda] environments. As a consequence, running this workflow has only 
a few individual dependencies. These are managed by the package manager Conda, 
which needs to be installed on your system before proceeding.

### Setting up the virtual environment

Create and activate the environment with necessary dependencies with Conda:

```bash
conda env create -f environment.yml
conda activate TEmerics
```





[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>

