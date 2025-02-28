## Table of Contents
1. [Run on HPC](#run-on-hpc)
2. [Run Locally](#run-locally)
3. [Pipeline Details](#pipeline-details)
    - [Profiles](#profiles)
    - [Configuration Files](#configuration-files)
    - [SLURM Configuration](#slurm-configuration)

## Run on HPC

This section provides a **quick setup** guide for running the pipeline on **HPC** with test data.

#### **1. Set Up Working Directory**
Navigate to your scratch space and create a folder for the pipeline:
```sh
cd $LSFM_CLUSTER_SCRATCH_USER_PATH
mkdir bio4beginners_nextflow
cd bio4beginners_nextflow
```
#### **2. Load Required Modules**
Load the necessary HPC modules:
```sh
module load USS/2022 gcc/9.4.0-pe5.34 miniconda3/4.12.0
```
If you are using conda for the first time, you need to initialize it:
```sh
conda init bash
```
After running this command, **close and reopen your terminal**. After that you will need to start from step 1.

#### **3. Create Conda Environment**
First, set up Bioconda according to the Bioconda documentation, notably setting up channels:
```sh
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now, create the conda environment:
```sh
conda create --name env_nf nextflow
```

And activate it:
```sh
conda activate env_nf
```

#### **4. Copy Source Code**
Copy the pipeline source folder and prepare the output directory:
```sh
cp -r /cfs/earth/scratch/shared/bioinfo4beginners/Genomics/bio4beginners_nextflow/source .
cd source
mkdir results
```

#### **5. Run the Pipeline**
Execute the Nextflow pipeline using test data:
```sh
nextflow run main.nf -profile conda -c slurm.config \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```
This command:
- Uses the **conda profile** to manage dependencies.
- Runs with **SLURM** using `slurm.config`.
- Converts all file paths to **absolute paths** to avoid path issues.

After execution, results will be stored in the `results/` directory.

## Run Locally

This section provides instructions for running the pipeline on a local machine.

#### **1. Download Source Code and Input Files**

Download .zip file from [google drive](https://drive.google.com/file/d/1A3vcVaQmviO27aJFSsh2r4mkBbFXPKTY/view?usp=sharing) and navigate to the project directory.

#### **2. Create Conda Environment**

First, set up Bioconda according to the Bioconda documentation:

```sh
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now, create the conda environment:

```sh
conda create --name env_nf nextflow
```

and activate it:

```sh
conda activate env_nf
```

#### **3. Run the Pipeline Locally**

Navigate to the source directory and execute the pipeline:

```sh
nextflow run main_local.nf -profile conda \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```

Pipeline will first create new conda enviroment from `Env_Genomics.yml` and it could take around 10 minutes.
After that it will execute in about 5 minutes. 

After execution, results will be stored in the `results/` directory.

#### **4. Run the Pipeline with docker**

If Conda fails, you can run it with any container, such as **Docker, Podman, or Singularity**. Just change the profile to `-profile docker`.  
If your container software requires **root privileges**, then you will also need to run Nextflow with elevated privileges:  

```sh
nextflow run main_local.nf -profile docker \
    --input $(realpath input/reads.fastq)  \
    --reference $(realpath input/chrM.fa) \
    --outdir $(realpath results)
```

## Pipeline Details

![Pipeline Diagram](img/mermaid-diagram-2025-02-13-145840.png)

### Profiles
Profiles describe the environment Nextflow will use to run processes, specified by the `-profile <conda>`. In this pipeline, `conda, docker, podman, singularity` profiles are supported.

The conda environment is described with key word `conda` in each process in `main.nf`. It could be set in two ways:
1. Using an existing conda environment. For example right now it set as:
   ```
   /cfs/earth/scratch/shared/bioinfo4beginners/Genomics/Env_Genomics
   ```
2. Using a `.yml` file to create a new environment or specifying a different conda path.

To modify the environment, update the conda path in `main.nf`.

### Configuration Files

Each Nextflow pipeline has a default configuration file called `nextflow.config`. This file allows defining parameters to be used in `main.nf`.

All input parameters (e.g., input files, reference genome, output directory) can be specified in `nextflow.config` using `params`. These parameters are accessible in `main.nf` through `params.<param_name>`.

For example, instead of specifying inputs in the command line, they can be provided in `nextflow.config`:

```nextflow
params {
    input = "input/reads.fastq"
    reference = "input/chrM.fa"
    outdir = "results"
}
```

This allows running the pipeline without explicitly passing parameters in the command line:

```sh
nextflow run main.nf -profile conda
```

If additional configuration files are needed (e.g., for running the pipeline on an HPC system), they can be included using the `-c` parameter when executing Nextflow.

For example:

```sh
nextflow run main.nf -profile conda -c slurm.config
```

This allows different configurations to be applied without modifying the default `nextflow.config` file.


### SLURM Configuration
The `slurm.config` file is used to run Nextflow with SLURM. It:
- Submits each process as a SLURM job
- Defines resource limits (CPU, memory, partition, etc.)
- Allows customization per process

Typically, to run a script on an HPC system with SLURM, we use `sbatch`. However, Nextflow supports SLURM as an executor and can submit jobs automatically. To enable this, set `process.executor` to `slurm` in the configuration file.

For Earth HPC, it is also required to specify the partition and resource limits. These details can be configured in `slurm.config`. It is also possible to set specific resource requirements for individual processes within this file.



