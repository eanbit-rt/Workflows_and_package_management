Using Bioconda for bioinformatics workflows
===

This is based on Andrew Perry's [tutorial](https://github.com/MonashBioinformaticsPlatform/bioconda-tutorial).

# Conda

The `conda` [package manager](https://en.wikipedia.org/wiki/Package_manager) was originally built by ContinuumIO to support managing dependencies for the [Anaconda Python](https://www.anaconda.com/download/#linux) distribution. It has since been adopted by community projects such as [conda forge](https://conda-forge.org/) and [bioconda](https://bioconda.github.io/) as an easy way to distribute software packages. The bioconda [package repository](https://bioconda.github.io/recipes.html) contains over 8000 bioinformatics packages and is also a crucial component of the [biocontainers](https://biocontainers.pro/) project.

## Advantages to using (bio)conda

* Manages self-contained environments, including dependencies. An "ordinary user" can install it for themselves.

* Large ecosystem of precompiled packages, organized as 'channels' (eg conda-forge, bioconda)

* Language agnostic (not only Python !)

* Creating new packages is straightforward compared with many systems (a recipe composed of a short build.sh and meta.yml).

* Interoperates gracefully with Python virtualenvs and pip.

* Bioconda project imposes quality control on packages with [best practices](http://bioconda.github.io/guidelines.html) and automated testing and package builds.

* Open Source tool (BSD), backed by a company (Continuum Analytics). Some compiled packages in the default 'channel' aren't Open Source.

## How to install conda

1. Get Miniconda from https://conda.io/miniconda.html (on Linux):

   ```
   curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```

   or

   ```
   wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```

   On MacOS:

   ```
   curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
   ```
 
2. Run the installer:
   ```
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. At the prompt agree to the license (type *yes* and press Enter) and install to the default location ($HOME/miniconda3). When asked, agree to add the Miniconda `bin` directory to the PATH, then log out and back in again.

4. Test that it is installed:
   ```
   $ conda help
   usage: conda [-h] [-V] command ...

   conda is a tool for managing and deploying applications, environments and packages.

   ...
   ```

## Install bioconda
1. Configure conda channels for bioconda
   ```
   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   ```

2. Check your conda config:
   ```
   $ cat ~/.condarc 
   channels:
     - conda-forge
     - bioconda
     - defaults
   ```

## Search packages

You can discover if a tool you are interested in is available as a conda package by executing `conda search`. For example, if you want to look for the seqtk tool, you can run:
```
$ conda search seqtk
```

Alternatively, you can search for it on the https://anaconda.org website, which is often faster.

## Install and run seqtk

1. Install the seqtk package:
```
$ conda install seqtk
Fetching package metadata ...............
Solving package specifications: .

Package plan for installation in environment /Users/perry/miniconda3:

The following NEW packages will be INSTALLED:

    seqtk:     1.2-0            bioconda

The following packages will be UPDATED:

    conda:     4.3.31-py36_0             --> 4.3.33-py36_0 conda-forge

The following packages will be SUPERSEDED by a higher-priority channel:

    conda-env: 2.6.0-h36134e3_0          --> 2.6.0-0       conda-forge

conda-env-2.6. 100% |################################| Time: 0:00:00 395.33 kB/s
seqtk-1.2-0.ta 100% |################################| Time: 0:00:00  63.60 kB/s
conda-4.3.33-p 100% |################################| Time: 0:00:01 279.73 kB/s
```

2. Download some sample data
   ```
   wget -c -O A.fastq.gz https://zenodo.org/record/1324070/files/wt_H3K4me3_read1.fastq.gz
   ```

   or

   ```
   curl https://zenodo.org/record/1324070/files/wt_H3K4me3_read1.fastq.gz >A.fastq.gz
   ```

3. Make a name list: `echo 'SRR5680996.11626500' >name.lst`
4. Extract a sequence: `seqtk subseq A.fastq.gz name.lst`
5. Reverse complement the sequence: `seqtk subseq A.fastq.gz name.lst | seqtk seq -r`

## Uninstall and re-install in 'environment'

The examples shown above have worked in the conda "root environment". For reproducible research it is good to keep dependencies needed for a specific task isolated from those needed for other tasks. Conda supports [environments](https://conda.io/docs/user-guide/tasks/manage-environments.html), which are isolated workspaces that can have their own versions of packages installed. There are two key concepts to environments:

1. Environments only store the packages explicitely installed by the user (and their dependencies). What is installed in an environment can be listed with `conda list`.

2. Different environments can have different versions of tools installed. For example one environment can have Python 2 installed whereas another can have Python 3.

Having learned about environments we can now uninstall `seqtk` from the root environment:

* Uninstall seqtk: `conda remove -y seqtk`

and create an environment named `seqtk_env` into which we have installed `seqtk`.

* Install in environment: `conda create -n seqtk_env -y seqtk`

Switching to an environment is called *activating* that environment:

*  Activate environment: `source activate seqtk_env`

And if you want to save or share the environment definition you can export it to a file:

* Export environment specification to a file: `conda env export >seqtk.yml`

Finally, that environment definition can be used to create a new environment in a conda install:

* Create a new environment from an exported specification: `conda env create -n seqtk_again --file seqtk.yaml`
