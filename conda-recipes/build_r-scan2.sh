#!/bin/bash

# Activate an environment where conda-build is installed.
conda activate build_env

# Build r-scan2.zoo
conda-build -c conda-forge -c bioconda -c ecaglayan conda-recipes/r-scan2.zoo

anaconda login # follow the prompt
anaconda upload --force \
    ~/workdir/PROGRAMS/miniconda3_py310/envs/build_env/conda-bld/linux-64/r-scan2.zoo-0.1-r44h3fd9d12_0.conda

