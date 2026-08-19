#!/bin/bash

# Before doing this, upload your local repo to your github's repo
# Then specify your github's repo in your meta.yaml file
# Then create a token from github because this command will prompt you for the password
# You may still get errors during the build and will most likely need to resolve them.
cd /home/ch253794/workdir/PROGRAMS/SCAN2_ZOO_R/scan2.zoo

git add conda-recipes & git commit -m "Update"
git push -u origin main

conda-build -c conda-forge -c bioconda -c ecaglayan conda-recipes/r-scan2.zoo > conda-build2.log 2>&1

anaconda login # follow the prompt
anaconda upload --force \
    /home/ch253794/workdir/PROGRAMS/miniconda3_py310/envs/build_env/conda-bld/linux-64/r-scan2.zoo-0.1-py314r44h3fd9d12_0.conda


