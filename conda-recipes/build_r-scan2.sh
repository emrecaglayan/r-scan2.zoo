#!/bin/bash

conda-build -c conda-forge -c bioconda -c jluquette r-scan2

anaconda upload  /home/ljl11/miniconda3_for_building/conda-bld/linux-64/r-scan2-1.0-py38r41h3fd9d12_6.tar.bz2


# Build and upload fastGHQuad so that it doesn't rely on r 3.5 (not used)
#conda skeleton cran fastGHQuad
#conda build -c conda-forge -c bioconda r-fastghquad/ --python 3.10

anaconda login # follow the prompt
anaconda upload --user ecaglayan \
    /home/ch253794/workdir/PROGRAMS/miniconda3_py310/envs/build_env/conda-bld/linux-64/r-fastghquad-1.0.1-r42hbad8821_0.conda



conda-build -c conda-forge -c bioconda -c ecaglayan conda-recipes/r-scan2.zoo

anaconda login # follow the prompt
anaconda upload --user ecaglayan \
    /home/ch253794/workdir/PROGRAMS/miniconda3_py310/envs/build_env/conda-bld/linux-64/r-fastghquad-1.0.1-r42hbad8821_0.conda

