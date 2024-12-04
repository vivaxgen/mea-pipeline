#!/usr/bin/bash

# installation script for vivaxgen ngs-pipeline [https://github.com/vivaxgen/ngs-pipeline]

# optional variable:
# - BASEDIR
# - OMIT

set -eu

# run the base.sh
# Detect the shell from which the script was called
parent=$(ps -o comm $PPID |tail -1)
parent=${parent#-}  # remove the leading dash that login shells have
case "$parent" in
  # shells supported by `micromamba shell init`
  bash|fish|xonsh|zsh)
    shell=$parent
    ;;
  *)
    # use the login shell (basename of $SHELL) as a fallback
    shell=${SHELL##*/}
    ;;
esac

# Parsing arguments
if [ -t 0 ] && [ -z "${BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./vvg-meapl] "
  read BASEDIR
fi

# default value
BASEDIR="${BASEDIR:-./vvg-meapl}"

uMAMBA_ENVNAME="${uMAMBA_ENVNAME:-mea-pl}"
PYVER="${PYVER:-3.12}"
SNAKEMAKEVER="${SNAKEMAKEVER:-9}"
source <(curl -L https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh)

echo "Installing latest htslib tools"
micromamba -y install "bcftools>=1.18" "samtools>=1.18" -c conda-forge -c bioconda -c defaults

echo "Installing vcftools"
micromamba -y install vcftools -c conda-forge -c bioconda

echo "Installing required Python modules"
pip3 install "snakemake<${SNAKEMAKEVER}"
pip3 install snakemake-executor-plugin-cluster-generic
pip3 install cyvcf2
pip3 install pysam
pip3 install pandas
pip3 install sgkit bio2zarr
pip3 install Pillow
pip3 install IPython
pip3 install matplotlib
pip3 install numba
pip3 install seaborn
pip3 install scipy

# pip3 install pycairo
# we use conda pycairo since pip pycairo does not have complete binary
# distribution
micromamba -y install pycairo -c conda-forge

pip3 install argcomplete
pip3 install openpyxl

echo "Cloning vivaxGEN MEA-Pipeline"
git clone https://github.com/vivaxgen/mea-pipeline.git ${ENVS_DIR}/mea-pipeline
ln -sr ${ENVS_DIR}/mea-pipeline/etc/bashrc.d/15-mea-pipeline ${BASHRC_DIR}/
ln -sr ${ENVS_DIR}/mea-pipeline/etc/bashrc.d/95-prompt-history ${BASHRC_DIR}/

echo "installing R and moimix"
micromamba -y install rpy2 r-devtools r-biocmanager r-tidyverse -c conda-forge -c bioconda
micromamba -y install r-matrixmodels r-mcmcpack r-modeltools r-flexmix -c conda-forge -c bioconda
micromamba -y install bioconductor-seqarray bioconductor-biocparallel bioconductor-biobase bioconductor-seqvartools -c conda-forge -c bioconda

R --no-save << EOF
BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
EOF

echo
echo "vivaxGEN MEA-Pipeline has been successfully installed. "
echo "To activate the MEA-Pipeline environment, execute the following command:"
echo
echo "    `realpath ${BINDIR}/activate`"
echo
echo "or source the activation file (eg. inside a script)"
echo
echo "    source `realpath ${BINDIR}/activate`"
echo

# EOF
