
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
pip3 install Pillow
pip3 install IPython

# pip3 install pycairo
# we use conda pycairo since pip pycairo does not have complete binary
# distribution
micromamba -y install pycairo -c conda-forge

pip3 install argcomplete
pip3 install openpyxl

pip3 install matplotlib
pip3 install seaborn

pip3 install sgkit[vcf] bio2zarr
pip3 install scikit-learn
pip3 install colorcet
pip3 install fastlmm
pip3 install statsmodels

micromamba -y install python-igraph -c conda-forge -c bioconda -c defaults
micromamba -y install plink plink2 -c conda-forge -c bioconda -c defaults

echo "installing R and moimix"
micromamba -y install rpy2 r-devtools r-biocmanager r-tidyverse r-optparse -c conda-forge -c bioconda
micromamba -y install r-ape -c conda-forge -c bioconda -c defaults
micromamba -y install r-matrixmodels r-mcmcpack r-modeltools r-flexmix -c conda-forge -c bioconda
micromamba -y install bioconductor-seqarray bioconductor-biocparallel bioconductor-biobase bioconductor-seqvartools -c conda-forge -c bioconda

# manual installation of moimix
Rscript -e "BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)"

# manual installation of hmmIBD

# manual installation of selscan2

# manual installation of GEMMA

# manual installation of hapQTL

# manual installation of BLINK

# EOF
