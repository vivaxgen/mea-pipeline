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
source <(curl -L https://raw.githubusercontent.com/vivaxgen/vvg-base/main/install.sh)

echo "Cloning vivaxGEN MEA-Pipeline"
git clone https://github.com/vivaxgen/mea-pipeline.git ${ENVS_DIR}/mea-pipeline
ln -sr ${ENVS_DIR}/mea-pipeline/etc/bashrc.d/15-mea-pipeline ${BASHRC_DIR}/
ln -sr ${ENVS_DIR}/mea-pipeline/etc/bashrc.d/95-prompt-history ${BASHRC_DIR}/

source ${ENVS_DIR}/mea-pipeline/etc/inst-scripts/inst-deps.sh

echo "mea-pipeline" >> ${ETC_DIR}/installed-repo.txt

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
