
vivaxGEN MEA-Pipeline Documentation
===================================

The MEA-Pipeline provides set of tools and Snakemake workflow for VCF processing,
population genetics and GWAS analysis.


Installation
------------

To install, execute the following command::

	"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/mea-pipeline/main/install.sh)


Examples
--------

The following are examples on how to evoke the pipeline for different analyses.

For the all following commands, it is assumed that the initial VCF files reside
in the ``base`` directory.

To filter VCF files for variant with sample missingness of 0.5, then variant
with MAF (minor allele frequency) of 0.1, and then filter for samples with
variant missingness of 0.5::

   $ mea-pl run-vcf-query -q base-V0.5-MAF0.1-S0.5

The following command filters for biallelic SNPs before appying the above
criteria::

   $ mea-pl run-vcf-query -q base-split-snv-dedup-V05-MAF0.1-S0.5

The following set the minimum depth to 10 and set the het calls if the
count and ratio (out of total) of the alternate allele reads is 5 and 0.1,
respectively, before applying the above criteria::

   $ mea-pl run-vcf-query -q base-HETd5r0.1-d10-split-snv-dedup-V0.5-MAF0.1-S0.5


Configuration
-------------

MEA-Pipeline requires a `config.yaml` file in the current working directory to run.
A minimal configuration file, using PvP01_v2 reference genome, is shown here:

.. code-block:: console

   regions:
     - PvP01_01_v2
     - PvP01_02_v2
     - PvP01_03_v2
     - PvP01_04_v2
     - PvP01_05_v2
     - PvP01_06_v2
     - PvP01_07_v2
     - PvP01_08_v2
     - PvP01_09_v2
     - PvP01_10_v2
     - PvP01_11_v2
     - PvP01_12_v2
     - PvP01_13_v2
     - PvP01_14_v2


Updating the Pipeline
---------------------

To update the pipeline, use the following command::

   $ $VVGBIN/update-box

The above command will update the pipeline and all necessary dependencies.


.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   userdocs/intro.rst
   userdocs/vcf_processing.rst
   userdocs/popgen_analysis.rst
   userdocs/gwas_analysis.rst


.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   develdocs/implementation.rst

