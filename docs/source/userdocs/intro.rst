Introduction to MEA-Pipeline
=============================

The MEA-Pipeline is a comprehensive framework designed for the analysis of genomic data.
It provides a set of tools and workflows for various tasks, including VCF processing,
population genetics, and GWAS analysis.

VCF Processing
--------------

The VCF processing workflow prioritizes retaining genotype data in VCF format
until the final step, converting to other formats only when required by
specific analysis tools.
Though slower to process, VCF files consolidate key information—genotypes,
allele depths, variant quality, and annotations—into a single, consistent
source.
This enables flexible, accurate processing of tasks like converting heterozygous
to major allele calls or adjusting minimum depth at any step.
Performance is optimized through multicore or multinode parallelization, with
processing divided by chromosome or genomic region as specified in config.yaml.


Population Genetics analysis
----------------------------

MEA-Pipeline provides common population genetics analysis for infectious
disease pathogens by extending the VCF processing to include additional
targets.
Almost most of population genetics analyses require tabular metadata file in a
tab-delimited format or feather format.
The metadata file should have at least a column named SAMPLE.
All column names should avoid space, as it can complicate the text matching.

GWAS Analysis
-------------

MEA-Pipeline support GWAS analysis with wrappers for FaST-LMM, WarpedLMM, GEMMA,
BLINK and hapQTL.
These wrappers are executed using the same query notation system.
VCF files are automatically converted to the input format required by the respected
software.
