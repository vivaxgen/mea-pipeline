VCF Processing
==============

This section describes how to process VCF files using the MEA-Pipeline.


Query Notation
--------------

.. list-table:: Query Notation
    :widths: 30 70
    :header-rows: 1
    
    * - Notation
      - Remark
    * - ``-variants``
      - Select the variants based on a bed file or a position file that is
        provided in the config.yaml with parameter ``variant_file``
    * - ``-samples``
      - Select the samples based on a headerless text file that is provided
        in the config.yaml with parameter ``sample_file``
    * - ``-d[INT]``
      - Set the minimum depth to [INT] for variant calls.
    * - ``-HETd[INT]r[FLOAT]``
      - Set alleles to het calls if the minimum alternate read count > INT and
        minimum alternate read ratio (out of total reads) > FLOAT.
    * - ``-HET2[OPT]``
      - Set het alleles to OPT, where OPT can be REF (reference alleles), ALT
        (alternate alleles), MISS (missing alleles) or MAJ (major alleles).
    * - ``-ANN``
      - Annotate variants using snpEff software, with provided snpEff database
        set in the config.yaml with parameter ``snpEff_config_file``,
        ``snpEff_data_dir``, and ``snpEff_db``.
    * - ``-BCSQ``
      - Annotate using  the ``samtools bcsq`` command, with gff file set in
        parameter ``gff_file``.
    * - ``-V[FLOAT]``
      - Select the variants that have sample missingness up to FLOAT value.
    * - ``-S[FLOAT]``
      - Select the samples that have variant missingness up to FLOAT value.
    * - ``-MAC[INT]``
      - Select variants with minimum allele count (MAC) of INT.
    * - ``-MAF[FLOAT]``
      - Select variants with minimum allele frequency (MAF) of FLOAT.
    * - ``-atom``
      - Decompose complex variants using the ``samtools norm`` command.
    * - ``-split``
      - Split multi-allelec variants into individual variants.
    * - ``-snv``
      - Select single-nucleotide variants.
    * - ``-dedup``
      - Deduplicate variants that have identical position by selecting those
        with the lowest sample missingness.
    * - ``-FWS[FLOAT]``
      - Select samples with Fws > FLOAT, as calculated by moimix.

