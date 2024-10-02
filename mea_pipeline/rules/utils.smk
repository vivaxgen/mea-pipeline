# mea_pipeline/rules/utils.smk
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

MEA_PIPELINE_BIN = config['MEA_PIPELINE_BASE'] + '/bin'

alphanum = r"\w+"
alphanumdash = r"[\w-]+"
alphanumdashdotplus = r"[\.+\w-]+"
alphanumdashdotplushashat = r"[@#\.+\w-]+"
digit = r"\d+"  # integers
digitdot = r"[\.\d]+"  # float, fraction


wildcard_constraints:
    reg=alphanumdash,  # chromosome region name
    complete_region=alphanumdash,
    population=alphanumdashdotplus,
    population_notation=alphanumdashdotplushashat,
    fname=r"[\.+\w-]+",  # standard file name
    MAC=digit,
    MAF=digitdot,
    mindepth=r"\d+",
    st=r"[\.\d]+",
    vt=r"[\.\d]+",
    qual=digitdot,
    fws=digitdot,
    var=r"[\w]+",  # variable name
    grp=alphanum,
    dhets=digit,
    rhets=digitdot,
    fn=alphanumdashdotplus,


# constants

cli = "mea-pl"

# parameters

REGIONS = config["regions"]
complete_region = config["complete_region"]
chromtranslation_file = config["chromtranslation_file"]

mindepth = config["mindepth"]


# utilities

color_palettes = {
    "xgfs_normal12": [
        "#ebac23",
        "#b80058",
        "#008cf9",
        "#006e00",
        "#00bbad",
        "#d163e6",
        "#b24502",
        "#ff9287",
        "#5954d6",
        "#00c6f8",
        "#878500",
        "#00a76c",
        "#bdbdbd"
      ]
}
# for color palette, see https://gist.github.com/xgfs/37436865b6616eebd09146007fea6c09


rule index_csi:
    # indexing VCF with .csi
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz",
    output:
        "{pfx}/{reg}.vcf.gz.csi",
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    # indexing VCF with .tbi
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz",
    output:
        "{pfx}/{reg}.vcf.gz.tbi",
    shell:
        "bcftools index --tbi {input}"


rule vcf2zarr:
    # convert VCF to Zarr
    threads: 3
    input:
        vcf = "{pfx}/{reg}.vcf.gz",
        idx = "{pfx}/{reg}.vcf.gz.tbi",
    output:
        zarr = directory("{pfx}/{reg}.vcz")
    shell:
        "vcf2zarr convert --no-progress --worker-processes {threads}"
        " {input.vcf} {output.zarr}"


# some litle functions

def count_file_lines(infile):
    c = 0
    for line in open(infile):
        if line.strip():
            c += 1
    return c


# EOF
