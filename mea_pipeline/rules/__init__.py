def path(snakefile: str):
    return __path__[0] + "/" + snakefile


vcf_query = path("vcf_query.smk")
utils = path("utils.smk")
vcf = path("vcf.smk")
popgen = path("popgen.smk")
selscan = path("selscan.smk")
fastlmm = path("fastlmm.smk")


# EOF
