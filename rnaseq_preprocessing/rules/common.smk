def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])
    
def get_fastq(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(**wildcards):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        return units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        #return [ f"{u.fq1}", f"{u.fq2}" ]

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand(
                "/output/trimmed/{sample}.{unit}.R{group}_val_{group}.fq.gz", group=[1, 2], **wildcards
            )
        # single end sample
        return "/output/trimmed/{sample}.{unit}_trimmed.fq.gz".format(**wildcards)


def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]

