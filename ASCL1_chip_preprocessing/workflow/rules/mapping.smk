rule bt2_align:
    input:
        sample=get_map_reads_input
    output:
        temp("/output/results/mapped/{basename}_{unit}.bam")
    log:
        "/output/logs/bowtie2/{basename}_{unit}.log"
    params:
        index=config["resources"]["ref"]["index"],  # path of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "0.71.1/bio/bowtie2/align"

rule merge_bams:
    input:
        lambda w: expand("/output/results/mapped/{basename}_{unit}.bam",
            unit = units.loc[units['basename'] == w.basename].unit.to_list(),
            basename = w.basename
        )
    output:
        temp("/output/results/merged/{basename}.bam")
    log:
        "/work/logs/picard/mergebamfiles/{basename}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT TMP_DIR=/temp/merge"
    wrapper:
        "0.64.0/bio/picard/mergesamfiles"


rule mark_merged_duplicates:
    input:
        "/output/results/merged/{basename}.bam"
    output:
        bam=temp("/output/tmp/picard_dedup/{basename}.bam"),
        metrics="/output/results/picard_dedup/{basename}.metrics.txt"
    log:
        "/work/logs/picard/picard_dedup/{basename}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT" #TMP_DIR=/output/tmp/picard_dedup"
    wrapper:
        "0.64.0/bio/picard/markduplicates"
