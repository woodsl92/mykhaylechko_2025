rule samtools_view:
    input:
        "/output/tmp/picard_dedup/{basename}.bam"
    output:
        temp("/output/tmp/sam-view/{basename}.bam")
    params:
        # TODO: move to config.yaml, including the explanatory comments below
        "-b -F 0x004 -G 0x009 -f 0x001"
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params. -F 1028 removes unmapped reads and duplicates
        # -f 0x001 retains paired reads
        # -G 0x009 removes reads where read is paired but mate is unmapped
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # please provide a respective bed file via "-L path/to/regions.bed"
    log:
        "/output/logs/samtools-view/{basename}.log"
    wrapper:
        "0.64.0/bio/samtools/view"


rule bamtools_filter_json:
    input:
        "/output/tmp/sam-view/{basename}.bam"
    output:
        temp("/output/results/filtered/{basename}.bam")
    params:
          # filters mismatches in all reads and filters pe-reads within a size range given in json-file
        json="/work/config/bamtools_filtering_rules.json"
    log:
        "/output/logs/filtered/{basename}.log"
    wrapper:
        "0.64.0/bio/bamtools/filter_json"

rule samtools_sort_merged:
    input:
        "/output/results/filtered/{basename}.bam"
    output:
        "/output/results/filtered/{basename}.sorted.bam"
    params:
        ""
    log:
        "/output/logs/samtools-sort-filtered/{basename}.log"
    threads:  # Samtools takes additional threads through its option -@
        8
    wrapper:
        "0.64.0/bio/samtools/sort"

rule merge_rep_bams:
    input:
        lambda wildcards: expand("/output/results/filtered/{sample}.sorted.bam", sample=groups[wildcards.group])
    output:
        temp("/output/results/merged_bam/{group}.bam")
    log:
        "/output/logs/merged_bam/{group}.log"
    params:
        "" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.78.0/bio/samtools/merge"

rule subsample_merged_bam:
    input:
        "/output/results/merged_bam/{group}.bam"
    output:
        "/output/results/merged_bam/{group}.subsample.bam"
    log:
        "/output/logs/merged_bam_subsample/{group}.log"
    shell:
        # Subsample to 160m reads
        """
        frac=$( samtools idxstats {input} | cut -f3 | awk 'BEGIN {{total=0}} {{total += $1}} END {{frac=160000000/total; if (frac > 1) {{print "1.0"}} else {{print frac}}}}' )
        samtools view -bs $frac {input} > {output}
        """

rule samtools_sort_merged_bam:
    input:
        "/output/results/merged_bam/{group}.subsample.bam"
    output:
        "/output/results/merged_bam/{group}.subsample.sorted.bam"
    params:
        ""
    log:
        "/output/logs/samtools-sort-filtered/{group}.log"
    threads:  # Samtools takes additional threads through its option -@
        8
    wrapper:
        "0.64.0/bio/samtools/sort"

rule index_merged_subsample_bam:
    input:
        "/output/results/merged_bam/{group}.subsample.sorted.bam"
    output:
        "/output/results/merged_bam/{group}.subsample.sorted.bam.bai"
    log:
        "/output/logs/merged_bam_index/{group}.subsample.log"
    wrapper:
        "0.64.0/bio/samtools/index"
