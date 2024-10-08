rule get_genome_test_data:
    output:
        "/work/resources/ref/genome.fasta"
    log:
        "/work/logs/ref/get-genome.log"
    params:
        build=config["resources"]["ref"]["build"]
    cache: True
    shell:
        "curl -L https://hgdownload.soe.ucsc.edu/goldenPath/{params.build}/bigZips/{params.build}.fa.gz |gzip -d > {output} 2> {log}"


rule get_annotation:
    output:
        "/work/resources/ref/annotation.gtf"
    params:
        build=config["resources"]["ref"]["build"]
    log:
        "/output/logs/ref/get_annotation.log"
    cache: True  # save space and time with between workflow caching (see docs)
    shell:
        "curl -L https://hgdownload.soe.ucsc.edu/goldenPath/{params.build}/bigZips/genes/{params.build}.refGene.gtf.gz |gzip -d > {output} 2> {log}"

rule gtf2bed:
    input:
        "/work/resources/ref/annotation.gtf"
    output:
        "/work/resources/ref/genome.bed"
    log:
        "/output/logs/ref/gtf2bed.log"
    shell:
        "/work/workflow/scripts/gtf2bed {input} > {output} 2> {log}"

rule genome_faidx:
    input:
        "/work/resources/ref/genome.fasta"
    output:
        "/work/resources/ref/genome.fasta.fai"
    log:
        "/output/logs/ref/genome-faidx.log"
    cache: True
    wrapper:
        "0.64.0/bio/samtools/faidx"

rule chromosome_size:
    input:
        "/work/resources/ref/genome.fasta.fai"
    output:
        "/work/resources/ref/genome.chrom.sizes"
    log:
        "/output/logs/ref/chromosome_size.log"
    shell:
        "cut -f 1,2 {input} > {output} 2> {log}"