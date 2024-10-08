rule samtools_flagstat:
    input:
        "/output/results/filtered/{basename}.sorted.bam"
    output:
        "/output/results/filtered/{basename}.flagstat"
    log:
        "/output/logs/samtools-flagstat/filtered/{basename}.log"
    wrapper:
        "0.64.0/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "/output/results/filtered/{basename}.sorted.bam",
        idx = "/output/results/filtered/{basename}.sorted.bam.bai"
    output:
        "/output/results/filtered/{basename}.filtered.idxstats"
    log:
        "/output/logs/samtools-idxstats/filtered/{basename}.filtered.log"
    wrapper:
        "0.64.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        "/output/results/filtered/{basename}.sorted.bam"
    output:
        "/output/results/filtered/{basename}.mapped.stats.txt"
    params:
        ""
    log:
        "/output/logs/samtools-stats/filtered/{basename}.mapped.log"
    wrapper:
        "0.64.0/bio/samtools/stats"