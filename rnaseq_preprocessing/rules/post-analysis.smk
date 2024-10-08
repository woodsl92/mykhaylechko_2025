rule samtools_flagstat:
    input:
        "/output/star/{sample}-{unit}/Aligned.out.bam"
    output:
        "/output/star/{sample}-{unit}.flagstat"
    log:
        "/output/logs/samtools-flagstat/{sample}-{unit}.log"
    wrapper:
        "0.64.0/bio/samtools/flagstat"

        