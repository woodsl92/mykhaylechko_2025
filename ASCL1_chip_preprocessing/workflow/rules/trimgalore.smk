rule trimgalore_pe:
    input:
        get_fastqs
    output:
        fastq=["/output/results/trimmed/{basename}_{unit}_R1_val_1.fq.gz", "/output/results/trimmed/{basename}_{unit}_R2_val_2.fq.gz"],
        report=["/output/results/trimmed/{basename}_{unit}_R1.fq.gz_trimming_report.txt", "/output/results/trimmed/{basename}_{unit}_R2.fq.gz_trimming_report.txt"]
    log:
        "/output/logs/trim_galore/{basename}_{unit}_pe.log"
    shell:
        "trim_galore -q 20 -fastqc --paired -o /output/results/trimmed/ {input}"

rule trimgalore_se:
    input:
        get_fastqs
    output:
        fastq="/output/results/trimmed/{basename}_{unit}_trimmed.fq.gz",
        report="/output/results/trimmed/{basename}_{unit}.fq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "/output/logs/trim_galore/{basename}_{unit}.log"
    wrapper:
        "0.71.1/bio/trim_galore/se"
