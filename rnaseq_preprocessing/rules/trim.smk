rule trimgalore_pe:
    input:
        get_fastq,
    output:
        fastq=["/output/trimmed/{sample}.{unit}.R1_val_1.fq.gz", "/output/trimmed/{sample}.{unit}.R2_val_2.fq.gz"],
        report=["/output/trimmed/{sample}.{unit}.R1.fq.gz_trimming_report.txt", "/output/trimmed/{sample}.{unit}.R2.fq.gz_trimming_report.txt"],
        zip = ["/output/trimmed/{sample}.{unit}.R1_val_1_fastqc.zip", "/output/trimmed/{sample}.{unit}.R2_val_2_fastqc.zip"]
    log:
        ["/output/logs/trim_galore/{sample}.{unit}_R1.log","/output/logs/trim_galore/{sample}.{unit}_R2.log"]
    shell:
        "trim_galore -q 20 -fastqc --paired -o /output/trimmed/ {input}"

rule trimgalore_se:
    input:
        get_fastq
    output:
        fastq="/output/trimmed/{sample}.{unit}_trimmed.fq.gz",
        report="/output/trimmed/{sample}.{unit}.fq.gz_trimming_report.txt"
    log:
        "/output/logs/trim_galore/{sample}_{unit}.log"
    shell:
        "trim_galore -q 20 -fastqc -o /output/trimmed/ {input}"

# rule cutadapt_pe:
#     input:
#         get_fastq,
#     output:
#         fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
#         fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
#         qc="trimmed/{sample}-{unit}.qc.txt",
#     params:
#         "-a {} {}".format(
#             config["trimming"]["adapter"], config["params"]["cutadapt-pe"]
#         ),
#     log:
#         "logs/cutadapt/{sample}-{unit}.log",
#     wrapper:
#         "0.17.4/bio/cutadapt/pe"


# rule cutadapt:
    # input:
    #     get_fastq,
    # output:
    #     fastq="trimmed/{sample}-{unit}.fastq.gz",
    #     qc="trimmed/{sample}-{unit}.qc.txt",
    # params:
    #     "-a {} {}".format(
    #         config["trimming"]["adapter"], config["params"]["cutadapt-se"]
    #     ),
    # log:
    #     "logs/cutadapt/{sample}-{unit}.log",
    # wrapper:
    #     "0.17.4/bio/cutadapt/se"
