rule fastqc:
    input:
        get_individual_fastq
    output:
        html="/output/results/fastqc/{basename}_{unit}_R{read}_fastqc.html",
        zip="/output/results/fastqc/{basename}_{unit}_R{read}_fastqc.zip"
    log:
        "/output/logs/fastqc/{basename}_{unit}_R{read}.log"
    shell:
        "fastqc {input} --outdir /output/results/fastqc"

# ToDo: add wrapper again and remove temporary script and env after wrapper release with matplotlib dependency
rule multiqc:
    input:
        get_multiqc_input
    output:
        "/output/results/multiqc/multiqc.html"
    log:
        "/output/logs/multiqc.log"
    wrapper:
        "0.64.0/bio/multiqc"