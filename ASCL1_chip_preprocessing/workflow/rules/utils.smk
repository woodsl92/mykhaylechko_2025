rule samtools_index:
    input:
        "/output/results/filtered/{basename}.sorted.bam"
    output:
        "/output/results/filtered/{basename}.sorted.bam.bai"
    params:
        "" # optional params string
    log:
        "/output/logs/samtools-index/filtered/{basename}.log"
    wrapper:
        "0.64.0/bio/samtools/index"

# rule samtools_index_tmp:
#     input:
#         "/output/tmp/{step}/{basename}.bam"
#     output:
#         "/output/tmp/{step}/{basename}.bam.bai"
#     params:
#         "" # optional params string
#     log:
#         "/output/logs/samtools-index/{step}/{basename}.log"
#     wrapper:
#         "0.64.0/bio/samtools/index"

# rule samtools_index_unit:
#     input:
#         "/output/results/mapped/{basename}.{unit}.bam"
#     output:
#         "/output/results/mapped/{basename}.{unit}.bam.bai"
#     params:
#         "" # optional params string
#     log:
#         "/output/logs/samtools-index/mapped/{basename}.{unit}.log"
#     wrapper:
#         "0.64.0/bio/samtools/index"
