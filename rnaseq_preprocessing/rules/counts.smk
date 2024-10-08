rule count_matrix:
    input:
        expand(
            "/output/star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "/output/counts/raw_counts.tsv",
    log:
        "/output/logs/count-matrix.log",
    params:
        samples=units["sample"].tolist(),
        strand=get_strandedness(units),
    script:
        "../scripts/count-matrix.py"
        