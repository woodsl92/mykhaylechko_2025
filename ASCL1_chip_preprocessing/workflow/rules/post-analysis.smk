rule preseq_lc_extrap:
    input:
        "/output/results/filtered/{basename}.bam"
    output:
        "/output/results/preseq/{basename}.lc_extrap"
    params:
        "-v -seed 1"
    log:
        "/output/logs/preseq/{basename}.log"
    wrapper:
        "0.64.0/bio/preseq/lc_extrap"

rule genomecov:
    input:
        bam="/output/results/filtered/{basename}.sorted.bam",
        flag="/output/results/filtered/{basename}.flagstat"
    output:
        "/output/results/bed_graph/{basename}.bedgraph"
    log:
        "/output/logs/bed_graph/{basename}.log"
    shell: #-fs option was not used because there are no single end reads any more # scale by library size
        "bedtools genomecov -bg -pc -ibam {input.bam} -scale $(grep 'mapped (' /output/results/filtered/{wildcards.basename}.flagstat | awk '{{print 1000000/$1}}') > {output}"


rule sort_genomecov:
    input:
        "/output/results/bed_graph/{basename}.bedgraph"
    output:
        temp("/output/results/bed_graph/{basename}.sorted.bedgraph")
    log:
        "/output/logs/sort_genomecov/{basename}.log"
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule bedGraphToBigWig:
    input:
        bedGraph="/output/results/bed_graph/{basename}.sorted.bedgraph",
        chromsizes="/work/resources/ref/genome.chrom.sizes"
    output:
        "/output/results/big_wig/{basename}.bigWig"
    log:
        "/output/logs/big_wig/{basename}.log"
    # params:
    #     ""
    # wrapper:
    #     "0.64.0/bio/ucsc/bedGraphToBigWig"
    shell:
        "bedGraphToBigWig {input.bedGraph} {input.chromsizes} {output}"

rule create_igv_bigwig:
    input:
        "/work/resources/ref/genome.bed",
        expand("/output/results/big_wig/{basename}.bigWig", basename = basenames.index)
    output:
        "/output/results/IGV/big_wig/merged_library.bigWig.igv.txt"
    log:
        "/output/logs/igv/create_igv_bigwig.log"
    shell:
        "find {input} -type f -name '*.bigWig' -exec echo -e '/output/results/IGV/big_wig/\"{{}}\"\t0,0,178' \;  > {output} 2> {log}"

rule compute_matrix:
    input:
         bed="/work/resources/ref/genome.bed",
         bigwig=expand("/output/results/big_wig/{basename}.bigWig", basename = basenames.index)
    output:
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=temp("/output/results/deeptools/matrix_files/matrix.gz"),
        matrix_tab=temp("/output/results/deeptools/matrix_files/matrix.tab")
    log:
        "/output/logs/deeptools/compute_matrix.log"
    threads: 20
    params:
        command="scale-regions",
        extra="--regionBodyLength 1000 "
              "--beforeRegionStartLength 3000 "
              "--afterRegionStartLength 3000 "
              "--skipZeros "
              "--smartLabels "
              "--numberOfProcessors 20 "
    wrapper:
        "0.64.0/bio/deeptools/computematrix"

rule plot_profile:
    input:
         "/output/results/deeptools/matrix_files/matrix.gz"
    output: #ToDo: add description to report caption
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotprofile.html.
        plot_img=report("/output/results/deeptools/plot_profile.pdf", caption="/work/report/plot_profile_deeptools.rst", category="GenomicRegions"),
        data=temp("/output/results/deeptools/plot_profile_data.tab")
    log:
        "/output/logs/deeptools/plot_profile.log"
    params:
        ""
    wrapper:
        "0.64.0/bio/deeptools/plotprofile"

rule plot_heatmap:
    input:
         "/output/results/deeptools/matrix_files/matrix.gz"
    output:  #ToDo: add description to report caption
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotheatmap.html.
        heatmap_img=report("/output/results/deeptools/heatmap.pdf", caption="/work/report/plot_heatmap_deeptools.rst", category="Heatmaps"),
        heatmap_matrix=temp("/output/results/deeptools/heatmap_matrix.tab")
    log:
        "/output/logs/deeptools/heatmap.log"
    params:
        ""
    wrapper:
        "0.64.0/bio/deeptools/plotheatmap"

rule phantompeakqualtools:
    input:
         "/output/results/filtered/{basename}.sorted.bam"
    output:  #ToDo: add description to report caption
        res_phantom="/output/results/phantompeakqualtools/{basename}.phantompeak.spp.out",
        r_data="/output/results/phantompeakqualtools/{basename}.phantompeak.Rdata",
        plot=report("/output/results/phantompeakqualtools/{basename}.phantompeak.pdf", caption="report/plot_phantompeak_phantompeakqualtools.rst", category="Phantompeak")
    threads:
        8
    log:
        "/output/logs/phantompeakqualtools/{basename}.phantompeak.log"
    shell:
        "( Rscript -e \"library(caTools); source('/work/workflow/scripts/run_spp.R')\" "
        "  -c={input} -savp={output.plot} -savd={output.r_data} "
        "  -out={output.res_phantom} -p={threads} 2>&1 ) >{log}"

rule phantompeak_correlation:
    input:
        data="/output/results/phantompeakqualtools/{basename}.phantompeak.Rdata",
        header="/work/workflow/header/spp_corr_header.txt"
    output:
        "/output/results/phantompeakqualtools/{basename}.spp_correlation_mqc.tsv"
    log:
        "/output/logs/phantompeakqualtools/correlation/{basename}.spp_corr.log"
    script:
        "/work/workflow/scripts/phantompeak_correlation.R"

rule phantompeak_multiqc:
    # NSC (Normalized strand cross-correlation) and RSC (relative strand cross-correlation) metrics use cross-correlation
    input:
        data="/output/results/phantompeakqualtools/{basename}.phantompeak.spp.out",
        nsc_header="/work/workflow/header/nsc_header.txt",
        rsc_header="/work/workflow/header/rsc_header.txt"
    output:
        nsc="/output/results/phantompeakqualtools/{basename}.spp_nsc_mqc.tsv",
        rsc="/output/results/phantompeakqualtools/{basename}.spp_rsc_mqc.tsv"
    log:
        "/output/logs/phantompeakqualtools/correlation/{basename}.nsc_rsc.log"
    shell:
        "( gawk -v OFS='\t' '{{print $1, $9}}' {input.data} | cat {input.nsc_header} - > {output.nsc} && "
        "  gawk -v OFS='\t' '{{print $1, $10}}' {input.data} | cat {input.rsc_header} - > {output.rsc} 2>&1 ) >{log}"