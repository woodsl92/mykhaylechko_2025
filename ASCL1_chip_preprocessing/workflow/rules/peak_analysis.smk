rule plot_fingerprint:
    input:
        bam_files=["/output/results/filtered/{basename}.sorted.bam", "/output/results/filtered/{control}.sorted.bam"],
        bam_idx=["/output/results/filtered/{basename}.sorted.bam.bai", "/output/results/filtered/{control}.sorted.bam.bai"],
        jsd_sample="/output/results/filtered/{control}.sorted.bam"
    output:  #ToDo: add description to report caption
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint=report("/output/results/deeptools/{basename}-{control}.plot_fingerprint.pdf", caption="report/plot_fingerprint_deeptools.rst", category="QC"),
        counts="/output/results/deeptools/{basename}-{control}.fingerprint_counts.txt",
        qc_metrics="/output/results/deeptools/{basename}-{control}.fingerprint_qcmetrics.txt"
    log:
        "/output/logs/deeptools/plot_fingerprint.{basename}-{control}.log"
    params:
        "--labels {basename} {control}",
        "--skipZeros ",
        "--numberOfSamples 500000 ", # ToDo: to config?
    threads:
        8
    wrapper:
        "0.66.0/bio/deeptools/plotfingerprint"

rule macs2_callpeak_broad:
    input:
        treatment="/output/results/filtered/{basename}.sorted.bam",
        control="/output/results/filtered/{control}.sorted.bam"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("/output/results/merged_macs2_callpeak/{basename}-{control}.broad",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                #  "_treat_pileup.bdg",
                #  "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "/output/logs/macs2/callpeak.{basename}-{control}.broad.log"
    params: # ToDo: move to config?
        "--broad --broad-cutoff 0.1 -f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
        # ToDo: Update wrapper to check for " --broad$" or " --broad " instead of only "--broad" (line 47),
        #  then "--broad" in params can be removed here in params
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule macs2_callpeak_narrow:
    input:
        # Need to get treatment/input pairs as input
        treatment="/output/results/filtered/{basename}.sorted.bam",
        control="/output/results/filtered/{control}.sorted.bam"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("/output/results/macs2_callpeak/{basename}-{control}.narrow",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                #  "_treat_pileup.bdg",
                #  "_control_lambda.bdg",
                 # these output extensions internally set the --narrow option:
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    log:
        "/output/logs/macs2/callpeak.{basename}-{control}.narrow.log"
    params: # ToDo: move to config?
        "-f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule plot_macs_qc:
    input:
        get_plot_macs_qc_input()
    output:  #ToDo: add description to report caption
        summmary="/output/results/macs2_callpeak/plots/plot_{peak}_peaks_macs2_summary.txt",
        plot=report("/output/results/macs2_callpeak/plots/plot_{peak}_peaks_macs2.pdf", caption="/work/workflow/report/plot_macs2_qc.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "/output/logs/macs2_callpeak/plot_{peak}_peaks_macs2.log"
    shell:
        "Rscript /work/workflow/scripts/plot_macs_qc.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

