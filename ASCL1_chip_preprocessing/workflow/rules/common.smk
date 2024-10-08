from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####

configfile: "/work/config/config.yaml"
validate(config, schema="/work/workflow/schemas/config.schema.yaml")

basenames = pd.read_csv(config["samples"], sep="\t", dtype = str).set_index("basename", drop=False)
basenames.index.names = ["basename_id"]
validate(basenames, schema="/work/workflow/schemas/samples.schema.yaml")
print(basenames)

conditions = pd.read_csv(config["samples"], sep="\t", dtype = str).set_index("condition", drop=False)
conditions.index.names = ["condition_id"]
validate(conditions, schema="/work/workflow/schemas/samples.schema.yaml")
print(conditions.index)
groups = {k: list(v) for k,v in conditions.groupby("condition")["basename"]}
print(groups)

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["basename", "unit"], drop=False)
units.index.names = ["basename_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="/work/workflow/schemas/units.schema.yaml")

report: "/work/workflow/report/workflow.rst"

##### wildcard constraints #####

wildcard_constraints:
    unit = "|".join(units["unit"]),
    basename = "|".join(basenames.index),
    condition = "|".join(conditions.index)


####### helpers ###########

def is_single_end(basename, unit):
    """Determine whether unit is single-end."""
    fq2_present = pd.isnull(units.loc[(basename, unit), "fq2"])
    if isinstance(fq2_present, pd.core.series.Series):
        # if this is the case, get_fastqs cannot work properly
        raise ValueError(
            f"Multiple fq2 entries found for basename-unit combination {basename}-{unit}.\n"
            "This is most likely due to a faulty units.tsv file, e.g. "
            "a unit name is used twice for the same basename.\n"
            "Try checking your units.tsv for duplicates."
        )
    return fq2_present

def get_individual_fastq(wildcards):
    """Get individual raw FASTQ files from unit sheet, based on a read (end) wildcard"""
    if ( wildcards.read == "0" or wildcards.read == "1" ):
        return units.loc[ (wildcards.basename, wildcards.unit), "fq1" ]
    elif wildcards.read == "2":
        return units.loc[ (wildcards.basename, wildcards.unit), "fq2" ]

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.basename, wildcards.unit):
        return units.loc[ (wildcards.basename, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.basename, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def is_control(basename):
    control = basenames.loc[basename]["control"]
    return pd.isna(control) or pd.isnull(control)

def is_control_condition(groups):
	# Get condition info
	control_condition = conditions.loc[groups]["condition"]
	# Does condition contain string "control", Return bool
	return "input" in control_condition

def get_rose_peak_input(groups):
    if not is_control_condition(groups):
        return "/output/results/macs2_merged/{group}.narrow_peaks.narrowPeak.bed"
def get_rose_bam_input(groups):
    if not is_control_condition(groups):
        return "/output/results/merged_bam/{group}.subsample.sorted.bam"

def get_sample_control_peak_combinations_list():
    sam_contr = []
    for basename in basenames.index:
        if not is_control(basename):
            sam_contr.extend(expand(["{basename}-{control}.{peak}"], basename = basename, control = basenames.loc[basename]["control"], peak = config["params"]["peak-analysis"]))
    return sam_contr

def get_peaks_count_plot_input():
    return expand(
        "/output/results/macs2_callpeak/peaks_count/{sam_contr_peak}.{peak}.peaks_count.tsv",
        sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
    )

def get_frip_score_input():
    return expand(
        "/output/results/intersect/{sam_contr_peak}.{peak}.peaks_frip.tsv",
        sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
    )

def get_plot_macs_qc_input():
    return expand(
        "/output/results/macs2_callpeak/{sam_contr_peak}_peaks.{peak}Peak",
        sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
    )

def get_plot_homer_annotatepeaks_input():
    return expand("/output/results/homer/annotate_peaks/{sam_contr_peak}.{peak}_peaks.annotatePeaks.txt",
        sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
    )

def get_map_reads_input(wildcards):
    if is_single_end(wildcards.basename, wildcards.unit):
        return "/output/results/trimmed/{basename}_{unit}_trimmed.fq.gz"
    return ["/output/results/trimmed/{basename}_{unit}_R1_val_1.fq.gz", "/output/results/trimmed/{basename}_{unit}_R2_val_2.fq.gz"]

def get_read_group(wildcards):
    """Denote basename name and platform in read group."""
    return r"-R '@RG\tID:{basename}-{unit}\tSM:{basename}-{unit}\tPL:{platform}'".format(
        basename=wildcards.basename,
        unit=wildcards.unit,
        platform=units.loc[(wildcards.basename, wildcards.unit), "platform"])

def get_multiqc_input(wildcards):
    multiqc_input = []
    for (basename, unit) in units.index:
        reads = [ "1", "2" ]
        if is_single_end(basename, unit):
            reads = [ "0" ]
            multiqc_input.extend(expand (
                [
                    "/output/logs/trim_galore/{basename}_{unit}.log",
                    "/output/results/trimmed/{basename}_{unit}_trimmed.fq.gz",
                    "/output/results/trimmed/{basename}_{unit}.fq.gz_trimming_report.txt"
                    ],
            basename = basename, unit = unit))
        else:
            multiqc_input.extend(expand (
                [
                    "/output/logs/trim_galore/{basename}_{unit}_pe.log",
                    "/output/results/fastqc/{basename}_{unit}_R1_fastqc.zip",
                    "/output/results/fastqc/{basename}_{unit}_R2_fastqc.zip",
                    "/output/results/fastqc/{basename}_{unit}_R1_fastqc.html",
                    "/output/results/fastqc/{basename}_{unit}_R2_fastqc.html"
                    ],
            basename = basename, unit = unit))

    for basename in basenames.index:
        multiqc_input.extend(
            expand (
                [
                    "/output/results/picard_dedup/{basename}.metrics.txt",
                    "/output/results/filtered/{basename}.flagstat",
                    "/output/results/filtered/{basename}.filtered.idxstats",
                    "/output/results/filtered/{basename}.mapped.stats.txt",
                    "/output/results/deeptools/plot_profile_data.tab",
                    "/output/results/phantompeakqualtools/{basename}.phantompeak.spp.out",
                    "/output/results/phantompeakqualtools/{basename}.spp_correlation_mqc.tsv",
                    "/output/results/phantompeakqualtools/{basename}.spp_nsc_mqc.tsv",
                    "/output/results/phantompeakqualtools/{basename}.spp_rsc_mqc.tsv"
                ],
                basename = basename
            )
        )
        if not is_control(basename):
            multiqc_input.extend(
                expand (
                    [
                        "/output/results/deeptools/{basename}-{control}.fingerprint_qcmetrics.txt",
                        "/output/results/deeptools/{basename}-{control}.fingerprint_counts.txt",
                        "/output/results/macs2_callpeak/{basename}-{control}.{peak}_peaks.xls"
                    ],
                basename = basename,
                control = basenames.loc[basename]["control"],
                peak = config["params"]["peak-analysis"]
            )
        )
        if config["params"]["lc_extrap"]:
                multiqc_input.extend( expand(["/output/results/preseq/{basename}.lc_extrap"], basename = basename))
    return multiqc_input

def all_input(wildcards):

    wanted_input = []

    # QC with fastQC and multiQC
    wanted_input.extend([
        "/output/results/multiqc/multiqc.html"
    ])

    # trimming reads
    for (basename, unit) in units.index:
        if is_single_end(basename, unit):
            wanted_input.extend(expand(
                    [
                        "/output/results/trimmed/{basename}_{unit}_trimmed.fq.gz",
                        "/output/results/trimmed/{basename}_{unit}.fq.gz_trimming_report.txt"
                    ],
                    basename = basename,
                    unit = unit
                )
            )
        else:
            wanted_input.extend(
                expand (
                    [
                        "/output/results/trimmed/{basename}_{unit}_R1_val_1.fq.gz",
                        "/output/results/trimmed/{basename}_{unit}_R2_val_2.fq.gz",
                        "/output/results/trimmed/{basename}_{unit}_R1.fq.gz_trimming_report.txt",
                        "/output/results/trimmed/{basename}_{unit}_R2.fq.gz_trimming_report.txt"
                    ],
                    basename = basename,
                    unit = unit
            )
        )

    # mapping, merging and filtering bam-files
    for basename in basenames.index:
        wanted_input.extend(
            expand (
                [
                    "/output/results/deeptools/plot_profile.pdf",
                    "/output/results/deeptools/heatmap.pdf",
                    "/output/results/deeptools/heatmap_matrix.tab",
                    "/output/results/phantompeakqualtools/{basename}.phantompeak.pdf",
                    "/output/results/filtered/{basename}.sorted.bam"
                ],
                basename = basename
            )
        )
        if not is_control(basename):
            wanted_input.extend(
                expand(
                    [
                        "/output/results/deeptools/{basename}-{control}.plot_fingerprint.pdf"
                    ],
                    basename = basename,
                    control = basenames.loc[basename]["control"],
                    peak = config["params"]["peak-analysis"]
                )
            )
            if config["params"]["peak-analysis"] == "broad":
                wanted_input.extend(
                    expand(
                        [
                            "/output/results/macs2_callpeak/{basename}-{control}.{peak}_peaks.gappedPeak"
                        ],
                        basename = basename,
                        control = basenames.loc[basename]["control"],
                        peak = config["params"]["peak-analysis"]
                    )
                )
            if config["params"]["peak-analysis"] == "narrow":
                wanted_input.extend(
                    expand(
                        [
                            "/output/results/macs2_callpeak/{basename}-{control}.{peak}_summits.bed"

                        ],
                        basename = basename,
                        control = basenames.loc[basename]["control"],
                        peak = config["params"]["peak-analysis"]
                    )
                )


    return wanted_input
