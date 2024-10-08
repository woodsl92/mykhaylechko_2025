Scripts used to perform the analysis described in "Comparative ASCL1 interactome analysis reveals CDK2-Cyclin A2 as suppressors of differentiation in MYCN-amplified neuroblastoma"
Mykhaylechko et al 2025

RNA-seq prepreprocessing - see guide at rnaseq_preprocessing/README.md

ASCL1 ChIP-seq prepreprocessing - see guide at ASCL1_chip_preprocessing/README.md
* Data aligned to hg38  
* Paired end data macs2 peak calling settings: -f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all
* Single-end data macs2 peak calling settings: -f BAM -g hs --SPMR --qvalue 0.05 --keep-dup all
* BigWig files generated with deeptools bamCoverage with CPM normalisation 

RNA-seq analysis - analysis_scripts/rnaseq.R

ASCL1 ChIP-seq analysis - analysis_scripts/script_ChIPseq.R, analysis_scripts/script_ChIP_peakannotations.R, analysis_scripts/script_chip_plots.txt

RNA-seq/ChIP-seq integration - analysis_scripts/script_combineddata_RNA_ChIP.R

Protein analysis - analysis_scripts/totalproteome.R, analysis_scripts/script_ASCL1_qPLEX-RIME.R



