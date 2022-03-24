# Analysis of temporal isoform usage in brown bears (*Ursus arctos*)

This repo contains code for all data processing and analysis used in our study of timeseries isoform usage associated with hibernation in brown bears.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
- Brown bear mRNA-seq data (from [Jansen et al. 2019](https://www.nature.com/articles/s42003-019-0574-4))
  - Raw data available at NCBI BioProject: [PRJNA413091](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA413091)
    - Three physiological stages: active, hyperphagia, hibernation
    - Three metabolic tissues: liver, muscle, adipose
- Iso-seq Transcriptome Assembly (from [Tseng et al. 2021](https://academic.oup.com/g3journal/article/12/3/jkab422/6472356))
  - Raw data available at NCBI BioProject: [PRJNA727613](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA727613)
  - Assembly and annotation available at: https://github.com/jokelley/brownbear-isoseq-act-hib



## Software
The following software was used processing and analysis:
- Kallisto v0.46.1 (isoform-level expression quantification)
- tappAS v1.0.7 (differential isoform usage and gene expression analyses)
- pfam_scan.pl v1.6 (identification of protein domains in isoform sequences)
- R v4.1.0 (various analyses and plotting)
- R packages:
  - ggplot2 v3.3.5 (plotting)
  - tidyverse v1.3.1 (data wrangling)
  - tximport v1.22.0 (importing Kallisto output)
  - Rrvgo v1.6.0 (reducing GO term redundancy)
  - ggupset v0.3.0 (upset plots)
  - ggvenn v0.1.9 (venn diagrams)
  - pheatmap v1.0.12 (heatmaps)
  - IsoformSwitchAnalyzeR v1.16.0 (isoform analysis and plotting)

## Contents
1. [Expression quantification](#1-isoform-level-expression-quantification-with-kallisto)
2. [Formatting tappAS input files](#2-formatting-kallisto-results-for-input-into-tappas)
3. [Differential isoform usage and gene expression analyses in tappAS](#3-analyses-of-differential-isoform-usage-and-gene-expression-in-tappas)
4. [Plotting differential isoform usage results](4-plotting-of-differential-isoform-usage-diu-results)
5. [Comparing differential isoform usage and gene-level expression](#5-comparisons-of-genes-with-differential-expression-and-differential-isoform-usage)
6. [GO term analysis](#6-gene-ontology-go-overrepresentation-analyses)
7. [Candidate gene analysis: insulin signaling genes](#7-analysis-of-insulin-signaling-candidate-genes)
8. [Characterizing protein domains in isoform switches](https://github.com/blairperry/BrownBear_Temporalsoforms#8-additional-characterization-of-isoform-switches)

---
### 1. Isoform-level expression quantification with Kallisto
These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.

#### Index transcriptome:
```bash
kallisto index -i transcripts.idx {transcriptome fasta}
```

#### Quantify isoform-level expression:
For each sample:
```bash
kallisto quant -i transcripts.idx --rf-stranded -o ./kallisto_quants_Feb2022/{sample name} -b 100 -t 5 {read1 file} {read2 file}
```
---
### 2. Formatting Kallisto results for input into tappAS
Kallisto result files were imported and combined in R using the tximport package. Raw count and experimental design files were then constructed for each tissue and exported for use with tappAS.

- Link to Rscript: [analyses/tappas/tappas_InputFileGeneration_01.25.22.R](analyses/tappas/tappas_InputFileGeneration_01.25.22.R)
---

### 3. Analyses of differential isoform usage and gene expression in tappAS
Raw count and experimental design files generated in the previous step were loaded into the [tappAS GUI](https://app.tappas.org/) to create a separate project for each tissue. Relevant import and analysis parameters are provided below.

- Import parameters:
  - Experiment type: Single series time-course
  - Normalize data by TMM method: Yes (default)
  - Expression values filter:
    - Low count values cutoff: 1.0 CPM (default)
    - Coefficient of variation cutoff: 100% (default)
  - Transcripts filter: Do not filter transcripts (default)
- Differential Isoform Usage (DIU) Analysis:
  - Data type: Transcripts (default)
  - Analysis method: maSigPro (default)
  - Significance level: 0.05 (default)
  - Polynomial degree: 2 (default)
  - Filtering
    - Filter minor isoforms: Fold filtering (default)
    - Fold expression difference: 2.0 (default)
- Differential Gene-level Expression Analysis (DEA):
  - Data type: genes (default)
  - Analysis method: maSigPro (default)
  - Polynomial degree: 2 (default)
  - Significance level: 0.05 (default)
  - Fold Change: 2 (default)
  - R^2 Cutoff: 0.7 (default)
  - Max K Clusters: 9 (default)
    - Use mclust for optimal number - Yes

---

### 4. Plotting of differential isoform usage (DIU) results
Result tables from DIU analysis were exported from tappAS and plotted using R. The following script contains code for upset plots comparing timepoint-specific isoform switches, venn diagram comparison of DIU genes across tissues, a custom plotting function for plotting isoform expression over time, and exporting gene lists for subsequent GO analyses (i.e., all DIU genes in each tissue, genes with timepoint specific switches in each tissue)
- Link to Rscript: [analyses/tappas/tappas_DIUResultParsing_02.24.22.R](analyses/tappas/tappas_DIUResultParsing_02.24.22.R)

---
### 5. Comparisons of genes with differential expression and differential isoform usage
Result tables from differential gene-level expression analses were exported from tappAS and imported into R. The following script contains code for importing, intersecting, and plotting genes found to be significant in differential isoform usage and/or differential gene-level expression analyses, as well as exporting gene lists for GO analyses (intersects genes singificant in both analyses per tissue), and a modified version of the custom plotting script from above that includes gene-level expression.
- Link to Rscript: [analyses/tappas/DEAvsDIU_03.18.22.R](analyses/tappas/DEAvsDIU_03.18.22.R)

---
### 6. Gene ontology (GO) overrepresentation analyses
GO term overrepresentation analyses were conducted with the online [PantherDB Gene List Analysis tool](pantherdb.org) with default parameters. For each target gene list, a list of all multi-isoform genes passing tappAS input filter criteria (described above) was used as the background.

The following script was used to import and plot statistically significant GO analysis results.
- Link to Rscript: [analyses/GOAnalyses/GOResultParsing_03.14.22.R](analyses/GOAnalyses/GOResultParsing_03.14.22.R)

---
### 7. Analysis of insulin signaling candidate genes
Genes involved in the insulin signaling pathway were downloaded from the online [UniProt](https://www.uniprot.org/) database using the following search query:
```
goa:("insulin receptor signaling pathway [8286]") AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"
```

The following script contains code used to investigate and plot results of differential isoform usage and differential gene-level expression for insulin signaling candidate genes:
- Link to Rscript: [analyses/candidate_Genes/InsulinSignaling_CandGeneExploration_03.07.22.R](analyses/candidate_Genes/InsulinSignaling_CandGeneExploration_03.07.22.R)

---

### 8. Additional characterization of isoform switches
To better understand potential consequences of isoform switches for key genes (i.e., gains and losses of protein domains), we performed additional isoform characterization using the IsoformSwitchAnalyzeR package. The following script contains code used to import and filter DIU isoforms from tappAS analysis results (we did not use IsoformSwitchAnalyzeR for isoform switch analysis), export isoform sequences for analysis with PFAM and import PFAM results, characterize various isoform features (i.e., sensitivity to non-sense mediated decay; NMD), and plot isoform structure and features.

- Link to Rscript: [analyses/isoformSwitchAnalyzeR/isoformSwitchAnalyzeR_03.07.22.R](analyses/isoformSwitchAnalyzeR/isoformSwitchAnalyzeR_03.07.22.R)

- The following commands were used to identify protein domains using PFAM:
```
pfam_scan.pl -fasta ../isoform_seqs/adipose/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./adipose/adipose_pfam_scan_res.txt
pfam_scan.pl -fasta ../isoform_seqs/liver/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./liver/liver_pfam_scan_res.txt
pfam_scan.pl -fasta ../isoform_seqs/muscle/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./muscle/muscle_pfam_scan_res.txt
```
  - NOTE: these commands all used the Pfam-A.hmm database downloaded [here](http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/).
