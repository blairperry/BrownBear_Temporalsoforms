
library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(tximport)


## NOTE: Using IsoformSwitchAnalyzeR to characterize isoform switches identified by Tappas.

# Import data and create switchAnalyzeRlist -------------------------------
# Following: http://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#importing-data-into-r

dir <- '/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/isoforms/kallisto_new_Feb2022/kallisto_quants_Feb2022'

sample_info <- readxl::read_xlsx('/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/prevRNAseq_info/prevRNA_SampleInfo.xlsx') %>%
  filter(Sample != 'CFA_S9')

filenames <- str_split_fixed(sample_info$Sample,'[_]',2)[,1]

files <- file.path(dir, filenames, "abundance.h5")
names(files) <- str_split_fixed(filenames,'[_]',2)[,1]

all(file.exists(files))

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

id.convert <- read_csv('/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/isoforms/z_old_preFeb2022/new_merge.combined.renamed.mapping.txt')

# Format counts and abundance data
all.counts <- as.data.frame(txi.kallisto$counts) %>%
  rownames_to_column(var='old_id') %>%
  left_join(id.convert) %>%
  dplyr::select(-old_id) %>%
  dplyr::select(tx_id=new_id,everything()) %>%
  as.data.frame() %>%
  column_to_rownames(var='tx_id')

all.abundance <- as.data.frame(txi.kallisto$abundance) %>%
  rownames_to_column(var='old_id') %>%
  left_join(id.convert) %>%
  dplyr::select(-old_id) %>%
  dplyr::select(tx_id=new_id,everything()) %>%
  as.data.frame() %>%
  column_to_rownames(var='tx_id')


# Format adipose, liver, and muscle-specific exp design tables ------------

adi.exp.design <- read_tsv('analyses/tappas/adipose_tappas_design.tsv') %>%
  dplyr::select(sampleID=1,condition=2) %>%
  mutate(condition=case_when(
    condition==1 ~ '1_Active',
    condition==2 ~ '2_Hyperphagia',
    condition==3 ~ '3_Hibernation'
    ))

liv.exp.design <- read_tsv('analyses/tappas/liver_tappas_design.tsv') %>%
  dplyr::select(sampleID=1,condition=2) %>%
  mutate(condition=case_when(
    condition==1 ~ '1_Active',
    condition==2 ~ '2_Hyperphagia',
    condition==3 ~ '3_Hibernation'
  ))

mus.exp.design <- read_tsv('analyses/tappas/muscle_tappas_design.tsv') %>%
  dplyr::select(sampleID=1,condition=2) %>%
  mutate(condition=case_when(
    condition==1 ~ '1_Active',
    condition==2 ~ '2_Hyperphagia',
    condition==3 ~ '3_Hibernation'
  ))


# make switchList objects -------------------------------------------------

adi.SwitchList <- importRdata(
  isoformCountMatrix   = all.counts,
  isoformRepExpression = all.abundance,
  designMatrix         = adi.exp.design,
  isoformExonAnnoation = "~/Downloads/new_merge.combined.renamed_corrected.cds.gtf",
  isoformNtFasta       = '~/Downloads/new_merge.pbIDs.fa' ,
  showProgress = T,
  fixStringTieAnnotationProblem = F
)


liv.SwitchList <- importRdata(
  isoformCountMatrix   = all.counts,
  isoformRepExpression = all.abundance,
  designMatrix         = liv.exp.design,
  isoformExonAnnoation = "~/Downloads/new_merge.combined.renamed_corrected.cds.gtf",
  isoformNtFasta       = '~/Downloads/new_merge.pbIDs.fa' ,
  showProgress = T,
  fixStringTieAnnotationProblem = F
)


mus.SwitchList <- importRdata(
  isoformCountMatrix   = all.counts,
  isoformRepExpression = all.abundance,
  designMatrix         = mus.exp.design,
  isoformExonAnnoation = "~/Downloads/new_merge.combined.renamed_corrected.cds.gtf",
  isoformNtFasta       = '~/Downloads/new_merge.pbIDs.fa' ,
  showProgress = T,
  fixStringTieAnnotationProblem = F
)


# Read in and add DIU analysis q-values from TAPPAS -----------------------

adip.tappas.res <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene=number_gene,q_value)

adip.tappas.txInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene,tx_id=number_transcript) %>%
  left_join(adip.tappas.res) %>%
  dplyr::select(2,3)

liv.tappas.res <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene=number_gene,q_value)

liv.tappas.txInfo <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene,tx_id=number_transcript) %>%
  left_join(liv.tappas.res) %>%
  dplyr::select(2,3)

mus.tappas.res <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene=number_gene,q_value)

mus.tappas.txInfo <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  dplyr::select(gene,tx_id=number_transcript) %>%
  left_join(mus.tappas.res) %>%
  dplyr::select(2,3)



# Set gene_switch_q_value using q-value from Tappas
adi.isoformFeatures.temp <- adi.SwitchList$isoformFeatures %>%
  left_join(adip.tappas.txInfo,by=c('isoform_id'='tx_id')) %>%
  mutate(gene_switch_q_value = q_value) %>%
  dplyr::select(-q_value)

adi.SwitchList$isoformFeatures <- adi.isoformFeatures.temp

liv.isoformFeatures.temp <- liv.SwitchList$isoformFeatures %>%
  left_join(liv.tappas.txInfo,by=c('isoform_id'='tx_id')) %>%
  mutate(gene_switch_q_value = q_value) %>%
  dplyr::select(-q_value)

liv.SwitchList$isoformFeatures <- liv.isoformFeatures.temp

mus.isoformFeatures.temp <- mus.SwitchList$isoformFeatures %>%
  left_join(mus.tappas.txInfo,by=c('isoform_id'='tx_id')) %>%
  mutate(gene_switch_q_value = q_value) %>%
  dplyr::select(-q_value)

mus.SwitchList$isoformFeatures <- mus.isoformFeatures.temp


# Filter to only include adipose DIU genes ------------------------------------------------------------------

adi.diu <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  filter(diu_result=='DIU' & major_isoform_switching=='YES')

liv.diu <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  filter(diu_result=='DIU' & major_isoform_switching=='YES')

mus.diu <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>%
  janitor::clean_names() %>%
  filter(diu_result=='DIU' & major_isoform_switching=='YES')


adi.SwitchListFiltered <- subsetSwitchAnalyzeRlist(
  adi.SwitchList,
  adi.SwitchList$isoformFeatures$gene_id %in% adi.diu$number_gene
)

liv.SwitchListFiltered <- subsetSwitchAnalyzeRlist(
  liv.SwitchList,
  liv.SwitchList$isoformFeatures$gene_id %in% liv.diu$number_gene
)

mus.SwitchListFiltered <- subsetSwitchAnalyzeRlist(
  mus.SwitchList,
  mus.SwitchList$isoformFeatures$gene_id %in% mus.diu$number_gene
)

temp <- mus.SwitchListFiltered$isoformFeatures


# # Extracting nucleotide and amino acid sequences -----------------------------------------

adi.SwitchListAnalyzed <- extractSequence(
  adi.SwitchListFiltered,
  pathToOutput = 'analyses/isoformSwitchAnalyzeR/isoform_seqs/adipose/',
  writeToFile=T,
  # writeToFile=F, # set to false after the initial fasta file export
  alsoSplitFastaFile = F,
  onlySwitchingGenes = F, # set to false since we already filtered to switching genes from Tappas analyses
  removeLongAAseq = F
)

liv.SwitchListAnalyzed <- extractSequence(
  liv.SwitchListFiltered,
  pathToOutput = 'analyses/isoformSwitchAnalyzeR/isoform_seqs/liver/',
  # writeToFile=T,
  writeToFile=F, # set to false after the initial fasta file export
  onlySwitchingGenes = F, # set to false since we already filtered to switching genes from Tappas analyses
  alsoSplitFastaFile = F,
  removeLongAAseq = F
)


mus.SwitchListAnalyzed <- extractSequence(
  mus.SwitchListFiltered,
  pathToOutput = 'analyses/isoformSwitchAnalyzeR/isoform_seqs/muscle/',
  # writeToFile=T,
  writeToFile=F, # set to false after the initial fasta file export
  alsoSplitFastaFile = F,
  onlySwitchingGenes = F, # set to false since we already filtered to switching genes from Tappas analyses
  removeLongAAseq = F
)

###### Run pfam in terminal as follows:
# pfam_scan.pl -fasta ../isoform_seqs/adipose/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./adipose/adipose_pfam_scan_res.txt
# pfam_scan.pl -fasta ../isoform_seqs/liver/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./liver/liver_pfam_scan_res.txt
# pfam_scan.pl -fasta ../isoform_seqs/muscle/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./muscle/muscle_pfam_scan_res.txt


# Import pfam analysis results --------------------------------------------

adi.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = adi.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analyses/isoformSwitchAnalyzeR/pfam_scans/adipose/adipose_pfam_scan_res.txt',
  showProgress=T
)

liv.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = liv.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analyses/isoformSwitchAnalyzeR/pfam_scans/liver/liver_pfam_scan_res.txt',
  showProgress=T
)

mus.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = mus.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analyses/isoformSwitchAnalyzeR/pfam_scans/muscle/muscle_pfam_scan_res.txt',
  showProgress=T
)


Analyze for alternative splicing ----------------------------------------

adi.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = adi.SwitchListAnalyzed,showProgress = T
)

liv.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = liv.SwitchListAnalyzed,showProgress = T
)

mus.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = mus.SwitchListAnalyzed,showProgress = T
)


# Analyze for switch consequences
adi.SwitchListAnalyzed <- analyzeSwitchConsequences(
  switchAnalyzeRlist = adi.SwitchListAnalyzed,showProgress = T,consequencesToAnalyze = c('domains_identified','isoform_length','exon_number'
                                                                                      ,'intron_structure','ORF_seq_similarity','ORF_genomic',
                                                                                      'ORF_length','5_utr_seq_similarity','5_utr_length',
                                                                                      '3_utr_seq_similarity','3_utr_length',
                                                                                      'genomic_domain_position','domain_length')
    )

liv.SwitchListAnalyzed <- analyzeSwitchConsequences(
  switchAnalyzeRlist = liv.SwitchListAnalyzed,showProgress = T,consequencesToAnalyze = c('domains_identified','isoform_length','exon_number'
                                                                                         ,'intron_structure','ORF_seq_similarity','ORF_genomic',
                                                                                         'ORF_length','5_utr_seq_similarity','5_utr_length',
                                                                                         '3_utr_seq_similarity','3_utr_length',
                                                                                         'genomic_domain_position','domain_length')
)

mus.SwitchListAnalyzed <- analyzeSwitchConsequences(
  switchAnalyzeRlist = mus.SwitchListAnalyzed,showProgress = T,consequencesToAnalyze = c('domains_identified','isoform_length','exon_number'
                                                                                         ,'intron_structure','ORF_seq_similarity','ORF_genomic',
                                                                                         'ORF_length','5_utr_seq_similarity','5_utr_length',
                                                                                         '3_utr_seq_similarity','3_utr_length',
                                                                                         'genomic_domain_position','domain_length')
)

# Analyze ORFs
adi.SwitchListAnalyzed <- analyzeORF(
  adi.SwitchListAnalyzed
)

liv.SwitchListAnalyzed <- analyzeORF(
  liv.SwitchListAnalyzed
)

mus.SwitchListAnalyzed <- analyzeORF(
  mus.SwitchListAnalyzed
)

#
# # Plot summaries
# extractSplicingSummary(
#   aSwitchListAnalyzed,
#   asFractionTotal = FALSE,
#   plotGenes=FALSE
# )
#
# extractConsequenceSummary(
#   aSwitchListAnalyzed,
#   consequencesToAnalyze='all',
#   plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
#   asFractionTotal = FALSE      # enables analysis of fraction of significant features
# )


# Save as Rdata for quicker import later on -------------------------------

#save(adi.SwitchListAnalyzed,liv.SwitchListAnalyzed,mus.SwitchListAnalyzed,file='analyses/isoformSwitchAnalyzeR/SwitchLists_03.14.22.Rdata')

# load('analyses/isoformSwitchAnalyzeR/SwitchLists_03.14.22.Rdata')

# Set gene name field == gene id (to work with automatic titles in isoform plots)
adi.SwitchListAnalyzed$isoformFeatures$gene_name <- adi.SwitchListAnalyzed$isoformFeatures$gene_id
liv.SwitchListAnalyzed$isoformFeatures$gene_name <- liv.SwitchListAnalyzed$isoformFeatures$gene_id
mus.SwitchListAnalyzed$isoformFeatures$gene_name <- mus.SwitchListAnalyzed$isoformFeatures$gene_id

# Plot specific isoform details

# GRB14
switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = c('PB.3206.1','PB.3206.5','PB.3206.10'))

switchPlotTranscript(adi.SwitchListAnalyzed, gene = 'TCIRG1')
switchPlotTranscript(adi.SwitchListAnalyzed, gene = 'AKT2')

# SHC1 (liver)
switchPlotTranscript(liv.SwitchListAnalyzed, gene = 'SHC1')

# ATK2 Adipose and muscle comparison
switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = c('PB.24366.2','PB.24366.5','PB.24366.14','PB.24366.20','PB.24366.23'))

switchPlotTranscript(mus.SwitchListAnalyzed, isoform_id = c('PB.24366.25','PB.24366.11','PB.24366.3','PB.24366.4'))

# Other isoforms from candidate set:

# Appl1 in liver and muscle
switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = c('PB.8915.4','PB.8915.2','PB.8915.16','PB.8915.10'))
switchPlotTranscript(mus.SwitchListAnalyzed, isoform_id = c('PB.8915.4','PB.8915.2','PB.8915.16','PB.8915.10'))

# SHC1 in liver
switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = c('PB.25956.7','PB.25956.6'))

# SORBS1 in liver
switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = c('PB.4750.3','PB.4750.24'))

#TCIRG1 in adipose
switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = c('PB.25002.4','PB.25002.1'))


# Supplementary Figure Building -------------------------------------------

#### NOTE: code below requires that the DEAvsDIU_03.01.22.R script (or later date if updated after this message) has been run 
#### through the definition of the plot_isoforms_and_gene and all associated data from that script is loaded in the R environment

## AKT2 adipose, muscle
p.akt.adi <- plot_isoforms_and_gene('AKT2',use_color=T,tissue='adipose',use_labels=T) + scale_y_log10()
p.akt.mus <- plot_isoforms_and_gene('AKT2',use_color=T,tissue='muscle',use_labels=T) + scale_y_log10()

p.akt.iso.adi <- switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = unique(p.akt.adi$data$number_transcript)) + ggtitle('')
p.akt.iso.mus <- switchPlotTranscript(mus.SwitchListAnalyzed, isoform_id = unique(p.akt.mus$data$number_transcript)) + ggtitle('')

p.akt.adi + p.akt.iso.adi + p.akt.mus + p.akt.iso.mus + plot_layout(ncol = 2,guides='collect',widths = c(1,2))


## APPL1 liver, muscle
p.appl1.liv <- plot_isoforms_and_gene('APPL1',use_color=T,tissue='liver',use_labels=T) + scale_y_log10()
p.appl1.mus <- plot_isoforms_and_gene('APPL1',use_color=T,tissue='muscle',use_labels=T) + scale_y_log10()

p.appl1.iso.liv <- switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = unique(p.appl1.liv$data$number_transcript)) + ggtitle('')
p.appl1.iso.mus <- switchPlotTranscript(mus.SwitchListAnalyzed, isoform_id = unique(p.appl1.mus$data$number_transcript)) + ggtitle('')

p.appl1.liv + p.appl1.iso.liv + p.appl1.mus + p.appl1.iso.mus + plot_layout(ncol = 2,guides='collect',widths = c(1,2))

## GRB14 adipose
p.grb14.adi <- plot_isoforms_and_gene('GRB14',use_color=T,tissue='adipose',use_labels=T) + scale_y_log10()
p.grb14.iso.adi <- switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = unique(p.grb14.adi$data$number_transcript)) + ggtitle('')
p.grb14.adi + p.grb14.iso.adi + plot_layout(ncol = 2,guides='collect',widths = c(1,2))

## SHC1 liver
p.shc1.liv <- plot_isoforms_and_gene('SHC1',use_color=T,tissue='liver',use_labels=T) + scale_y_log10()
p.shc1.iso.liv <- switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = unique(p.shc1.liv$data$number_transcript)) + ggtitle('')
p.shc1.liv + p.shc1.iso.liv + plot_layout(ncol = 2,guides='collect',widths = c(1,2))

## SORBS1 liver
p.sorbs1.liv <- plot_isoforms_and_gene('SORBS1',use_color=T,tissue='liver',use_labels=T) + scale_y_log10()
p.sorbs1.iso.liv <- switchPlotTranscript(liv.SwitchListAnalyzed, isoform_id = unique(p.sorbs1.liv$data$number_transcript)) + ggtitle('')
p.sorbs1.liv + p.sorbs1.iso.liv + plot_layout(ncol = 2,guides='collect',widths = c(1,2))

## TCIRG1 adipose
p.tcirg1.adi <- plot_isoforms_and_gene('TCIRG1',use_color=T,tissue='adipose',use_labels=T) + scale_y_log10()
p.tcirg1.iso.adi <- switchPlotTranscript(adi.SwitchListAnalyzed, isoform_id = unique(p.tcirg1.adi$data$number_transcript)) + ggtitle('')
p.tcirg1.adi + p.tcirg1.iso.adi + plot_layout(ncol = 2,guides='collect',widths = c(1,2))



