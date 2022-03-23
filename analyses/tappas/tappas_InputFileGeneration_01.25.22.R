

library(tidyverse)
library(tximport)

# Read in Kallisto results and combine ------------------------------------

dir <- '/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/isoforms/kallisto_new_Feb2022/kallisto_quants_Feb2022'

sample_info <- readxl::read_xlsx('/Volumes/WorkingDrive_BWP/_bearNecropsyRNAseq_Sept2021/data/prevRNAseq_info/prevRNA_SampleInfo.xlsx') %>% 
  filter(Sample != 'CFA_S9') # Note: sample was excluded - previously shown to contain hair follicle 

filenames <- str_split_fixed(sample_info$Sample,'[_]',2)[,1]

files <- file.path(dir, filenames, "abundance.h5")
names(files) <- str_split_fixed(filenames,'[_]',2)[,1]

all(file.exists(files))

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

head(txi.kallisto$counts)

# Read in ID conversion table (BBEAR to PB IDs) ---------------------------

id.convert <- read_csv('data/new_merge.combined.renamed.mapping.txt') # Read in table to convert between old and new tx ids

allTissue.counts <- as.data.frame(txi.kallisto$counts) %>% # Reading in the count estimates from Kallisto
  rownames_to_column(var='old_id') %>% 
  left_join(id.convert) %>% 
  select(-old_id) %>% 
  select(tx_id=new_id,everything()) %>% 
  as.data.frame() %>% 
  column_to_rownames(var='tx_id') %>% 
  rownames_to_column()



#############################################################################
######################    Adipose    ########################################
#############################################################################

# Make experiment design table
adip.exp <- sample_info %>% 
  select(sampleid = 1, timepoint = 3, group = 2) %>% 
  mutate(sampleid = str_split_fixed(sampleid,'[_]',2)[,1],
         timepoint = ifelse(timepoint == 'Active',1,ifelse(timepoint=='Hyperphagia',2,3)),
         group = ifelse(group == 'Adipose',1,ifelse(group =='Liver',2,3))) %>% 
  filter(group==1) %>% 
  dplyr::select(Sample=1,Time=2,Group=3) %>% 
  arrange(Time)
adip.exp


adip.counts <- as.data.frame(allTissue.counts[,which(names(allTissue.counts) %in% union('rowname',adip.exp$Sample))]) %>% 
  column_to_rownames('rowname')


write.table(adip.counts,'analyses/tappas/adipose_tappas_expMatrix_02.04.22.tsv',sep = '\t',quote = F,row.names = T)
write_tsv(adip.exp,'analyses/tappas/adipose_tappas_design.tsv')


#############################################################################
######################    Liver    ########################################
#############################################################################

# Make experiment design table
liv.exp <- sample_info %>% 
  select(sampleid = 1, timepoint = 3, group = 2) %>% 
  mutate(sampleid = str_split_fixed(sampleid,'[_]',2)[,1],
         timepoint = ifelse(timepoint == 'Active',1,ifelse(timepoint=='Hyperphagia',2,3)),
         group = ifelse(group == 'Adipose',1,ifelse(group =='Liver',2,3))) %>% 
  filter(group==2) %>% 
  dplyr::select(Sample=1,Time=2,Group=3) %>% 
  arrange(Time)
liv.exp


liv.counts <- as.data.frame(allTissue.counts[,which(names(allTissue.counts) %in% union('rowname',liv.exp$Sample))]) %>% 
  column_to_rownames('rowname')



write.table(liv.counts,'analyses/tappas/liver_tappas_expMatrix_02.04.22.tsv',sep = '\t',quote = F,row.names = T)
write_tsv(liv.exp,'analyses/tappas/liver_tappas_design.tsv')

#############################################################################
######################    Muscle    ########################################
#############################################################################

# Make experiment design table
mus.exp <- sample_info %>% 
  select(sampleid = 1, timepoint = 3, group = 2) %>% 
  mutate(sampleid = str_split_fixed(sampleid,'[_]',2)[,1],
         timepoint = ifelse(timepoint == 'Active',1,ifelse(timepoint=='Hyperphagia',2,3)),
         group = ifelse(group == 'Adipose',1,ifelse(group =='Liver',2,3))) %>% 
  filter(group==3) %>% 
  dplyr::select(Sample=1,Time=2,Group=3) %>% 
  arrange(Time)
mus.exp

mus.counts <- as.data.frame(allTissue.counts[,which(names(allTissue.counts) %in% union('rowname',mus.exp$Sample))]) %>% 
  column_to_rownames('rowname')

write.table(mus.counts,'analyses/tappas/muscle_tappas_expMatrix_02.04.22.tsv',sep = '\t',quote = F,row.names = T)
write_tsv(mus.exp,'analyses/tappas/muscle_tappas_design.tsv')

