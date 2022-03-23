
library(tidyverse)

### Want to find genes that are NOT differentially expressed (DEA) but have significant DIU+Isoform Switching
### Particularly interested in muscle and liver, which showed little DE in previous studies


# Read in gene info tables to get isoform number --------------------------

adi.geneInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  dplyr::select(number_gene,multiIso)

adi.geneInfo %>% filter(multiIso==T) %>% nrow()

liv.geneInfo <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  dplyr::select(number_gene,multiIso)

liv.geneInfo %>% filter(multiIso==T) %>% nrow()

mus.geneInfo <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  dplyr::select(number_gene,multiIso)

mus.geneInfo %>% filter(multiIso==T) %>% nrow()


# Read in DEA results ----------------------------------------------------

adi.dea <- read.table('analyses/tappas/adipose_res/Adipose_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(tissue = 'Adipose',result='Differential gene expression') %>% 
  left_join(adi.geneInfo,by=c('gene_id'='number_gene'))

adi.dea %>% filter(multiIso==T) %>% nrow()

liv.dea <- read.table('analyses/tappas/liver_res/Liver_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(tissue = 'Liver',result='Differential gene expression')%>% 
  left_join(liv.geneInfo,by=c('gene_id'='number_gene'))

liv.dea %>% filter(multiIso==T) %>% nrow()

mus.dea <- read.table('analyses/tappas/muscle_res/Muscle_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(tissue = 'Muscle',result='Differential gene expression')%>% 
  left_join(mus.geneInfo,by=c('gene_id'='number_gene'))

mus.dea %>% filter(multiIso==T) %>% nrow()


all.dea <- adi.dea %>% 
  bind_rows(liv.dea,mus.dea ) %>% 
  dplyr::select(gene = gene_id,tissue,result,multiIso)

# DE venn

vennSets.de <- list('Adipose'= adi.dea$gene_id,
                 'Liver'  = liv.dea$gene_id,
                 'Muscle' = mus.dea$gene_id)

ggvenn::ggvenn(vennSets.de,fill_color = c('#E69F00','#56B4E9','#009E73'))



# Read in DIU results ----------------------------------------------------

adi.diu <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  # filter(diu_result == 'DIU') %>% 
  mutate(dea = ifelse(number_gene %in% adi.dea$gene_id,T,F)) %>% 
  mutate(tissue = 'Adipose',result='Differential isoform usage')

liv.diu <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  # filter(diu_result == 'DIU') %>% 
  mutate(dea = ifelse(number_gene %in% liv.dea$gene_id,T,F)) %>% 
  mutate(tissue = 'Liver',result='Differential isoform usage')

mus.diu <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  # filter(diu_result == 'DIU') %>% 
  mutate(dea = ifelse(number_gene %in% mus.dea$gene_id,T,F)) %>% 
  mutate(tissue = 'Muscle',result='Differential isoform usage')

all.diu <- adi.diu %>% 
  bind_rows(liv.diu,mus.diu) %>% 
  dplyr::select(gene=number_gene,tissue,result)

# DIU venn

adi.res.diu <- adi.diu %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')
liv.res.diu <- liv.diu %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')
mus.res.diu <- mus.diu %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')

vennSets <- list('Adipose'= adi.res.diu$number_gene,
                 'Liver'  = liv.res.diu$number_gene,
                 'Muscle' = mus.res.diu$number_gene)

ggvenn::ggvenn(vennSets,fill_color = c('#E69F00','#56B4E9','#009E73'))


# Plot overlap between DEA and DIU sets in each tissue --------------------

all.combined <- all.dea %>%
  filter(multiIso==T) %>% 
  dplyr::select(-multiIso) %>% 
  bind_rows(all.diu) %>% 
  group_by(gene,tissue) %>% 
  add_count() %>% 
  mutate(result_simple = ifelse(n == 2, 'Both', result)) %>% 
  dplyr::select(gene,tissue,result_simple) %>% 
  unique() %>% 
  ungroup() %>% 
  group_by(tissue,result_simple) %>% 
  tally() %>% 
  mutate(result_simple = factor(result_simple,levels=c('Differential gene expression','Differential isoform usage','Both')))

p1 <- ggplot(all.combined,aes(x=tissue,y=n,fill=result_simple,label=n)) +
  geom_bar(stat='identity',position = position_dodge()) +
  geom_text(position = position_dodge(.9),vjust= -0.5,hjust=0.5) +
  labs(y='Number of Genes',x='',fill='') +
  scale_fill_manual(values=c('Differential gene expression'='darkblue','Differential isoform usage'='Seagreen','Both'='goldenrod')) +
  theme_classic(base_size = 16) + theme(legend.position = 'bottom')

p1

# ggsave('_ms/figures/fig_pieces/Fig2_DEAvsDIU_03.07.22.pdf')



# Relevant gene lists for GO characterization -----------------------------

## For now, doing intersection of DIU and DEA analyses 
## and background as all multi-isoform genes that went into DIU analysis (subset of those that went into DEA)

adi.intersect <- all.dea %>%
  filter(multiIso==T) %>% 
  dplyr::select(-multiIso) %>% 
  bind_rows(all.diu) %>% 
  group_by(gene,tissue) %>% 
  add_count() %>% 
  mutate(result_simple = ifelse(n == 2, 'Both', result)) %>% 
  dplyr::select(gene,tissue,result_simple) %>% 
  unique() %>% 
  ungroup() %>% 
  filter(tissue == 'Adipose',result_simple == 'Both') %>% 
  dplyr::select(gene)

liv.intersect <- all.dea %>%
  filter(multiIso==T) %>% 
  dplyr::select(-multiIso) %>% 
  bind_rows(all.diu) %>% 
  group_by(gene,tissue) %>% 
  add_count() %>% 
  mutate(result_simple = ifelse(n == 2, 'Both', result)) %>% 
  dplyr::select(gene,tissue,result_simple) %>% 
  unique() %>% 
  ungroup() %>% 
  filter(tissue == 'Liver',result_simple == 'Both')%>% 
  dplyr::select(gene)

mus.intersect <- all.dea %>%
  filter(multiIso==T) %>% 
  dplyr::select(-multiIso) %>% 
  bind_rows(all.diu) %>% 
  group_by(gene,tissue) %>% 
  add_count() %>% 
  mutate(result_simple = ifelse(n == 2, 'Both', result)) %>% 
  dplyr::select(gene,tissue,result_simple) %>% 
  unique() %>% 
  ungroup() %>% 
  filter(tissue == 'Muscle',result_simple == 'Both')%>% 
  dplyr::select(gene)

adi.multiIsoBackground <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(number_gene)

liv.multiIsoBackground <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(number_gene)

mus.multiIsoBackground <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(number_gene)

# write_tsv(adi.intersect,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Adipose_DeAndDiuGenes.txt',col_names = F)
# write_tsv(liv.intersect,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Liver_DeAndDiuGenes.txt',col_names = F)
# write_tsv(mus.intersect,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Muscle_DeAndDiuGenes.txt',col_names = F)

# write_tsv(adi.multiIsoBackground,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Adipose_Background_allInputGenes.txt',col_names = F)
# write_tsv(liv.multiIsoBackground,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Liver_Background_allInputGenes.txt',col_names = F)
# write_tsv(mus.multiIsoBackground,'analyses/tappas/_DIUvsDEA_GeneSetsForGO/Muscle_Background_allInputGenes.txt',col_names = F)


### Plot a few genes from intersect

insr.geneList <- read_tsv('analyses/candidate_Genes/InsulinSignaling_Uniprot_CandGenes_AllSymbols.txt',col_names = 'gene')

adi.intersect %>% filter(gene %in% insr.geneList$gene)

adi.res <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Adipose') %>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

adi.txExp <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

adi.geneExp <-  read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  mutate(overall_direction = ifelse(x1_t3_mean_exp_level < x1_t1_mean_exp_level,'Down','Up'))

liv.intersect %>% filter(gene %in% insr.geneList$gene)

liv.res <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Liver') %>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

liv.txExp <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

liv.geneExp <-  read_tsv('analyses/tappas/liver_res/Liver_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  mutate(overall_direction = ifelse(x2_t3_mean_exp_level < x2_t1_mean_exp_level,'Down','Up'))

mus.res <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Liver') %>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

mus.txExp <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

mus.geneExp <-  read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(multiIso = ifelse(isoforms == 1, F, T)) %>% 
  mutate(overall_direction = ifelse(x3_t3_mean_exp_level < x3_t1_mean_exp_level,'Down','Up'))

temp <- adi.geneExp %>% 
  filter(name_description %in% adi.intersect$gene)

temp %>% group_by(overall_direction) %>% tally() # Majority of genes in intersect (93 out of 118) have lower expression in hibernation compared to active

plot_isoforms_and_gene <- function(gene_of_interest,tissue,use_color,use_labels=F) {
  
  res_file <- NULL
  exp_file <- NULL
  gene_exp_file <- NULL
  
  if (tissue == 'Adipose' | tissue == 'adipose') {
    res_file <- adi.res
    exp_file <- adi.txExp
    gene_exp_file <- adi.geneExp
  }
  
  if (tissue == 'Liver' | tissue == 'liver') {
    res_file <- liv.res
    exp_file <- liv.txExp
    gene_exp_file <-liv.geneExp
    
  }
  
  if (tissue == 'Muscle' | tissue == 'muscle') {
    res_file <- mus.res
    exp_file <- mus.txExp
    gene_exp_file <-mus.geneExp
    
  }
  
  checks <- res_file %>% 
    filter(gene_description == gene_of_interest)
  
  checks2 <- gene_of_interest %in% res_file$gene_description
  if (checks2 == F) {
    print(paste(gene_of_interest,' is not expressed in ',tissue,sep = ''))
    return()
  }
  
  if (checks$diu_result == 'Not DIU') {
    print(paste(gene_of_interest,' is not DIU',sep = ''))
    return()
  }
  
  switches <- res_file %>% 
    filter(gene_description == gene_of_interest) %>% 
    dplyr::select(switching_times) %>% 
    mutate(switch_names = case_when(
      switching_times=='1' ~ 'Active',
      switching_times=='12' ~ 'Active and Hyperphagia',
      switching_times=='13' ~ 'Active and Hibernation',
      switching_times=='23' ~ 'Hyperphagia and Hibernation',
      switching_times=='2' ~ 'Hyperphagia',
      switching_times=='3' ~ 'Hibernation',
      switching_times=='123' ~ 'all timepoints'
      
    ))
  
  switch_names <-  switches$switch_names
  
  color_scale <-  c('Major'='Black','Minor'='grey50')
  if (use_color==T) {
    color_scale <-  c('Major'='firebrick3','Minor'='grey50')
  }
  
  label_size <- 0
  if (use_labels==T) {
    label_size=4
  }
  plot.data <- exp_file %>%
    filter(gene_description == gene_of_interest) %>%
    dplyr::select(1,2,3,4,6,Active = 7, Hyperphagia = 8, Hibernation = 9) %>%
    pivot_longer(-c(1,2,3,4,5),names_to = 'timepoint',values_to = 'avg_exp') %>%
    mutate(timepoint = factor(timepoint,levels=c('Active','Hyperphagia','Hibernation'))) %>%
    group_by(number_transcript) %>%
    mutate(sum = sum(avg_exp)) %>%
    ungroup() %>%
    mutate(iso_type = ifelse(sum == max(sum),'Major','Minor'))
  
  tx_base <- paste(str_split_fixed(plot.data[1,1],'[.]',3)[,1],str_split_fixed(plot.data[1,1],'[.]',3)[,2],sep='.')
  
  gene.plot.data <- gene_exp_file %>% filter(name_description == gene_of_interest) %>%
    dplyr::select(2,Active = 4, Hyperphagia = 5, Hibernation = 6) %>%
    pivot_longer(-1,names_to = 'timepoint',values_to = 'avg_exp') %>%
    mutate(timepoint = factor(timepoint,levels=c('Active','Hyperphagia','Hibernation'))) 
  
  ggplot(plot.data,aes(x=timepoint,
                       y=avg_exp,
                       group=number_transcript,
                       lty=coding,
                       color=iso_type,
                       size=iso_type)) +
    geom_path() +
    geom_path(data=gene.plot.data,inherit.aes = F, 
              aes(x=timepoint,y=avg_exp,group=name_description),
              color='blue',lwd=2,alpha=0.5) +
    geom_point(size=1.5) +
    geom_point(data=gene.plot.data,inherit.aes = F, 
              aes(x=timepoint,y=avg_exp,group=name_description),
              color='blue',size=2.5,alpha=1) +
    ggrepel::geom_text_repel(data=subset(plot.data,timepoint=='Hibernation'),
                             aes(label=str_split_fixed(number_transcript,'[.]',3)[,3]),
                             size=label_size,nudge_x = 0.1,min.segment.length = 1) +
    scale_linetype_manual(values=c('YES'=1,'NO'=2)) +
    scale_color_manual(values=color_scale) +
    scale_size_manual(values=c('Major'=1,'Minor'=0.5)) +
    # geom_point(size=2) +
    labs(x=' ',
         y='Average Expression (TPM)',
         color='Isoform Type:',
         linetype = 'Coding:',
         size='Isoform Type:',
         title=paste(gene_of_interest, ' (',tx_base,')' ,' in ',str_to_sentence(tissue), '\nIsoform Switches at ',switch_names,sep='')) +
    theme_linedraw(base_size = 12) + 
    theme(plot.title.position = 'plot')
}


plot_isoforms_and_gene('GRB14',use_color=T,tissue='adipose',use_labels=T)
plot_isoforms_and_gene('TCIRG1',use_color=T,tissue='adipose',use_labels=T)

# AKT2 comparison
p.akt.adi <- plot_isoforms_and_gene('AKT2',use_color=T,tissue='adipose',use_labels=T) + scale_y_log10()
p.akt.mus <- plot_isoforms_and_gene('AKT2',use_color=T,tissue='muscle',use_labels=T) + scale_y_log10()

p.akt.adi + p.akt.mus + plot_layout(guides='collect')

plot_isoforms_and_gene('APPL1',use_color=T,tissue='liver',use_labels=T) + scale_y_log10()
plot_isoforms_and_gene('APPL1',use_color=T,tissue='muscle',use_labels=T) + scale_y_log10()

plot_isoforms_and_gene('SHC1',use_color = T,tissue = 'liver',use_labels = T)
plot_isoforms_and_gene('SORBS1',use_color = T,tissue = 'liver',use_labels = T)+ scale_y_log10()


plot_isoforms_and_gene('TCIRG1',use_color = T,tissue = 'adipose',use_labels = T)

# ggsave('_ms/figures/fig_pieces/Fig3_GRB14_ExpPanel_03.09.22.pdf')


plot_isoforms_and_gene('AKT2',use_color=T,tissue='adipose',use_labels=T) + scale_y_log10()

plot_isoforms_and_gene('ATF1',use_color=T,tissue='adipose',use_labels=T)


