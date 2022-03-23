
library(tidyverse)
library(patchwork)
library(ggupset)

# Read in Tappas DIU analysis results for each tissue ---------------------

adi.res <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Adipose') %>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

head(adi.res)

liv.res <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Liver')%>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

mus.res <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(tissue = 'Muscle')%>% 
  dplyr::select(1,2,3,4,5,6,t1_mean_exp_level=7,t2_mean_exp_level=8,t3_mean_exp_level=9,10)

# Combined into one full table --------------------------------------------

allTissue.res <- adi.res %>% 
  bind_rows(liv.res,mus.res)


# Define tissue color palette ---------------------------------------------

tissue_colors <- c('Adipose'='#E69F00','Liver'='#56B4E9','Muscle'='#009E73')

# Exploratory plotting ----------------------------------------------------

allTissue.res %>% 
  filter(diu_result=='DIU') %>% 
  group_by(tissue,major_isoform_switching) %>% 
  tally() %>% 
  ggplot(aes(x=tissue,y=n,fill=major_isoform_switching)) +  
  geom_text(aes(label=n),position = position_dodge(.9),vjust= -0.5,hjust=0.5) +
  geom_bar(stat='identity',position=position_dodge()) +
  labs(x='Tissue',y='Number of Genes',fill='Major Isoform Switching?',title = 'Number of DIU Genes w/ One or More Major Isoform Switches') +
  theme_linedraw() + theme(panel.grid = element_blank(),plot.title.position = 'plot')

allTissue.res %>% 
  filter(diu_result=='DIU') %>% 
  group_by(tissue,major_isoform_switching) %>% 
  tally() %>% 
  ggplot(aes(x=tissue,y=n,fill=major_isoform_switching)) +
  geom_bar(stat='identity',position = 'fill') +
  labs(x='Tissue',y='Percent of Genes',fill='Major Isoform Switching?',title = 'Percent of DIU Genes w/ One or More Major Isoform Switches') +
  theme_linedraw() + theme(panel.grid = element_blank(),plot.title.position = 'plot')


# Reformat table for GGUpset plots ----------------------------------------

allTissue.res.tidy <- allTissue.res %>% 
  filter(diu_result=='DIU') %>% 
  filter(major_isoform_switching=='YES') %>% 
  mutate(switching_times = 
           case_when(
             switching_times == 1 ~ 'Active',
             switching_times == 2 ~ 'Hyperphagia',
             switching_times == 3 ~ 'Hibernation',
             switching_times == 12 ~ 'Active;Hyperphagia',
             switching_times == 13 ~ 'Active;Hibernation',
             switching_times == 23 ~ 'Hyperphagia;Hibernation',
             switching_times == 123 ~ 'Active;Hyperphagia;Hibernation'
           )) %>% 
  mutate(switching_times = str_split(switching_times,';'))

p.upset <- allTissue.res.tidy %>%
  ggplot(aes(x=switching_times,fill=tissue)) +
  geom_bar(show.legend = F) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-.5,size=3) +
  scale_x_upset(sets = c('Active','Hyperphagia','Hibernation'),order_by = 'freq') +
  facet_wrap(~tissue,ncol=1) +
  scale_y_continuous(limits=c(0,450)) +
  scale_fill_manual(values=tissue_colors) +
  labs(x='Timing of Major Isoform Switches',y='\n\nNumber of Genes\n',
       title='Summary of Isoform Switching') +
  theme_linedraw() + theme(panel.grid = element_blank(),
                           plot.title.position = 'plot',
                           strip.background = element_rect(fill='grey22')) +
  theme_combmatrix(combmatrix.panel.point.color.fill = "grey22",
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.make_space = FALSE)

p.upset


# Read in adipose transcript expression table --------

adi.txExp <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

liv.txExp <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

mus.txExp <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names()

# Build plotting function 
plot_isoforms <- function(gene_of_interest,tissue,use_color,use_labels=F) {
  
  res_file <- NULL
  exp_file <- NULL
  
  if (tissue == 'Adipose' | tissue == 'adipose') {
    res_file <- adi.res
    exp_file <- adi.txExp
  }
  
  if (tissue == 'Liver' | tissue == 'liver') {
    res_file <- liv.res
    exp_file <- liv.txExp
  }
  
  if (tissue == 'Muscle' | tissue == 'muscle') {
    res_file <- mus.res
    exp_file <- mus.txExp
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
  
  ggplot(plot.data,aes(x=timepoint,
                       y=avg_exp,
                       group=number_transcript,
                       lty=coding,
                       color=iso_type,
                       size=iso_type)) +
    geom_path() +
    geom_point(size=1.5) +
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

plot_isoforms('ACYP1',use_color=T,tissue='adipose',use_labels=T)

plot_isoforms('ALDH4A1',use_color=T,tissue='adipose')

plot_isoforms('ITGB3BP',tissue='adipose',use_color=T)


p1 <- plot_isoforms('ACYP1',use_color=T,tissue='adipose',use_labels = T)
p2 <- plot_isoforms('ACCS',use_color=T,tissue='adipose',use_labels = T)
p3 <- plot_isoforms('ANKZF1',use_color=T,tissue='adipose',use_labels = T)

p123 <- p2 + p1 + p3 + plot_layout(guides='collect')

p.upset
p123



# Venn/comparison of DIU genes across tissues -----------------------------

# Broad level - genes with any DIU across tissues (ignoring switch times)

adi.res.diu <- adi.res %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')
liv.res.diu <- liv.res %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')
mus.res.diu <- mus.res %>% filter(diu_result=='DIU' & major_isoform_switching == 'YES')

vennSets <- list('Adipose'= adi.res.diu$number_gene,
                 'Liver'  = liv.res.diu$number_gene,
                 'Muscle' = mus.res.diu$number_gene)

ggvenn::ggvenn(vennSets,fill_color = c('#E69F00','#56B4E9','#009E73'))

all.tissue.intersect <- intersect(adi.res.diu$number_gene,intersect(liv.res.diu$number_gene,mus.res.diu$number_gene))

PPFIBP1.1 <- plot_isoforms('PPFIBP1',tissue='adipose',use_color = T,use_labels = T)
PPFIBP1.2 <- plot_isoforms('PPFIBP1',tissue='muscle',use_color = T,use_labels = T)
PPFIBP1.3 <- plot_isoforms('PPFIBP1',tissue='liver',use_color = T,use_labels = T)


SLC18B1.1 <- plot_isoforms('SLC18B1',tissue='adipose',use_color = T,use_labels = T)
SLC18B1.2 <- plot_isoforms('SLC18B1',tissue='muscle',use_color = T,use_labels = T)
SLC18B1.3 <- plot_isoforms('SLC18B1',tissue='liver',use_color = T,use_labels = T)


PPFIBP1.all <- PPFIBP1.1 + PPFIBP1.3 + PPFIBP1.2 + plot_layout(guides = 'collect')

SLC18B1.all <- SLC18B1.1 + SLC18B1.3 + SLC18B1.2 + plot_layout(guides = 'collect')

PPFIBP1.all / SLC18B1.all

### Venn of genes with timepoint-specific switches (more overlap between particular tissues at particular timepoint?)

adi.active.spec <- adi.res.diu %>% 
  filter(switching_times==1)
liv.active.spec <- liv.res.diu %>% 
  filter(switching_times==1)
mus.active.spec <- mus.res.diu %>% 
  filter(switching_times==1)

active.spec.vennSets <- list('Adipose'= adi.active.spec$number_gene,
                             'Liver'  = liv.active.spec$number_gene,
                             'Muscle' = mus.active.spec$number_gene)

ggvenn::ggvenn(active.spec.vennSets,fill_color = c('#E69F00','#56B4E9','#009E73')) + ggtitle('Acitve-specific switches')

adi.hyper.spec <- adi.res.diu %>% 
  filter(switching_times==2)
liv.hyper.spec <- liv.res.diu %>% 
  filter(switching_times==2)
mus.hyper.spec <- mus.res.diu %>% 
  filter(switching_times==2)

hyper.spec.vennSets <- list('Adipose'= adi.hyper.spec$number_gene,
                             'Liver'  = liv.hyper.spec$number_gene,
                             'Muscle' = mus.hyper.spec$number_gene)

ggvenn::ggvenn(hyper.spec.vennSets,fill_color = c('#E69F00','#56B4E9','#009E73')) + ggtitle('Hyperphagia-specific switches')

adi.hib.spec <- adi.res.diu %>% 
  filter(switching_times==3)
liv.hib.spec <- liv.res.diu %>% 
  filter(switching_times==3)
mus.hib.spec <- mus.res.diu %>% 
  filter(switching_times==3)

hib.spec.vennSets <- list('Adipose'= adi.hib.spec$number_gene,
                            'Liver'  = liv.hib.spec$number_gene,
                            'Muscle' = mus.hib.spec$number_gene)

ggvenn::ggvenn(hib.spec.vennSets,fill_color = c('#E69F00','#56B4E9','#009E73')) + ggtitle('Hibernation-specific switches')


# Subsetting relevant gene sets for GO and other characterization ---------

# Relevant sets:
# All DIU per tissue
# Active, hyper, hib specific switches for each tissue
# (HAVE NOT YET DONE THIS) Major isoform specific genes - i.e., genes where the major isoform is only highest expressed in a particular time points
# Background - all detected genes with more than 1 isoform in a given tissue
# Tissue-sepcific (based on venn)

adi.all.diu <- adi.res.diu$number_gene
liv.all.diu <- liv.res.diu$number_gene
mus.all.diu <- mus.res.diu$number_gene

adi.actSpec.geneList <- adi.active.spec$number_gene 
adi.hypSpec.geneList <- adi.hyper.spec$number_gene 
adi.hibSpec.geneList <- adi.hib.spec$number_gene 

liv.actSpec.geneList <- liv.active.spec$number_gene 
liv.hypSpec.geneList <- liv.hyper.spec$number_gene 
liv.hibSpec.geneList <- liv.hib.spec$number_gene 

mus.actSpec.geneList <- mus.active.spec$number_gene 
mus.hypSpec.geneList <- mus.hyper.spec$number_gene 
mus.hibSpec.geneList <- mus.hib.spec$number_gene 

adi.background <- adi.res$number_gene
liv.background <- liv.res$number_gene
mus.background <- mus.res$number_gene

adi.tissSpec <- adi.res.diu %>% 
  filter(!(number_gene %in% union(liv.res.diu$number_gene,mus.res.diu$number_gene)))
adi.tissSpec.geneList <- adi.tissSpec$number_gene

liv.tissSpec <- liv.res.diu %>% 
  filter(!(number_gene %in% union(adi.res.diu$number_gene,mus.res.diu$number_gene)))
liv.tissSpec.geneList <- liv.tissSpec$number_gene

mus.tissSpec <- mus.res.diu %>% 
  filter(!(number_gene %in% union(adi.res.diu$number_gene,liv.res.diu$number_gene)))
mus.tissSpec.geneList <- mus.tissSpec$number_gene

tissSpec.background <- union(adi.res.diu$number_gene,union(liv.res.diu$number_gene,mus.res.diu$number_gene))

write_tsv(as.data.frame(adi.all.diu),'analyses/tappas/_DIU_GeneSetsForGO/Adipose_AllDIUGenes.txt',col_names = F)
write_tsv(as.data.frame(liv.all.diu),'analyses/tappas/_DIU_GeneSetsForGO/Liver_AllDIUGenes.txt',col_names = F)
write_tsv(as.data.frame(mus.all.diu),'analyses/tappas/_DIU_GeneSetsForGO/Muscle_AllDIUGenes.txt',col_names = F)

write_tsv(as.data.frame(adi.actSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Adipose_Active_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(adi.hypSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Adipose_Hyperphagia_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(adi.hibSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Adipose_Hibernation_SpecificSwitches.txt',col_names = F)

write_tsv(as.data.frame(liv.actSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Liver_Active_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(liv.hypSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Liver_Hyperphagia_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(liv.hibSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Liver_Hibernation_SpecificSwitches.txt',col_names = F)

write_tsv(as.data.frame(mus.actSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Muscle_Active_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(mus.hypSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Muscle_Hyperphagia_SpecificSwitches.txt',col_names = F)
write_tsv(as.data.frame(mus.hibSpec.geneList),'analyses/tappas/_DIU_GeneSetsForGO/Muscle_Hibernation_SpecificSwitches.txt',col_names = F)

write_tsv(as.data.frame(adi.background),'analyses/tappas/_DIU_GeneSetsForGO/Adipose_Background1_AllMultiIsoGenes.txt',col_names = F)
write_tsv(as.data.frame(liv.background),'analyses/tappas/_DIU_GeneSetsForGO/Liver_Background1_AllMultiIsoGenes.txt',col_names = F)
write_tsv(as.data.frame(mus.background),'analyses/tappas/_DIU_GeneSetsForGO/Muscle_Background1_AllMultiIsoGenes.txt',col_names = F)

