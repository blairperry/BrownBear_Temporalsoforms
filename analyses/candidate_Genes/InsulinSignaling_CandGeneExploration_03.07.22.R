
library(tidyverse)
library(colorspace)

insr.candList <- read_tsv('analyses/candidate_Genes/InsulinSignaling_Uniprot_CandGenes_AllSymbols.txt',col_names = 'gene')

tissue_colors <- c('Adipose'='#E69F00','Liver'='#56B4E9','Muscle'='#009E73')


# Intersect with ALL genes in transcriptome -------------------------------

gff.geneInfo <- read_tsv('data/new_merge.combined.renamed_tappAS_annot_from_SQANTI3.geneEntriesOnly.txt',col_names = F) %>% 
  dplyr::select(desc = 9) %>% 
  mutate(desc = str_split_fixed(desc,'[;]',2)[,1] %>% str_remove_all('ID=')) %>% 
  mutate(present_in_bear = T) %>% 
  unique()

insr.candList.full <- read_tsv('analyses/candidate_Genes/insulinReceptorSignalingPathway_Uniprot_10.20.21.tsv') %>% 
  janitor::clean_names() %>% 
  select(entry_name,gene_names) %>% 
  mutate(entry_name = str_remove_all(entry_name,'_HUMAN')) %>% 
  mutate(gene_names = str_replace_all(gene_names,' ',';')) %>% 
  mutate(gene_names = paste(entry_name,';',gene_names,sep = '')) %>% 
  # mutate(gene_names = str_split(gene_names,';')) %>% 
  separate_rows(gene_names,sep = ';') %>% 
  dplyr::select(synonym = gene_names, main_gene = entry_name)
  
insr.candList.full.overlap <- insr.candList.full %>% 
  left_join(gff.geneInfo,by=c('synonym' = 'desc')) %>% 
  filter(present_in_bear == T) %>% 
  unique()

sum(insr.candList.full.overlap$main_gene %in% insr.candList.full$main_gene) / length(unique(insr.candList.full$main_gene))








# Intersect w/ DIU results ------------------------------------------------

adi.diu.insr <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  # filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  filter(diu_result == 'DIU') %>%
  mutate(insrGene = ifelse(number_gene %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Adipose',result='Differential isoform usage')

liv.diu.insr <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  # filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  filter(diu_result == 'DIU') %>%
  mutate(insrGene = ifelse(number_gene %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Liver',result='Differential isoform usage')

mus.diu.insr <- read_tsv('analyses/tappas/muscle_res/Muscle_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  # filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  filter(diu_result == 'DIU') %>%
  mutate(insrGene = ifelse(number_gene %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Muscle',result='Differential isoform usage')

all.diu.insr <- adi.diu.insr %>% 
  bind_rows(liv.diu.insr,mus.diu.insr) %>% 
  dplyr::select(gene=number_gene,tissue,insrGene,switching_times) %>% 
  filter(insrGene==T) %>% 
  arrange(gene) %>% 
  mutate(gene=factor(gene,levels=rev(unique(.$gene))),
         switching_times = 
                  case_when(
                    is.na(switching_times) ~ 'None',
                    switching_times == 1 ~ 'Active',
                    switching_times == 2 ~ 'Hyperphagia',
                    switching_times == 3 ~ 'Hibernation',
                    switching_times == 12 ~ 'Active and Hyperphagia',
                    switching_times == 13 ~ 'Active and Hibernation',
                    switching_times == 23 ~ 'Hyperphagia and Hibernation',
                    switching_times == 123 ~ 'All timepoints'
                  )) %>% 
  mutate(switching_times = factor(switching_times,levels=c('Active','Hyperphagia','Hibernation','Active and Hyperphagia','Active and Hibernation','Hyperphagia and Hibernation','All timepoints','None')))


ggplot(all.diu.insr,aes(x=tissue,y=gene)) +
  # geom_point(aes(fill=switching_times),pch=23,size=6,show.legend = T) +
  geom_point(aes(pch=switching_times),size=6,show.legend = T) +
  # scale_fill_discrete_qualitative('Dark3') +  
  labs(x='Tissue',y='Candidate Gene',title='Insulin Signaling Candidate Genes: Significantly DIU',fill='Switching times:') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')
# ggsave('analyses/candidate_Genes/Fig3_InsrDIUGenes_03.08.22.pdf')

# Intersect w/ DEA results ------------------------------------------------

adi.dea.insr <- read.table('analyses/tappas/adipose_res/Adipose_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(insrGene = ifelse(gene_id %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Adipose',result='Differential gene expression')

adi.dea.insrONLY <- adi.dea.insr %>% filter(insrGene==T)

liv.dea.insr <- read.table('analyses/tappas/liver_res/Liver_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(insrGene = ifelse(gene_id %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Liver',result='Differential gene expression')

liv.dea.insrONLY <- liv.dea.insr %>% filter(insrGene==T)

mus.dea.insr <- read.table('analyses/tappas/muscle_res/Muscle_DEA_result_gene.tsv') %>% 
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  mutate(insrGene = ifelse(gene_id %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Muscle',result='Differential gene expression')

mus.dea.insrONLY <- mus.dea.insr %>% filter(insrGene==T)

all.dea.insr <- adi.dea.insr %>% 
  bind_rows(liv.dea.insr,mus.dea.insr) %>% 
  dplyr::select(gene=gene_id,tissue,insrGene) %>% 
  filter(insrGene==T) %>% 
  arrange(gene) %>% 
  mutate(gene=factor(gene,levels=rev(unique(.$gene))))

# ggplot(all.dea.insr,aes(x=tissue,y=gene,fill=tissue)) +
#   geom_point(inherit.aes = F, aes(x='Muscle',y='INSR'),alpha=0) +
#   geom_point(pch=23,size=6,show.legend = F) +
#   scale_fill_manual(values=tissue_colors) +
#   labs(x='',y='Candidate Gene',title='Insulin Signaling Candidate Genes: Significantly DE') +
#   theme_linedraw() + theme(plot.title.position = 'plot')

# Plot gene expression profiles

adi.geneInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  filter(number_gene %in% adi.dea.insrONLY$gene_id) %>% 
  dplyr::select(gene=1,Active=4,Hyperphagia=5,Hibernation=6,Active2 = 4) %>% 
  pivot_longer(-1,names_to = 'season',values_to = 'meanExp') %>% 
  mutate(season = factor(season,levels=c('Active','Hyperphagia','Hibernation','Active2')),
         tissue='Adipose')

liv.geneInfo <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  filter(number_gene %in% liv.dea.insrONLY$gene_id) %>% 
  dplyr::select(gene=1,Active=4,Hyperphagia=5,Hibernation=6,Active2 = 4) %>% 
  pivot_longer(-1,names_to = 'season',values_to = 'meanExp') %>% 
  mutate(season = factor(season,levels=c('Active','Hyperphagia','Hibernation','Active2')),
         tissue='Liver')

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
} # function to scale exp data within groups

all.geneInfo <- adi.geneInfo %>% 
  bind_rows(liv.geneInfo) %>% 
  group_by(gene,tissue) %>% 
  mutate(scaledExp = scale_this(meanExp))

ggplot(all.geneInfo,aes(x=season,y=meanExp,group=gene,color=gene)) +
  geom_path(alpha=0) +
  geom_path(data=subset(all.geneInfo,season != 'Active2'),lwd=1,show.legend = F) +
  geom_path(data=subset(all.geneInfo,season %in% c('Hibernation','Active2')),lty=2,show.legend = F,lwd=1) +
  ggrepel::geom_text_repel(data=subset(all.geneInfo,season=='Active'),aes(label=gene),size=3,nudge_x = -2) +
  facet_wrap(~tissue,ncol=2) +
  scale_y_log10() +
  scale_x_discrete(labels=c('Active','Hyperphagia','Hibernation','Active')) +
  labs(x='',y='Mean Gene-Level Expression',title='Insuling Signaling: Significantly DE Genes') +
  theme_linedraw()

# ggsave('analyses/candidate_Genes/Fig3_InsrGeneDEAProfiles_03.08.22.pdf')

# Read in data to generate heatmap of candidate genes

adi.txInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(tx_id = 1, gene_id = 3)

adi.exp <- read.table('analyses/tappas/adipose_res/Adipose_NormTranscriptMatrix.tsv') %>% 
  as.data.frame() %>% 
  rownames_to_column('tx_id') %>% 
  left_join(adi.txInfo,by='tx_id') %>% 
  filter(!is.na(gene_id)) %>% 
  filter(gene_id %in% adi.geneInfo$gene) %>% 
  dplyr::select(tx_id,gene_id,everything())

adi.exp.byGene <- adi.exp %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'norm_count') %>% 
  group_by(gene_id,sample) %>% 
  summarise(gene_count = sum(norm_count)) %>% 
  mutate(time = case_when(
    str_detect(sample,'Hy') ~ 'Hyperphagia',
    str_detect(sample,'Hi') ~ 'Hibernation',
    str_detect(sample,'A') ~ 'Active'
  )) %>% 
  mutate(time = factor(time, levels=c('Active','Hyperphagia','Hibernation'))) %>% 
  arrange(time) %>% 
  dplyr::select(-time) %>% 
  pivot_wider(names_from=sample,values_from = gene_count) %>% 
  column_to_rownames('gene_id')

pheatmap::pheatmap(adi.exp.byGene,
                   scale='row',
                   cellheight = 15,cellwidth = 15,
                   cluster_cols = F,gaps_col = c(5,11),
                   cutree_rows = 3,
                   border_color = 'black',
                   # color=diverge_hcl(n=50,palette = 'Blue-Red2')
                   # color=sequential_hcl(n=50,palette = 'Oslo',fixup = T,rev = T),
                   color=viridis::viridis(50),
                   # filename = '_ms/figures/fig_pieces/Fig3_InsrDEAgenes_AdipHeat_03.15.22.pdf'
                   )


liv.txInfo <- read_tsv('analyses/tappas/liver_res/Liver_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(tx_id = 1, gene_id = 3)

liv.exp <- read.table('analyses/tappas/liver_res/Liver_NormTranscriptMatrix.tsv') %>% 
  as.data.frame() %>% 
  rownames_to_column('tx_id') %>% 
  left_join(liv.txInfo,by='tx_id') %>% 
  filter(!is.na(gene_id)) %>% 
  filter(gene_id %in% liv.geneInfo$gene) %>% 
  dplyr::select(tx_id,gene_id,everything())

liv.exp.byGene <- liv.exp %>% 
  pivot_longer(-c(1,2),names_to = 'sample',values_to = 'norm_count') %>% 
  group_by(gene_id,sample) %>% 
  summarise(gene_count = sum(norm_count)) %>% 
  mutate(time = case_when(
    str_detect(sample,'Hy') ~ 'Hyperphagia',
    str_detect(sample,'Hi') ~ 'Hibernation',
    str_detect(sample,'A') ~ 'Active'
  )) %>% 
  mutate(time = factor(time, levels=c('Active','Hyperphagia','Hibernation'))) %>% 
  arrange(time) %>% 
  dplyr::select(-time) %>% 
  pivot_wider(names_from=sample,values_from = gene_count) %>% 
  column_to_rownames('gene_id')

pheatmap::pheatmap(liv.exp.byGene,
                   scale='row',
                   cellheight = 15,cellwidth = 15,
                   cluster_cols = F,gaps_col = c(5,11),
                   cutree_rows = 3,
                   border_color = 'black',
                   # color=diverge_hcl(n=50,palette = 'Blue-Red2')
                   # color=sequential_hcl(n=50,palette = 'Oslo',fixup = T,rev = T),
                   color=viridis::viridis(50),
                   # filename = '_ms/figures/fig_pieces/Fig3_InsrDEAgenes_LivpHeat_03.15.22.pdf'
)







