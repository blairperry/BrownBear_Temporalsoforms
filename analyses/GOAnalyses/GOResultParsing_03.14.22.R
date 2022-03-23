
library(tidyverse)
# BiocManager::install("rrvgo")
library(rrvgo)

# Adipose results ---------------------------------------------------------

adi.diuAndDE.bp.res <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-BP-GO.txt',skip = 12,
                             col_names = c('go_term','num_ref','num_target','expected','over_under','fold_enrich','p','fdr')) %>% 
  mutate(category='biological_process')

adi.diuAndDE.mf.res <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-MF-GO.txt',skip = 12,
                                col_names = c('go_term','num_ref','num_target','expected','over_under','fold_enrich','p','fdr')) %>% 
  mutate(category='molecular_function')

adi.diuAndDE.cc.res <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-CC-GO.txt',skip = 12,
                                col_names = c('go_term','num_ref','num_target','expected','over_under','fold_enrich','p','fdr')) %>% 
  mutate(category='cellular_component')

adi.diuAndDE.res.all <- adi.diuAndDE.cc.res %>% 
  bind_rows(adi.diuAndDE.bp.res,adi.diuAndDE.mf.res) %>% 
  filter(str_detect(fold_enrich, '<', negate = T)) %>% 
  mutate(fold_enrich = as.numeric(fold_enrich))

# Quick barplot 
ggplot(adi.diuAndDE.res.all %>% filter(fdr < 0.05),aes(fill=category,x=fold_enrich,y=reorder(go_term,fold_enrich))) +
  geom_bar(stat='identity',orientation = 'y',color='black') +
  scale_fill_manual(values=c('cellular_component'='grey','biological_process'='grey15')) +
  labs(x='Fold enrichment',y='GO Term',fill='GO Category:') +
  theme_linedraw()


# Format for GO term simplication w/ Rrvgo (GO term and FDR)

adi.diuAndDE.res.simple <- adi.diuAndDE.res.all %>% 
  filter(fdr < 0.05) %>% 
  mutate(go_id = str_split_fixed(go_term,'[(]',2)[,2] %>% str_remove_all('[)]')) %>% 
  select(go_id,fdr) %>% 
  mutate(score = -log10(fdr))

# Calculate similarity between terms
simMatrix.bp <- calculateSimMatrix(adi.diuAndDE.res.simple$go_id,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                # ont='CC',
                                method="Rel")

# Reduce to less redundant parent terms
scores.bp <- setNames(adi.diuAndDE.res.simple$score,adi.diuAndDE.res.simple$go_id)

reducedTerms.bp <- reduceSimMatrix(simMatrix, 
                                scores, # group representative is that with highest score (here -log10(q))
                                threshold = 0.7,
                                orgdb="org.Hs.eg.db")


simMatrix.cc <- calculateSimMatrix(adi.diuAndDE.res.simple$go_id,
                                   orgdb="org.Hs.eg.db",
                                   # ont="BP",
                                   ont='CC',
                                   method="Rel")

# Reduce to less redundant parent terms
scores.cc <- setNames(adi.diuAndDE.res.simple$score,adi.diuAndDE.res.simple$go_id)

reducedTerms.cc <- reduceSimMatrix(simMatrix.cc, 
                                   scores, # group representative is that with highest score (here -log10(q))
                                   threshold = 0.7,
                                   orgdb="org.Hs.eg.db")


reducedTerms.all <- reducedTerms.bp %>% 
  bind_rows(reducedTerms.cc)

# Barplot colored by rrvgo terms

adi.diuAndDE.res.all %>% 
  filter(fdr < 0.05) %>% 
  mutate(term = str_split_fixed(go_term,'[(]',2)[,1] %>% trimws(which = 'both')) %>%  
  left_join(reducedTerms.all,by='term') %>% 
  ggplot(aes(fill=parentTerm,x=fold_enrich,y=reorder(go_term,fold_enrich))) +
  geom_bar(stat='identity',orientation = 'y') +
  # scale_fill_manual(values=c('cellular_component'='grey','biological_process'='grey15')) +
  labs(x='Fold enrichment',y='GO Term',fill='GO Category:') +
  theme_linedraw()



# Isoform switch and overall gene expression for GO genes -----------------
library(ggupset)

# read in adipose DIU+DE GO results

adi.diuAndDE.bp.res <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-BP-GO.txt',skip = 12,
                                col_names = c('go_term','num_ref','num_target','expected','over_under','fold_enrich','p','fdr')) %>% 
  mutate(category='biological_process')

adi.diuAndDE.cc.res <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-CC-GO.txt',skip = 12,
                                col_names = c('go_term','num_ref','num_target','expected','over_under','fold_enrich','p','fdr')) %>% 
  mutate(category='cellular_component')

adi.diuAndDE.res.sig <- adi.diuAndDE.cc.res %>% 
  bind_rows(adi.diuAndDE.bp.res) %>% 
  filter(str_detect(fold_enrich, '<', negate = T)) %>% 
  mutate(fold_enrich = as.numeric(fold_enrich)) %>% 
  filter(fdr < 0.05) %>% 
  mutate(term = str_split_fixed(go_term,'[(]',2)[,1] %>% trimws(which = 'both')) %>% 
  mutate(go_id = str_split_fixed(go_term,'[(]',2)[,2] %>% str_remove_all('[)]'))


go.geneAndTerm <- read_tsv('analyses/GOAnalyses/Go-Results-Full-Output/adipose/Adipose-DeandDiu-GO_GenesAndTerms.txt',col_names = c('geneID','term')) %>% 
  separate_rows(geneID,sep=',') %>% 
  mutate(term = trimws(term,which = 'both')) %>% 
  filter(term %in% adi.diuAndDE.res.sig$term)

adi.diu <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_DIUGene_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  filter(diu_result == 'DIU' & major_isoform_switching=='YES') %>%
  mutate(tissue = 'Adipose',result='Differential isoform usage')

adi.diu.goSet <- adi.diu %>% 
  filter(number_gene %in% go.geneAndTerm$geneID)

adi.diu.goSet.tidy <- adi.diu.goSet %>% 
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

adi.diu.goSet.tidy %>%
  ggplot(aes(x=switching_times)) +
  geom_bar(show.legend = F) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-.5,size=3) +
  scale_x_upset(sets = c('Active','Hyperphagia','Hibernation'),order_by = 'freq') +
  # facet_wrap(~tissue,ncol=1) +
  # scale_y_continuous(limits=c(0,450)) +
  # scale_fill_manual(values=tissue_colors) +
  labs(x='Timing of Major Isoform Switches',y='\n\nNumber of Genes\n',
       title='Summary of Isoform Switching') +
  theme_linedraw() + theme(panel.grid = element_blank(),
                           plot.title.position = 'plot',
                           strip.background = element_rect(fill='grey22')) +
  theme_combmatrix(combmatrix.panel.point.color.fill = "grey22",
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.make_space = FALSE)



# Gene heatmap
adi.txInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_Transcripts.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(tx_id = 1, gene_id = 3)

adi.geneInfo <- read_tsv('analyses/tappas/adipose_res/Adipose_tappAS_GenesInfo.tsv') %>% 
  janitor::clean_names() %>% 
  filter(number_gene %in% adi.diu.goSet$number_gene) %>% 
  dplyr::select(gene=1,Active=4,Hyperphagia=5,Hibernation=6,Active2 = 4) %>% 
  pivot_longer(-1,names_to = 'season',values_to = 'meanExp') %>% 
  mutate(season = factor(season,levels=c('Active','Hyperphagia','Hibernation','Active2')),
         tissue='Adipose')

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

adi.diu.goSet.annot <- go.geneAndTerm %>% 
  left_join(reducedTerms.all,by='term') %>% 
  select(geneID,parentTerm) %>% 
  mutate(value = 1) %>% 
  unique() %>% 
  pivot_wider(names_from = parentTerm,values_from = value,values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames('geneID')


pheatmap::pheatmap(adi.exp.byGene,
                   scale='row',
                   # cellheight = 15,cellwidth = 15,
                   annotation_row = adi.diu.goSet.annot,
                   cluster_cols = F,gaps_col = c(5,11),
                   cutree_rows = 3,
                   border_color = 'black',
                   color=viridis::viridis(50),
                   # filename = '_ms/figures/fig_pieces/Fig3_InsrDEAgenes_AdipHeat_03.15.22.pdf'
)








