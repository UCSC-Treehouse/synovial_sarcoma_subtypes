#Yvonne Vasquez
#Last Edited Date: 7/30/19
#title: ss_cluster_outlier_genes_v7

library(readr)
library(dplyr)
library(tidyverse)
library(xlsx)


base_dir = "/Users/yvonnevasquez/Desktop/Synovial_Sarcoma_Project/"
setwd(base_dir)


path_to_data = "/Users/yvonnevasquez/Desktop/Synovial_Sarcoma_Project/syn_sarcoma_outliers/"

###
# sample info
###

fusion_results <- read_tsv(file.path(base_dir, "ss_sample_fusion_data - Sheet1.tsv"))

cluster_assignments <- read_tsv(file.path(base_dir,"synovial-hydra-assignments.tsv"))

histological_types <- read_tsv(file.path(base_dir, "SS_sample_histological_type_info - Sheet1.tsv"))

#create data table with 
sample_info <- full_join(fusion_results, cluster_assignments, by=c("sample_id"="sample"))

#add histology info
sample_info_full <- full_join(sample_info, histological_types, by=c('sample_id'='th_sampleid')) 

#calculate how many samples are in each SS subgroup based on fusion partner (SSX1 vs SSX2 vs NA)
fusion_partner_sample_sizes <- sample_info_full %>% group_by(SS18_fusion_partner) %>% summarize(SS18_fusion_partner_group_size = n())

#calculate how many samples are in each SS subtype cluster 
cluster_sample_sizes <- sample_info_full %>% group_by(cluster) %>% summarize(cluster_group_size = n())


###
# gene info
###


#column specifications: specifies data type found in each column
col_spec=cols(
  Gene = col_character(),
  sample = col_double(),
  is_top_5 = col_character(),
  pc_low = col_double(),
  pc_median = col_double(),
  pc_high = col_double(),
  pc_outlier = col_character(),
  pc_is_filtered = col_character(),
  pd_low = col_character(),
  pd_median = col_character(),
  pd_high = col_character(),
  pd_outlier = col_character(),
  pc_percentile = col_integer()
)

#combine all SS sample expression data files
expression_data_raw<-lapply(sample_info_full$sample_id, function(sample_id) {
  # sample_id = sample_info_full$sample_id[1]
  file_name <-  file.path(path_to_data, paste0("outlier_results_", sample_id))
  if (file.exists(file_name)) {
    raw_outlier_candidates <- read_tsv(file_name, col_types=col_spec) %>%
      rename(expression_in_log2tpm1 =sample, gene=Gene) %>%
      mutate(
        sample_id=sample_id,
      ) 
  }
}) %>% bind_rows
#Filter for genes that are pan-cancer (pc) or pan-disease (pd) up-outliers
expression_data_filtered <- filter(expression_data_raw,pc_outlier=="pc_up" | pd_outlier == "pd_up") 

#total number of synovial sarcoma genes
total_ss_genes <- 
  expression_data_filtered %>% distinct(gene, .keep_all = TRUE)

expression_data_filtered_genes <- select(expression_data_filtered, gene,pc_outlier, pd_outlier, sample_id)

## Add sample info
expression_data_filtered_genes_anno <- left_join(expression_data_filtered_genes, sample_info_full,  by="sample_id")

#Add cancer gene and druggable gene info
cancer_gene_names <- scan(file = file.path(base_dir, "aggregatedCancerGenes_2018-01-04_12.20.15PM.txt"), what = 'list' )
druggable_gene_names <- read.table(file = file.path(base_dir, "treehouseDruggableGenes_2019-06-12.txt"), header=TRUE, sep='\t')

filtered_cancer_druggable_genes <- expression_data_filtered_genes_anno %>%
  mutate('Is_a_cancer_gene?' = gene %in% cancer_gene_names,  'Is a TH druggable gene?' = gene %in% druggable_gene_names$gene) 

#find number of total genes in each cluster
num_cluster_genes <- 
  filtered_cancer_druggable_genes %>% 
  group_by(cluster) %>%
  distinct(gene, .keep_all=TRUE) %>%
  summarize(n()) 

#find number of genes in each cluster are both pc and pd up outliers 
num_pc_and_pd_cluster_genes <-
  filtered_cancer_druggable_genes %>%
  filter(pd_outlier =='pd_up' & pc_outlier == 'pc_up') %>%
  group_by(cluster) %>%
  distinct(gene, .keep_all=TRUE) %>%
  summarize(num_pc_and_pd_genes=n()) 
  
#find number of genes in each cluster that are only pc up outliers 
num_pc_up_cluster_genes <-
  filtered_cancer_druggable_genes %>%
  filter(pc_outlier == 'pc_up' & is.na(pd_outlier)) %>%
  group_by(cluster) %>%
  distinct(gene, .keep_all=TRUE) %>%
  summarize(num_pc_up_genes = n()) 

#find number of genes in each cluster that are only pd up outliers
num_pd_up_cluster_genes <-
  filtered_cancer_druggable_genes %>%
  filter(pd_outlier == 'pd_up' & is.na(pc_outlier)) %>% 
  group_by(cluster) %>%
  distinct(gene, .keep_all=TRUE) %>%
  summarize(num_pd_up_genes = n()) 

#returns table showing how many genes are pc_only, pd_only, and both pc_and_pd in each SS cluster
cluster_gene_type_summary <-
  merge(num_pc_and_pd_cluster_genes, num_pd_up_cluster_genes) %>% 
  merge(num_pc_up_cluster_genes, by = 'cluster', all = TRUE) 

###top 10% pd-up genes
cluster_pd_up_only_genes <-
  
  filtered_cancer_druggable_genes %>%
  group_by(gene) %>% 
  mutate(num_both_pc_and_pd_samples = sum(pd_outlier =='pd_up' & pc_outlier == 'pc_up', NA, na.rm = TRUE),
         num_pc_only_samples = sum(pc_outlier == 'pc_up' & is.na(pd_outlier)),
         num_pd_only_samples = sum(pd_outlier == 'pd_up' & is.na(pc_outlier))) %>% 
  filter(pd_outlier == 'pd_up' & is.na(pc_outlier)) 

write_tsv(cluster_pd_up_genes, path = 'top_10_percent_cluster_pd_up_genes.tsv')


#Filter for genes that show up as both pc_up and pd_up outliers
genes_with_genes_both_pc_pd_up <-
  filtered_cancer_druggable_genes %>%
  group_by(gene) %>% 
  mutate(num_both_pc_and_pd_samples = sum(pd_outlier =='pd_up' & pc_outlier == 'pc_up', NA, na.rm = TRUE),
         num_pc_only_samples = sum(pc_outlier == 'pc_up' & is.na(pd_outlier)),
         num_pd_only_samples = sum(pd_outlier == 'pd_up' & is.na(pc_outlier))) %>% 
  filter(pd_outlier == 'pd_up' & pc_outlier == 'pc_up') 

#Filter for genes that show up as pd_up outliers
genes_with_pd_up_genes <-
  
  filtered_cancer_druggable_genes %>%
  group_by(gene) %>% 
  mutate(num_samples = sum(pd_outlier == 'pd_up', NA, na.rm = TRUE), 
         num_both_pc_and_pd_samples = sum(pd_outlier =='pd_up' & pc_outlier == 'pc_up', NA, na.rm = TRUE),
         num_pc_only_samples = sum(pc_outlier == 'pc_up' & is.na(pd_outlier)),
         num_pd_only_samples = sum(pd_outlier == 'pd_up' & is.na(pc_outlier))) %>% 
  filter(pd_outlier == 'pd_up') %>% 
  distinct(gene, .keep_all=TRUE) %>% 
  select(gene, num_samples, num_both_pc_and_pd_samples, num_pd_only_samples, num_pc_only_samples) %>% View()


#Filter pc_up outlier genes in syvnovial sarcoma samples
genes_with_pc_up_genes <-
  
  filtered_cancer_druggable_genes %>%
  group_by(gene) %>% 
  mutate(gene_is_outlier = sum(pc_outlier == 'pc_up', NA, na.rm = TRUE), 
         gene_is_pc_and_pd_outlier = sum(pd_outlier =='pd_up' & pc_outlier == 'pc_up', NA, na.rm = TRUE),
         gene_is_pc_outlier = sum(pc_outlier == 'pc_up' & is.na(pd_outlier)),
         gene_is_pd_outlier = sum(pd_outlier == 'pd_up' & is.na(pc_outlier))) %>% 
  filter(pc_outlier == 'pc_up') %>% 
  distinct(gene, .keep_all=TRUE) %>% 
  select(gene, gene_is_outlier, gene_is_pc_and_pd_outlier, gene_is_pc_outlier, gene_is_pd_outlier) %>% View()

#pc_up outlier genes in syvnovial sarcoma sample clusters
cluster_pc_outliers <-
  
  filtered_cancer_druggable_genes %>%
  filter(pc_outlier=='pc_up') %>% 
  group_by(gene, cluster) %>% 
  summarize(n_samples_with_outlier=n()) %>% View()







