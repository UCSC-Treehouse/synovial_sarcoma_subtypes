#Yvonne Vasquez
#Date: 7/30/19
#title: ss_cluster_outlier_genes_v6

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

#create data table with 
sample_info <- full_join(fusion_results, cluster_assignments, by=c("sample_id"="sample"))

# sample tile plot
# ggplot(sample_info) + geom_tile(aes(x=sample_id, y=cluster, fill=SS18_fusion_partner)) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

#calculate how many samples are in each SS subgroup based on fusion partner (SSX1 vs SSX2 vs NA)
fusion_partner_sample_sizes <- sample_info %>% group_by(SS18_fusion_partner) %>% summarize(SS18_fusion_partner_group_size = n())

#calculate how many samples are in each SS subtype cluster 
cluster_sample_sizes <- sample_info %>% group_by(cluster) %>% summarize(cluster_group_size = n())


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
expression_data_raw<-lapply(sample_info$sample_id, function(sample_id) {
  # sample_id = sample_info$sample_id[1]
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
expression_data_filtered = filter(expression_data_raw,pc_outlier=="pc_up" | pd_outlier == "pd_up") 

expression_data_filtered_genes <- select(expression_data_filtered, gene,pc_outlier, pd_outlier, sample_id)

## Add sample info
expression_data_filtered_genes_anno <- left_join(expression_data_filtered_genes, sample_info,  by="sample_id")

#Add cancer gene and druggable gene info
cancer_gene_names <- scan(file = file.path(base_dir, "aggregatedCancerGenes_2018-01-04_12.20.15PM.txt"), what = 'list' )
druggable_gene_names <- read.table(file = file.path(base_dir, "treehouseDruggableGenes_2019-06-12.txt"), header=TRUE, sep='\t')

filtered_cancer_druggable_genes <- expression_data_filtered_genes_anno %>%
  mutate('Is_a_cancer_gene?' = gene %in% cancer_gene_names,  'Is a TH druggable gene?' = gene %in% druggable_gene_names$gene) 


#Table of cluster specific outlier gene info
cluster_sample_outliers <- 
  
  filtered_cancer_druggable_genes %>%
  group_by(gene, cluster) %>% 
  summarize(n_outlier_instances_in_all_clusters=n()) %>% 
  group_by(gene) %>%
  mutate(n_clusters_with_outlier = n()) %>% View()
  filter(n_clusters_with_outlier == 1) %>% #outliers in 1 cluster, can change to be in 2 clusters, or 3
  arrange(desc(n_outlier_instances_in_all_clusters)) %>% 
  left_join(cluster_sample_sizes, by = "cluster") %>% 
  mutate(percent_cluster_with_outlier = (n_outlier_instances/cluster_group_size)*100) %>%View()

#
fusion_partner_sample_outliers <- 
  
  
  filtered_cancer_druggable_genes %>%
  filter(!is.na(SS18_fusion_partner)) %>% 
  group_by(gene, SS18_fusion_partner) %>% 
  summarize(n_outlier_instances_in_all_clusters=n()) %>% 
  group_by(gene) %>% 
  mutate(n_clusters_with_outlier = n()) %>% 
  filter(n_clusters_with_outlier == 1) %>% #outliers in 1 cluster, can change to be in 2 clusters, or 3
  arrange(desc(n_outlier_instances_in_all_clusters)) %>%  
  left_join(fusion_partner_sample_sizes, by = "SS18_fusion_partner") %>% View()
  mutate(percent_fusion_with_outlier = (n_outlier_instances/cluster_group_size)*100) %>% View()






