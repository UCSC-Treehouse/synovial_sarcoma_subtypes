#Yvonne Vasquez
#Last Edited Date: 8/31/19
#title: ss_cluster_outlier_genes_v9

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


###pc_up outlier in synovial sarcoma samples analysis

###shows pc_up genes in cluster 1 with different fusions
cluster1_pc_outliers <-
  
  filtered_cancer_druggable_genes %>%
  filter(cluster == '1', pc_outlier=='pc_up') %>% 
  group_by(gene) %>%
  mutate(gene_is_pc_outlier = n()) %>%
  filter(gene_is_pc_outlier >1) %>% 
  group_by(gene, SS18_fusion_partner) %>%
  mutate(samples_with_fusion = n()) %>%
  select(gene, sample_id, SS18_fusion_partner) %>% 
  group_by(gene) %>%
  summarize(
    n_samples_with_outlier = n(),
    number_SSX1_samples = sum(SS18_fusion_partner == "SSX1"),
    pct_SSX1_samples = round(100*number_SSX1_samples/n())
    )
  


pc_up_outliers <- filtered_cancer_druggable_genes %>%
  select(-pd_outlier) %>%
  filter(pc_outlier == "pc_up")

#pc_up outliers in greater than 10% of cluster samples
pc_up_outliers_gt_10pct_cluster_samples <- 
  pc_up_outliers %>%
  left_join(cluster_sample_sizes) %>%
  group_by(gene, cluster) %>% 
  mutate(
    n_samples_in_cluster_with_upoutlier = n(),
    pct_samples_in_cluster_with_upoutlier = n_samples_in_cluster_with_upoutlier/cluster_group_size
  ) %>%
  filter(pct_samples_in_cluster_with_upoutlier>0.1,
         n_samples_in_cluster_with_upoutlier > 1) 



# present in one cluster only 
#this finds genes unique to each cluster, can be changed to find genes in 2 or 3 clusters
#present in 10 % of cluster samples and >1 sample
pc_up_outliers_in_one_cluster_gt_10pct <- 
  
  pc_up_outliers_gt_10pct_cluster_samples %>%
  group_by(gene, cluster) %>%
  summarize(n_outlier_instances_in_all_clusters=n()) %>% 
  group_by(gene) %>%
  mutate(n_clusters_with_outlier = n()) %>% 
  filter(n_clusters_with_outlier == 1) %>% #outliers in 1 cluster, can change to be in 2 clusters, or 3
  rename(n_samples_in_cluster_with_outlier = n_outlier_instances_in_all_clusters) %>%
  arrange(desc(n_samples_in_cluster_with_outlier)) 

write_tsv(pc_up_outliers_in_one_cluster_gt_10pct, path='pc_up_outliers_in_one_cluster_gt_10pct.tsv')



###add sample info to pc_up_outliers_in_one_cluster_gt_10pct
pc_outliers_in_one_cluster_gt_10pct_sample_info <-
  
  pc_up_outliers_in_one_cluster_gt_10pct %>%
  left_join(pc_up_outliers_gt_10pct_cluster_samples) %>% 
  select(-n_samples_in_cluster_with_upoutlier, 
         -pc_outlier) %>% 
  arrange(desc(n_samples_in_cluster_with_outlier))

write_tsv(pc_outliers_in_one_cluster_gt_10pct_sample_info, path = 'pc_outliers_in_one_cluster_gt_10pct_sample_info.tsv')


# omit samples without fusion info
#calculate fusion occurrence for each outlier
fusion_info_for_pc_outliers_in_one_cluster<-
  
  pc_outliers_in_one_cluster_gt_10pct_sample_info %>%
  filter(! is.na(SS18_fusion_partner)) %>%
  group_by(gene, cluster, n_samples_in_cluster_with_outlier) %>%
  summarize(
  number_SSX1_samples = sum(SS18_fusion_partner == "SSX1"),
pct_SSX1_samples = round(100*number_SSX1_samples/n())
) %>% arrange(desc(n_samples_in_cluster_with_outlier)) 
write_tsv(fusion_info_for_pc_outliers_in_one_cluster, path = 'fusion_info_for_pc_outliers_in_two_clusters.tsv')
  

# make a histogram
# how many samples a gene appears in
# and facet by cluster

plotme1 <- pc_outliers_in_one_cluster_gt_10pct_sample_info %>%
  left_join(cluster_sample_sizes, by = "cluster") %>% 
  mutate(cluster_label = paste0(cluster, ", n=", cluster_group_size.y))
  
outliers_by_cluster <- plotme1 %>%
  group_by(gene, cluster, cluster_label) %>%
  summarize(n_outliers = n()) %>% View()

ggplot(plotme1) + 
  geom_histogram(aes(x=n_samples_in_cluster_with_outlier), binwidth=0.5) + 
  facet_wrap (~ cluster_label, ncol = 1) +
  labs(y='n_outliers')  

ggsave("sample outlier count based on cluster ID.png")
