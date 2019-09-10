#Yvonne Vasquez
#Last Edited Date: 8/15/19
#title: gli1_expression_in_synovial_sarcoma_samples

library(readr)
library(dplyr)
library(tidyverse)
library(xlsx)

# base_dir = "/Users/yvonnevasquez/Desktop/Synovial_Sarcoma_Project/"
base_dir = "/Users/hbeale/Documents/Dropbox/ucsc/projects/gitCode/synovial_sarcoma_subtypes/"
setwd(base_dir)


# path_to_data = "/Users/yvonnevasquez/Desktop/Synovial_Sarcoma_Project/syn_sarcoma_outliers/"
path_to_data = "/Users/hbeale/Documents/Dropbox/ucsc/projects/gitCode/synovial_sarcoma_subtypes/synovial_sarcoma_outliers"
###
# sample info
###

#open fusion partner info
fusion_results <- read_tsv(file.path(base_dir, "ss_sample_fusion_data - Sheet1.tsv"))

#open hydra cluster assignmnets
cluster_assignments <- read_tsv(file.path(base_dir,"synovial-hydra-assignments.tsv"))

#open histology info
histological_types <- read_tsv(file.path(base_dir, "SS_sample_histological_type_info - Sheet1.tsv"))

#create data table with fusion partner & hydra cluster assignment info
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


#total number of synovial sarcoma genes
total_ss_genes <- 
  expression_data_raw %>% distinct(gene, .keep_all = TRUE)

expression_data_filtered_genes <- select(expression_data_raw, gene,pc_outlier, pd_outlier, sample_id, expression_in_log2tpm1, pc_high)

## Add sample info
expression_data_filtered_genes_anno <- left_join(expression_data_filtered_genes, sample_info_full,  by="sample_id")

#Add cancer gene and druggable gene info
cancer_gene_names <- scan(file = file.path(base_dir, "aggregatedCancerGenes_2018-01-04_12.20.15PM.txt"), what = 'list' )
druggable_gene_names <- read.table(file = file.path(base_dir, "treehouseDruggableGenes_2019-06-12.txt"), header=TRUE, sep='\t')

filtered_cancer_druggable_genes <- expression_data_filtered_genes_anno %>%
  mutate('Is_a_cancer_gene?' = gene %in% cancer_gene_names,  'Is a TH druggable gene?' = gene %in% druggable_gene_names$gene) 

gli1_expression_in_samples <- filtered_cancer_druggable_genes %>%
  filter(gene == "GLI1") %>%
  arrange(desc(expression_in_log2tpm1))


###
#  Plot Gli1 expression in synovial sarcoma samples
###

#bar plot showing Gli1 expression in synovial sarcoma samples
ggplot(gli1_expression_in_samples) + 
  geom_bar(aes(x=(reorder(sample_id, -expression_in_log2tpm1)), y=expression_in_log2tpm1), colour='black', width=1, stat='identity') +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  labs(title= "Gli1 expression in synovial sarcoma samples", x="synovial_sarcoma_sample_id") +
  geom_hline(yintercept = 3.123675, colour='red') + 
  geom_text(aes(0,3.123675,label = 'pan cancer up-outlier threshold = 3.123675', vjust = -1), size=3,position= position_nudge(x=6, y=0), colour='white')

ggsave("Gli1 expression in synovial sarcoma samples ID.png")


#histogram showing Gli1 expression levels and pc up-outlier cut-off
ggplot(gli1_expression_in_samples) +
  geom_histogram(aes(x=expression_in_log2tpm1), colour='black') +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  labs(title= "Number of synovial sarcoma samples with Gli1 expression", x="Gli1 expression in log2tpm1", y="num samples with expression value") +
  geom_vline(xintercept = 3.123675, colour='red') +
  geom_text(aes(3.123675,0, label = 'pan cancer up-outlier threshold = 3.123675'), position = position_nudge(x=0.9, y=4.2), colour='red')

ggsave("Number of synovial sarcoma samples with Gli1 expression ID.png")





  