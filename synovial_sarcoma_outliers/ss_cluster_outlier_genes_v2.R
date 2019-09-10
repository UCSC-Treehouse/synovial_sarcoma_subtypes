#Yvonne Vasquez
#Date: 7/15/19
#title: ss_cluster_outlier_genes_v1

library(readr)
library(dplyr)
library(tidyverse)
library(xlsx)

base_dir = "/Users/yvonnevasquez/Documents/GitHub/synovial_sarcoma_subtypes/synovial_sarcoma_outliers/"
setwd(base_dir)

path_to_data = "/Users/yvonnevasquez/Documents/GitHub/synovial_sarcoma_subtypes/synovial_sarcoma_outliers/syn_sarcoma_outliers/"

###
# sample info
###

fusion_results <- read_tsv(file.path(base_dir, "ss_sample_fusion_data - Sheet1.tsv"))

cluster_assignments <- read_tsv(file.path(base_dir,"synovial-hydra-assignments.tsv"))

sample_info <- full_join(fusion_results, cluster_assignments, by=c("sample_id"="sample"))


###
# gene info
###


#specifies data type found in each column
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

#combine all expression data files
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

#Create data table with all rows that have a gene that is an upoutlier for that sample
expression_data_filtered = filter(expression_data_raw,pc_outlier=="pc_up" | pd_outlier == "pd_up") 

#Create data table with only gene & sample_id columns
expression_data_filtered_genes <- select(expression_data_filtered, gene,pc_outlier, pd_outlier, sample_id)

## Add sample info
expression_data_filtered_genes_anno <- left_join(expression_data_filtered_genes, sample_info,  by="sample_id")

#Creates new table with all genes that are an up outlier in more than one sample and for which samples it is an up outlier
multi_sample_outliers <- expression_data_filtered_genes_anno %>%
  group_by(gene) %>%
  mutate(num_samples = n()) %>%
  filter(num_samples > 1)

outlier_pct_samples_present_in <- multi_sample_outliers %>%
  mutate(pct_samples_outlier_in = (num_samples/36)*100)

gene_summary_table <- multi_sample_outliers %>%
  group_by(gene, num_samples) %>%
  summarize(list_of_samples = paste(sample_id, collapse = ", ")) 


min_samples_with_outlier <- 10

outliers_to_plot <- multi_sample_outliers %>%
  filter(num_samples > min_samples_with_outlier) %>%
  arrange(desc(num_samples)) %>%
  ungroup %>%
  mutate(gene = factor(gene, levels = unique(gene)))
         

ggplot(outliers_to_plot) + geom_histogram(aes(x=gene), stat="count") +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  facet_wrap ( ~ SS18_fusion_partner, 
               ncol = 1)


ggplot(outliers_to_plot) + 
  stat_count(mapping = aes(x=gene, y=..prop.., group=1)) +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  facet_wrap ( ~ SS18_fusion_partner, 
               ncol = 1)

ggplot(outliers_to_plot) + geom_histogram(aes(x=gene), stat="prop") +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  facet_wrap ( ~ SS18_fusion_partner, 
               ncol = 1)

require(scales)
ggplot(outliers_to_plot, aes(x = gene)) +  
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels = percent) + 
  facet_wrap ( ~ SS18_fusion_partner, 
               ncol = 1) +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) 

ggplot(outliers_to_plot) + geom_histogram(aes(x=gene), stat="count") +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  facet_wrap ( ~ cluster, 
               ncol = 1)


         
  
ggplot(outliers_to_plot) + geom_point(aes(x=gene, y=num_samples)) +  
  theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  facet_grid (SS18_fusion_partner ~ . )
  
    
    ggplot(outliers_to_plot) + geom_point(aes(x=gene, y=num_samples)) +
    
  
  # show genes that are outliers in multiple TCGA samples
  # multi_sample_outliers %>% filter(grepl("TCGA", sample_id)) %>%
  #   group_by(gene) %>%
  #   summarize(n_tcga_samples_with_this_outlier = n())
  
  #fix this, paste(cluster)
  write_tsv(gene_summary_table, path=paste(cluster,'.tsv'))
  write.xlsx(x=gene_summary_table,file=paste(cluster,'.xlsx'))
  
  return(gene_summary_table)
  
})
