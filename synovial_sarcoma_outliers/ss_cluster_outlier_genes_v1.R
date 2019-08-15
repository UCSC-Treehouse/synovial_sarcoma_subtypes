#Yvonne Vasquez
#Date: 7/15/19
#title: ss_cluster_outlier_genes_v1

library(readr)
library(dplyr)
library(tidyverse)
library(xlsx)

path_to_data = "/Users/yvonnevasquez/Desktop/Synovial_Sarcoma_Project/syn_sarcoma_outliers/"

# make a list of the three synovial sarcoma (SS) cluster expression files
cluster_list <- c('SS_cluster1', 'SS_cluster2', 'SS_cluster3')
#Contains up outlier gene summary for each SS cluster
gene_summary_table_all <- lapply(cluster_list, function(cluster) {
# cluster = cluster_list[1]
  outlier_results_files<-tibble(file_name=list.files(paste0(path_to_data, cluster)))
  
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
  expression_data_raw<-lapply(outlier_results_files$file_name, function(file_name) {
    this_sample_id=gsub("^.*outlier_results_", "", file_name)
    raw_outlier_candidates <- read_tsv(paste(path_to_data, cluster, file_name, sep="/"), col_types=col_spec) %>%
      rename(expression_in_log2tpm1 =sample, gene=Gene) %>%
      mutate(
        sample_id=this_sample_id,
        
      ) 
  }) %>% bind_rows
  
  #Create data table with all rows that have a gene that is an upoutlier for that sample
  expression_data_filtered = filter(expression_data_raw,pc_outlier=="pc_up")
  
  #Create data table with only gene & sample_id columns
  expression_data_filtered_genes <- select(expression_data_filtered, gene,sample_id)
  
  #order data table by gene
  ordered_genes <- arrange(expression_data_filtered_genes, gene)

  #Creates new table with all genes that are an up outlier in more than one sample and for which samples it is an up outlier
  gene_summary_table <- ordered_genes %>%
    group_by(gene) %>%
    summarize(num_samples = n(), list_of_samples = paste(sample_id, collapse = ", ")) %>%
    filter(num_samples > 1)
  
  multi_sample_outliers <- ordered_genes %>%
    group_by(gene) %>%
    mutate(num_samples = n()) %>%
    filter(num_samples > 1)
  
  gene_summary_table <- multi_sample_outliers %>%
    group_by(gene, num_samples) %>%
    summarize(list_of_samples = paste(sample_id, collapse = ", ")) 
    
  
# show genes that are outliers in multiple TCGA samples
  # multi_sample_outliers %>% filter(grepl("TCGA", sample_id)) %>%
  #   group_by(gene) %>%
  #   summarize(n_tcga_samples_with_this_outlier = n())
  
    #fix this, paste(cluster)
    write_tsv(gene_summary_table, path=paste(cluster,'.tsv'))
    write.xlsx(x=gene_summary_table,file=paste(cluster,'.xlsx'))
    
  return(gene_summary_table)

})


