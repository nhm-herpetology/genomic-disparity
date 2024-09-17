##########################################################################
######################     STEPS in R   ################################
###########################################################################

##install package 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BUSpaRse")

### Load packages ################
library("BUSpaRse")

###############################  extract genes annotated on the specific chromosome assembly ##########################################################
transcripts <- tr2g_gff3("cricetulus_griseus_chr6.gff3",
                         #Genome = "blablabla.fa",
                         get_transcriptome = FALSE,
                         out_path = ".",
                         write_tr2g = TRUE,
                         transcript_id = "transcript_id",
                         gene_id = "gene_id",
                         gene_name = "Name",
                         transcript_version = "version",
                         gene_version = "version",
                         version_sep = ".",
                         transcript_biotype_col = "biotype",
                         gene_biotype_col = "biotype",
                         transcript_biotype_use = "all",
                         gene_biotype_use = "all",
                         chrs_only = TRUE,
                         compress_fa = FALSE,
                         save_filtered_gff = TRUE,
                         overwrite = FALSE,
                         source = c("ensembl")

########     the output file name is tr2g and contains all the gene info we need! open this in excel sheet for further analyses   ################################

##################################################################################################################################################################
################################   Gene_ID files have been prepped to remove duplicate gene names and multiple NAs         #######################################
##################################################################################################################################################################

setwd("path/to/folder")

#### Load the required libraries
library(dplyr)

######### Jaccard's similiarity function which is not case sensitive to gene names e.g. Rag1 versus RAG1  #######################
#################################################################################################################################

jaccard_new <- function(a, b) {
  a <- tolower(a)
  b <- tolower(b)
  
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


###### Read your data from the CSV file 
data <- read.csv("chromosome_set_2.csv")


#####  Get the column names
column_names <- colnames(data)
column_names

####### Create an empty data frame to store the results
result_df <- data.frame("Column1" = character(0), "Column2" = character(0), "Jaccard_Similarity" = numeric(0))

######### Loop through all pairs of columns
for (i in 1:(length(column_names) - 1)) {
  for (j in (i + 1):length(column_names)) {
    #### Get the column names for the current pair
    col1_name <- column_names[i]
    col2_name <- column_names[j]
    
    ##### Extract the data for the current pair of columns
    col1_data <- data[[col1_name]]
    col2_data <- data[[col2_name]]
    
    ##### Calculate the Jaccard similarity
    jaccard_sim <- jaccard_new(col1_data, col2_data)
    
    ###### Add the results to the data frame
    result_df <- rbind(result_df, data.frame("Column1" = col1_name, "Column2" = col2_name, "Jaccard_Similarity" = jaccard_sim))
  }
}

####   Write the results to an output file
write.csv(result_df, "jaccard_new_similarity_set2.csv", row.names = FALSE)





