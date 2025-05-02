# Structural Disparity Analysis: Sex-linked genes tutorial

![Amniote_header](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Amniote-busco/Amniote_header.jpg)

>This tutorial was developed as part of the ULTRAMOD project https://github.com/AshwiniVM/ULTRAMOD

**Dependencies**

* GFF3 file
* BWA
* samtools
* R statistical software
  

## Step 1: Download generic feature format (gff3) annotation file
<details>
  <summary>Click to expand content!</summary>

>In this tutorial we will use 16 species of amniote vertebrates from Lovell et al. [2022](https://elifesciences.org/articles/78526) and BUSCO landmarks assocaited with the human X chromosome. The groups include representatives from placental mammals, marsupial mammals, monotreme mammals, birds, and squamate reptiles.  

by Poikela et al. [2024](https://academic.oup.com/gbe/article/16/3/evae024/7628473)

Species | Name (Lovell et al.) | Group  | Chromosome
------------ | -------------  | ------------- | -------------
_Mus musculus_	| mouse | Placental | X 
_Choloepus hoffmannis_	| sloth | Placental | X 
_Homo sapiens_	| human | Placental | X 
_Tursiops truncatus_	| dolphin | Placental | X
_Rhinolophus ferrumequinum_	| horseshoeBat | Placental | X
_Sarcophilus harrisii_	| tasmaniandevil | Marsupial| X 
_Trichosurus vulpecula_	| brushtailPossum | Marsupial | X 
_Monodelphis domestica_	| opossum | Marsupial | X 
_Tachyglossus aculeatus_ | echidna | Monotreme | 6  
_Ornithorhynchus anatinus_	| platypus | Monotreme | 6
_Taeniopygia guttata_ | zebrafinch | Bird | 4A
_Cygnus olor_	| swan | Bird | 13
_Calypte anna_ | hummingbird | Bird | 4
_Gallus gallus_ | chicken | Bird | 4
_Lacerta agilis_ | sandLizard | Squamate | Z 
_Thamnophis elegans_ | garterSnake | Squamate | 12

**We will download the gff3 file from Lovell et al. (2022)** 

```  
wget https://github.com/jtlovell/GENESPACE_data/raw/refs/heads/master/vertebrates/gffWithOgs.txt.gz
```

Next, we will expand the file so that we can work with it in R

```  
gunzip gffWithOgs.txt.gz
```

There should now be a file called 'gffWithOgs.txt' in your directory. If not already in your R working directory, please move this file to the working directory before moving to Step 2.

</details>

## Step 2: Select landmark genes of interest
<details>
  <summary>Click to expand content!</summary>

In this example we have selected human X-linked BUSCO genes for use as landmarks. To identify these in the gff3 file we will need to use several R packages.

```  
library(tidyr)
library(dplyr)
library(vegan)
library(stringr)
```

Now let's load the gff3 file into R. We then need to make a function that will ensure gene_id name formats will not create errors during landmark identification. 

```  
original_data <- read.csv("gffWithOgs.csv", header = TRUE)

normalize_gene_id <- function(gene_id) {
  gene_id <- tolower(gene_id)      # Convert to lowercase
  gene_id <- str_replace_all(gene_id, "[-_]", "")  # Remove hyphens and underscores
  return(gene_id)
}

```

Next, we will extract and normalize all BUSCO genes that are located on the human X chromosome

```
human_X_genes <- original_data %>%
  filter(genome == "human" & chr == "X") %>%
  mutate(id = normalize_gene_id(id)) %>%  # Apply normalization function
  select(id)
count(human_X_genes)

```

The count read command should report there are 1838 BUSCO genes on the human X chromosome. Now we will move to Step 3 where find occurences of these genes in other speceis and then create a curated list of landmarks found in a single genomic region across all species. 

  </details>

  ## Step 3: Identify orthologs on chromosomes
<details>
  <summary>Click to expand content!</summary>

Now that we have the genes of interest identified, we can locate them on other genome asseblies in the gff3 annotation file. First, we need to normalize the gene IDs and then find occurrences of them in non-human species included in the file.

```  
original_data <- original_data %>%
  mutate(id = normalize_gene_id(id))

gene_occurrences <- original_data %>%
  filter(id %in% human_X_genes$id)

```

Next, we will count the number of unique chromosomes in each species that contain human X-linked genes.

```

chr_counts <- gene_occurrences %>%
  group_by(genome) %>%
  summarise(matching_chromosomes = n_distinct(chr))  # Count distinct chromosomes per species

chr_counts

```

The final command 'chr_counts' should produce this table: 

genome | matching_chromosomes
------------ | ------------- 
brushtailPossum	| 16 
chicken	| 23 
dolphin | 4
echidna | 4
garterSnake | 9
hoseshowBat | 3
human | 1
hummingbird | 6
mouse | 13
opossum | 15
platypus | 4
sandLizard | 8
sloth | 14
swan | 7
tasmaniandevil | 6
zebrafinch | 7

Instead of listing them in a table, we can confirm that all 16 species have matching chromosomes using the following commands. 

```  
write.csv(chr_counts, "matching_chromosomes_tohumanX.csv", row.names = FALSE)

cat("Number of species with matching chromosomes:", nrow(chr_counts), "\n")

print(chr_counts) 

```

This should result in R telling you 'Number of species with matching chromosomes: 16'. The next series of commands will identify a set of 16 chromosomes (one for each species) with common human X-linked BUSCO genes on them. First, we count the number of genes appearing on each species chromosomes.

```

total_genes_in_human_X <- nrow(human_X_genes)

gene_counts_per_chr <- gene_occurrences %>%
  group_by(genome, chr) %>%
  summarise(matching_gene_count = n(), .groups = "drop")  # Count genes per species' chromosome

```

Next, we will filter chromosomes from all of the species so that we only keep chromosomes with at least 1/8 of the human X-linked genes. This is an arbitraty threshold that you may want to experiment with when using this appraoch with other datasets.  

```

gene_proportions <- gene_counts_per_chr %>%
  mutate(proportion = matching_gene_count / total_genes_in_human_X) %>%
  filter(proportion >= 1/8)

print(gene_proportions, n = 23)

```

The last command here should produce the following table. 

genome | matching_chromosomes | matching_gene_count | proportion
------------ | ------------- | ------------- | ------------- 
brushtailPossum | 2           |        137  |    0.163
brushtailPossum | X             |        193  |   0.230
chicken     |    1             |        136  |    0.162
chicken     |    4             |        249  |    0.297
dolphin     |    X             |        568  |    0.678
echidna     |    6             |        352  |    0.420
garterSnake |    12            |       212  |    0.253
horseshoeBat |    1             |        591  |    0.705
human       |    X             |        838  |    1    
hummingbird |    1             |        127  |    0.152
hummingbird |    4             |        236  |    0.282
mouse       |    X             |        588  |    0.702
opossum     |   X             |        291  |    0.347
platypus    |    6             |        355  |    0.424
sandLizard  |    4             |        122  |    0.146
sandLizard  |    Z             |        221  |    0.264
sloth       |    X             |        419  |    0.5  
swan        |    1             |        130  |    0.155
swan         |   13            |        251  |    0.300
tasmaniandevil | 3             |        147  |    0.175
tasmaniandevil | X             |        332  |    0.396
zebrafinch     | 1             |        128  |    0.153
zebrafinch     | 4A            |        243  |    0.290

This indicates which chromosomes share the highest proportion of X-linked BUSCO landmarks with the human X chromosome. We can see taht most distantly related species have two chromosomes with matching genes whereas closely realted speceis (e.g. other placental mammals) have a single matched chromosome (also the X chromosome). Note: horseshoeBat chromosome 1 is equivalent to the X chromosome in this species. 

Now we will save these results and extract the species and chromosomes that pass filtering

```
write.csv(gene_proportions, "filtered_proportion_matching_humanX.csv", row.names = FALSE)

filtered_chromosomes <- gene_proportions %>%
  select(genome, chr) 

filtered_gene_data <- original_data %>%
  inner_join(filtered_chromosomes, by = c("genome", "chr"))

gene_counts <- filtered_gene_data %>%
  group_by(id) %>%
  summarise(species_count = n_distinct(genome), .groups = "drop")

total_species <- n_distinct(filtered_gene_data$genome)

conserved_genes <- gene_counts %>%
  filter(species_count == total_species) %>%
  select(id)

nrow(conserved_genes)

```

The final command should report that there are 53 genes across a single set of chromosomes for the 16 species. This represents a maximized number of landmarks on a single chromosome. Now we merge the gene IDs back with the position data and export a positions spreadsheet for further preparation. 

```
conserved_gene_data <- filtered_gene_data %>%
  filter(id %in% conserved_genes$id)

write.csv(conserved_gene_data, "conserved_landmarks_humanX.csv", row.names = FALSE)

```

A final checking step was developed to confirm that one chromosome per species is included in your exported spreadsheet. 

```
conserved_gene_data <- read.csv("conserved_landmarks_humanX.csv")

chromosome_count_per_species <- conserved_gene_data %>%
  group_by(genome) %>%
  summarise(unique_chromosomes = n_distinct(chr), .groups = "drop")

print(chromosome_count_per_species)

```

The last command should result in the following table which confirms that there is only one chromosome per species in the filtered dataset.

genome | unique_chromosomes
------------ | ------------- 
brushtailPossum    |              1
chicken            |              1
dolphin            |              1
echidna            |              1
garterSnake        |              1
horseshoeBat       |              1
human              |              1
hummingbird        |              1
mouse              |              1
opossum            |             1
platypus           |             1
sandLizard         |              1
sloth              |              1
swan               |              1
tasmaniandevil     |              1
zebrafinch         |              1

Now we are ready to proceed to Step 4 where we will prepare our spreadsheet for structural disparity analysis. 

  </details>
  
  ## Step 4: Preparing landmarked chromosomes for downstream analysis
<details>
  <summary>Click to expand content!</summary>

To prepare our filtered list of chromosomes for downstream analysis we first need to combine species names with chromosome identities and pivot the table. 


```

conserved_gene_data <- conserved_gene_data %>%
  mutate(genome_chr = paste(genome, chr, sep = "_"))

gene_position_data <- conserved_gene_data %>%
  select(genome_chr, id, start)

gene_position_matrix <- gene_position_data %>%
  pivot_wider(names_from = id, values_from = start)

write.csv(gene_position_matrix, "gene_position_matrix.csv", row.names = FALSE)

```

After this step you will have a CSV file called ```gene_position_matrix.csv```. To improve our estimates of structural disparity, we will need to correct for any orientation issues and bound the landmarks. The rationale for these correctiosn is discussed in Mohan et al. (Under review).  

**Orientation issues**

One way to see evidence of orientation issues is to visualize the placement of the landmarks as analyzed in Lovell et al. (2022) on the chromosomes. In this example of the comparisons between the human X chromosome and the bat X/1 chromosome we can see evidence of a large inversion of landmarks. In the image below BUSCO landmarks are indicated in brown on the grey chromosomes, homologous BUSCO landmarks are connected with gold lines. 

![human_bat_orig](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Amniote-busco/human_bat_original.jpg)

This may be indicative of a real inversion that took place, however, if we take the reverse complement of all the landmark positions in the human chromosome we see that most landmark positions are now matched (syntenically) with the bat chromosomes. This implies that the direction of chromosome sequencing likely explains the putative inversion observed in the first plot (and not a large chromosomal inversion that occurred since humans and bats diverged from a common ancestor). Correcting these type of orientation issues improves the accuracy of structural disparity analysis.  

![human_bat_flipped](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Amniote-busco/human_bat_flipped.jpg)

When chromosomes are identified with likely orientation issues, we can use R to correct (flip) the landmark positions based on the size of each chromosome. Below are the various chromosome sizes from the Lovell et al. (2022) dataset. 

genome | chromosome_size
------------ | ------------- 
brushtailPossum    |              60706338
chicken            |              90861225
dolphin            |             128693445 
echidna            |              61686051
garterSnake        |              60626299
horseshoeBat       |              124933378
human              |              156040895
hummingbird        |              18597117
mouse              |              169476592
opossum            |             89414197
platypus           |             51493492
sandLizard         |              47440541
sloth              |              193839925
swan               |              21471338
tasmaniandevil     |              83081154
zebrafinch         |              19491698

Based on examining the maximized synteny plot from Lovell et al. (2022), we identified seven species that needed to be flipped: (1) human, (2) sloth, (3) brushtail possum, (4) echidna, (5) swan, (6) hummingbird, and (7) zebrafinch.  

```
landmarkflip <- read.csv("gene_position_matrix.csv", header =T, row.names = 1)

landmarkflip <- t(landmarkflip)

landmarkflip <-as.data.frame(landmarkflip)

fun1 <- function(x) {156040895-x+1}
fun2 <- function(x) {193839925-x+1}
fun3 <- function(x) {60706338-x+1}
fun4 <- function(x) {61686051-x+1}
fun5 <- function(x) {21471338-x+1}
fun6 <- function(x) {18597117-x+1}
fun7 <- function(x) {19491698-x+1}

Hum <-lapply(landmarkflip$human_X, fun1)
Slo <-lapply(landmarkflip$sloth_X, fun2)
Bru <-lapply(landmarkflip$brushtailPossum_X, fun3)
Ech <-lapply(landmarkflip$echidna_6, fun4)
Swa <-lapply(landmarkflip$swan_13, fun5)
Hub <-lapply(landmarkflip$hummingbird_4, fun6)
Zeb <-lapply(landmarkflip$zebrafinch_4A, fun7)

landmarkflip$human_X <- Hum
landmarkflip$sloth_X <- Slo
landmarkflip$brushtailPossum_X <- Bru
landmarkflip$echidna_6 <- Ech
landmarkflip$swan_13 <- Swa
landmarkflip$hummingbird_4 <- Hub
landmarkflip$zebrafinch_4A <- Zeb


A <-as.numeric(landmarkflip$human_X)
B <-as.numeric(landmarkflip$sloth_X)
C <-as.numeric(landmarkflip$brushtailPossum_X)
D <-as.numeric(landmarkflip$echidna_6)
E <-as.numeric(landmarkflip$swan_13)
F <-as.numeric(landmarkflip$hummingbird_4)
G <-as.numeric(landmarkflip$zebrafinch_4A)

landmarkflip$human_X <- A
landmarkflip$sloth_X <- B
landmarkflip$brushtailPossum_X <- C
landmarkflip$echidna_6 <- D
landmarkflip$swan_13 <- E
landmarkflip$hummingbird_4 <- F
landmarkflip$zebrafinch_4A <- G

landmarkflip <- t(landmarkflip)

write.csv(landmarkflip, file = "gene_position_matrix_flipped.csv")


```

Now that we have the chromosome synteny maximized we can proceed to bounding procedure which will make disparity estimates less influenced by clustered landmarks. 

**Bounding the landmarks**


  </details>

  ## Step 5: Principal component analysis of landmark disparity
<details>
  <summary>Click to expand content!</summary>

Now that the landmarks have been flipped (when necessary) and bounded, we can conduct the disparity analysis. 

```
gene_position_matrix <- read.csv("gene_position_matrix_flipped_bounded.csv")
head(gene_position_matrix)
# Prepare the data for PCA
# Extract the species information (genome_chr) and gene positions (gene IDs)
pca_data <- gene_position_matrix
species_info <- pca_data$genome_chr  # Save the genome_chr column as the row names

# Remove the genome_chr column, as we won't use it in PCA
pca_data <- pca_data %>% select(-genome_chr)

# Perform PCA (without scaling)
pca_result <- prcomp(pca_data, center = TRUE, scale. = FALSE)
head(pca_result$x)

# Prepare PCA results for plotting
pca_df <- data.frame(pca_result$x)  
head(pca_df)
write.csv(pca_df, file="PCA_humanX.csv")
pca_df$genome_chr <- species_info  # Add the genome_chr column back to the PCA result for coloring
head(pca_df$genome_chr)
head(pca_df)
write.csv(pca_df, file="PCA_humanX_bounded.csv")

```


  </details>
