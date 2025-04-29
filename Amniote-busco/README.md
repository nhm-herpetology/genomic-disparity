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

Species | Group  | Chromosome
------------ | -------------  | -------------
_Mus musculus_	| Placental | X 
_Choloepus hoffmannis_	| Placental | X 
_Homo sapiens_	| Placental | X 
_Tursiops truncatus_	| Placental | X
_Rhinolophus ferrumequinum_	| Placental | X
_Sarcophilus harrisii_	| Marsupial| X 
_Trichosurus vulpecula_	| Marsupial | X 
_Monodelphis domestica_	| Marsupial | X 
_Tachyglossus aculeatus_ | Monotreme | 6  
_Ornithorhynchus anatinus_	| Monotreme | 6
_Taeniopygia guttata_ | Bird | 4A
_Cygnus olor_	| Bird | 13
_Calypte anna_ | Bird | 4
_Gallus gallus_ | Bird | 4
_Lacerta agilis_ | Squamate | Z 
_Thamnophis elegans_ | Squamate | 12

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
brushtailPossum | 2   |          |        137  |    0.163
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


# Save the results

  </details>
  
  ## Step 4: Preparing landmarked chromosomes for downstream analysis
<details>
  <summary>Click to expand content!</summary>

  </details>

  ## Step 5: Principal component analysis of landmark disparity
<details>
  <summary>Click to expand content!</summary>

  </details>
