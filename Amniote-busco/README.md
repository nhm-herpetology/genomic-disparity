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

  </details>
  
  ## Step 4: Preparing landmarked chromosomes for downstream analysis
<details>
  <summary>Click to expand content!</summary>

  </details>

  ## Step 5: Principal component analysis of landmark disparity
<details>
  <summary>Click to expand content!</summary>

  </details>
