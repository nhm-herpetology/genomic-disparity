# Genomic Disparity: Placental Mammal UCE Tutorial

![Mammal-uces-PC2-PC3](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/Mammal-uces-PC2-PC3.jpg)

>This tutorial was developed as part of the ULTRAMOD project https://github.com/AshwiniVM/ULTRAMOD

**Dependencies**

* NCBI Entrez Direct UNIX E-utilities
* BWA
* samtools
* R statistical software
  

## Step 1: Download chromosome sequences from genome assemblies
<details>
  <summary>Click to expand content!</summary>

>In this tutorial we will use 26 placental mammal species belonging to five different orders. We will focus on the available autosomes from each assembly (ignoring sex chromosomes and unplaced scaffolds). 

Species | Order  | Autosomes | GenBank 
------------ | -------------  | ------------- | ------------- 
_Bos indicus_	| Artiodactyla | 29 | CM003021.1-CM003049.1
_Bos taurus_	| Artiodactyla | 29 | CM008168.2-CM008196.2
_Bubalus bubalis_	| Artiodactyla | 24 | CM034271.1-CM034294.1
_Capra aegagrus_	| Artiodactyla | 29 | CM003214.1-CM003243.1
_Capra hircus_	| Artiodactyla | 29 | CM004562.1-	CM004590.1
_Giraffa tippelskirchi_	| Artiodactyla | 14 | CM018103.1-CM018116.1 
_Ovis aries_	| Artiodactyla | 26 | CM028704.1-CM028729.1 
_Ceratotherium simum_	| Perissodactyla | 40 | CM043809.1-CM043848.1
_Equus asinus_ | Perissodactyla | 30 | CM027690.2-CM027719.2 
_Equus caballus_	| Perissodactyla | 31 | CM009148.1-CM009178.1
_Felis catus_ | Carnivora | 18 | CM031412.1-CM031429.1
_Neomonachus schauinslandi_	| Carnivora | 16 | CM035896.1-CM035912.1 
_Panthera tigris_ | Carnivora | 18 | 	CM031431.1-CM031448.1 
_Cricetulus griseus_ | Rodentia | 9 | CM023436.1-CM023444.1 
_Mus caroli_ | Rodentia | 19 | LT608242.1-LT608248.1 
_Mus musculus_ | Rodentia | 19 | CM000994.3-CM001012.3
_Mus pahari_ | Rodentia | 23 | LT608296.1-LT608309.1
_Mus spretus_	| Rodentia | 19 | OW971678.1-OW971697.1
_Rattus norvegicus_	| Rodentia | 20 | CM070391.1-	CM070410.1 
_Peromyscus maniculatus_ | Rodentia | 23 | CM010879.2-CM010901.2
_Gorilla gorilla_	| Primates | 23 | CM055446.2-CM068951.1
_Macaca fascicularis_ | Primates | 20 | BLPH02000001.1-BLPH02000020.1 
_Macaca mulatta_	| Primates | 20 | CM014336.1-CM014355.1
_Pan troglodytes_ | Primates | 23 | CM054434.2-	CM068906.1 
_Papio anubis_ | Primates | 19 | CM018180.1-CM018198.1
_Piliocolobus tephrosceles_ | Primates | 21 | 	CM019240.1-CM019260.1

**We will download the chromosome data from genome assemblies using Entrez Direct UNIX E-utilities** 

First, we need to install the E-utilities:

```  
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

Next, we will download the chromosome sequences using the configuration file ```chromosome_download.config``` and the shell script  ```chromosome_download.sh```. The configuration file is formatted with one species or individual per row, with the first two columns being sample info, and subsequent columns being accessions numbers of chromosome assemblies to use with the pipeline. For example: 

```
Taxon_1 "CM017296.1" "CM017297.1" "CM017298.1" "CM017299.1" "CM017300.1" "CM017301.1" "CM017302.1" "CM017303.1" "CM017304.1" "CM017305.1" "CM017306.1" "CM017307.1" "CM017308.1" "CM017309.1" "CM017310.1" "CM017311.1" "CM017312.1" "CM017313.1" "CM017314.1" "CM017315.1" "CM017316.1" "CM017317.1" "CM017318.1" "CM017319.1" "CM017320.1" 
Taxon_2 "CM017296.1" "CM017297.1" "CM017298.1" "CM017299.1" "CM017300.1" "CM017301.1" "CM017302.1" "CM017303.1" "CM017304.1" "CM017305.1" "CM017306.1" "CM017307.1" "CM017308.1" "CM017309.1" "CM017310.1" "CM017311.1" "CM017312.1" "CM017313.1" "CM017314.1"  
Taxon_3 "CM017296.1" "CM017297.1" "CM017298.1" "CM017299.1" "CM017300.1" "CM017301.1" "CM017302.1" "CM017303.1" "CM017304.1" "CM017305.1" "CM017306.1" "CM017307.1" "CM017308.1" "CM017309.1" "CM017310.1" "CM017311.1"  
Taxon_4 "CM017296.1" "CM017297.1" "CM017298.1" "CM017299.1" "CM017300.1" "CM017301.1" "CM017302.1" "CM017303.1" "CM017304.1" "CM017305.1" "CM017306.1" "CM017307.1" "CM017308.1" "CM017309.1" 
```

After the configuration file is ready we make the downloading shell script executable and then run it. Note: for the script to work the Entrez Direct UNIX E-utilities needs to be in your $PATH.

```
chmod +x chromosome_download.sh
```

```
./chromosome_download.sh
```

Depending on the number of taxa you are using, this may download a substantial amount of data. It may take up a while to complete but progress updates will be sent from the script as each taxon is processed for you to track progress. The tutorial mammal dataset (26 species) takes XX hours to download using 8 CPUS and 40G of RAM.   

</details>


## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>

  **Landmarks can be any single-copy, conserved sequence that can be aligned to genomes in your dataset, but we have used ultraconserved elements (UCEs) as an example.**  

Download the Tetrapod 5k probe sequences (these will be used to UCE locations on the chromosomes). The probe set can also be downloaded [here](https://www.ultraconserved.org/)
  
```
wget https://raw.githubusercontent.com/nhm-herpetology/museum-NGS-training/main/Unit_03/Bioinformatics_Lab/Tetrapods-UCE-5Kv1.fasta
``` 
</details>

## Step 3: Cluster chromosomes based on landmark similarity

<details>
  <summary>Click to expand content!</summary>

  MMDS in R statistical software is used to identify which chromosomes likley contain homologous blocks of genomes (i.e. supergenes, Marian fragments etc.). 

Following this procedure you should have identified the following chromosomes as belonging to **Chromosome Set 1**

Species | Order  | Chromosome (in MMDS set 1) 
------------ | -------------  | ------------- 
_Bos indicus_	| Artiodactyla | 2
_Bos taurus_	| Artiodactyla | 2
_Bubalus bubalis_	| Artiodactyla | 2
_Capra aegagrus_	| Artiodactyla | 2
_Capra hircus_	| Artiodactyla | 2
_Giraffa tippelskirchi_	| Artiodactyla | 3
_Ovis aries_	| Artiodactyla | 2
_Ceratotherium simum_	| Perissodactyla | 9
_Equus asinus_ | Perissodactyla | 4
_Equus caballus_	| Perissodactyla | 18
_Felis catus_ | Carnivora | C1
_Neomonachus schauinslandi_	| Carnivora | 3
_Panthera tigris_ | Carnivora | C1
_Cricetulus griseus_ | Rodentia | 6
_Mus caroli_ | Rodentia | 2
_Mus musculus_ | Rodentia | 2
_Mus pahari_ | Rodentia | 3
_Mus spretus_	| Rodentia | 2
_Rattus norvegicus_	| Rodentia | 3
_Peromyscus maniculatus_ | Rodentia | 4
_Gorilla gorilla_	| Primates | 3
_Macaca fascicularis_ | Primates | 12
_Macaca mulatta_	| Primates | 12
_Pan troglodytes_ | Primates | 2B
_Papio anubis_ | Primates | 10
_Piliocolobus tephrosceles_ | Primates | 11

</details>

## Step 4: Remove non-homologous landmarks from chromosomes

<details>
  <summary>Click to expand content!</summary>

  R statistical software is used to remove the landmarks that do not match in the refined chromosome set.  

</details>

## Step 5: Principal components analysis of landmark disparity

<details>
  <summary>Click to expand content!</summary>

  R statistical software is used to make exciting plots!

</details>
