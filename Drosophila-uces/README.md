# Genomic Disparity Analysis: _Drosophila virilis_ group tutorial


>This tutorial was developed as part of the ULTRAMOD project https://github.com/AshwiniVM/ULTRAMOD

Chromosomal inversions are thought to be involved in speciation, the formation of supergenes and other evolutionary phenomena. In this tutorial we will use Genomic Disparity Analysis to identify chromosomal inversions that have occurred in chromosome 5 of the _Drosophila virilis_ group. 

**Dependencies**

* NCBI Entrez Direct UNIX E-utilities
* BWA
* samtools
* R statistical software
  

## Step 1: Download chromosome 5 sequences
<details>
  <summary>Click to expand content!</summary>

>In this tutorial we will use six examples of _Drosophila_ chromosome 5 which has known chromosomal inversions identified by Poikela et al. (2024) in their study of _Drosophila_ speciation.

Species | Autosome | GenBank/Reference 
------------ | -------------  | -------------
_Drosophila americana_	| 5 | CM061086.1 
_Drosophila flavomontana_	| 5 | Poikela et al. (2024)
_Drosophila montana_	|  5 | Poikela et al. (2024)
_Drosophila novamexicana_	| 5 | CM061080.1
_Drosophila virilia_	| 5 | CM017608.2
_Drosophila virilia_	| 5 | CM061075.1

>The two chromosome 5 assemblies from Poikela et al. (2024) were obtained from https://zenodo.org/records/10635471.

To obtain the FASTA data for chromosome 5 we will download files from NCBI and this GitHub repository:

```
esearch -db nucleotide -query CM061086.1 | efetch -format fasta > D_americana_c5.fasta
esearch -db nucleotide -query  CM061080.1 | efetch -format fasta > D_novamexicana_c5.fasta
esearch -db nucleotide -query CM017608.2 | efetch -format fasta > D_virilia_c5_A.fasta.fasta
esearch -db nucleotide -query CM061075.1 | efetch -format fasta > D_virilia_c5_B.fasta.fasta

wget https://raw.githubusercontent.com/nhm-herpetology/genomic-disparity/main/Drosophila-uces/D_montana_c5.zip
wget https://raw.githubusercontent.com/nhm-herpetology/genomic-disparity/main/Drosophila-uces/D_flavomontana_c5.zip

unzip D_montana_c5.zip
unzip D_flavomontana_c5.zip

mv D_montata_c5/D_montata_c5.fasta ./
mv D_flavomontata_c5/D_flavomontata_c5.fasta ./
```
>Note:NCBI Entrez Direct UNIX E-utilities need to be in your $PATH for the above commands to work correctly. 

Now that all chromosome 5 sequences for the six taxa are in your working directory, we can move on to mapping the UCE landmarks on each chromosome. 

</details>

## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>
  
>Landmarks can be any single-copy, conserved sequence that can be aligned to chromosomes in your dataset, but we have used ultraconserved elements (UCEs) in this tutorial as an example. More information about dipteran UCEs can be found [here](https://www.ultraconserved.org/)

We need to download the ```Diptera-UCE-2.7K-v1.fasta``` file to map to the different chromosomes. 

```
wget https://raw.githubusercontent.com/nhm-herpetology/genomic-disparity/main/Drosophila-uces/Diptera-UCE-2.7K-v1.fasta
```

We can use a looping command to conduct BWA mapping on all of the chromosomes in the directory. 

```
for i in $(ls *.fasta)
do
	bwa index $i
	bwa mem $i Diptera-UCE-2.7K-v1.fasta -t 6 > bwa_mem_align_UCEs_$i.sam 
	samtools view -S -b bwa_mem_align_UCEs_$i.sam > UCE_$i.bam 
	samtools sort UCE_$i.bam  -o UCE_$i.sorted.bam 
	samtools index UCE_$i.sorted.bam 
	samtools view -F 4 UCE_$i.bam > mapped_$i.sam 
	wc -l mapped_$i.sam  | tr "," "." > UCEcount_$i.csv
done

rm mapped_Diptera-UCE-2.7K-v1.fasta.sam
rm UCW_Diptera-UCE-2.7K-v1.fasta.bam
rm UCE_Diptera-UCE-2.7K-v1.fasta.sorted.bam
rm UCE_Diptera-UCE-2.7K-v1.fasta.sorted.bam.bai
rm UCEcount_Diptera-UCE-2.7K-v1.fasta.csv

for f in *.sam; do
    mv "$f" "$(basename "$f" .sam).tsv"
done

mkdir chromosome_5
mv mapped_*.tsv chromosome_5/
cd chromosome_5/

cut -f1,3 *.tsv > Landmarks_merged.tsv
sed 's/\t/,/g' Landmarks_merged.tsv > Landmarks_merged.csv
sed ' 1 s/.*/chromosomes,uces/' Landmarks_merged.csv > sample_input_pres_abs.csv
```

You will now have a file called ```sample_input_pres_abs.csv``` that reports which landmarks were mapped to which chromosomes. We will use this information to prepare the chromosomes for Genomic Disparity Analysis. 

</details>

## Step 3: Preparing landmarked chromosomes for downstream analysis

<details>
  <summary>Click to expand content!</summary>

>We need to remove any landmarks that the chromosomes do not share, check the directionality/orientation of the chromosomes, and merge UCE landmarks before Genomic Disparity Analysis.

We will use the ```landmark_pres_abs.csv``` file from **Step 2** to remove UCEs that the chromosomes do not share in R:

```
library(reshape2)
Input_pres_abs <- as.data.frame(read.csv("sample_input_pres_abs.csv", stringsAsFactors = F))
matrix <-dcast(Input_pres_abs, chromosomes ~ uces, length)
write.csv(matrix, file = 'landmark_pres_abs.csv')

data1 <- read.csv("landmark_pres_abs.csv", row.names = 1)
chr <- data1
chr_matchedUCEs <- apply(chr, 1, function(row) all(row != 0))
chr_clean <- chr[chr_matchedUCEs,]
write.csv(chr_clean, file = "present_landmarks.csv")
```
>This procedure should result in the identification of 413 UCE landmarks that the _Drosophila_ have in common on chromosome 5.

Next, we will extract position information for the 413 UCE landmarks from the BWA-mapping using R: 

```
library(matrixStats)
library(dplyr)

homologousUCE <- read.csv("present_landmarks.csv", row.names = 1)

folder_path <- "./chromosome_5"
file_list <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

matrices <- list()
for (file in file_list) {
  species <- gsub(".tsv", "", basename(file))
  data <- read.table(file, header = FALSE, sep = "\t")
  matrix <- as.matrix(data)
  matrices[[species]] <- matrix
}

for (species in names(matrices)) {
  df_name <- species
  df <- as.data.frame(matrices[[species]])
  df <- df[, c(1, 2, 4)]
  df <- data.frame(lapply(df, function(x) {gsub("-", ".", x)}))
  colnames(df)[colnames(df) == "V1"] <- "chromosomes"
  colnames(df)[colnames(df) == "V4"] <- df_name
  homologousUCE <- left_join(homologousUCE, df, by = "chromosomes")
}

columns_to_remove <- colnames(homologousUCE)[apply(homologousUCE == 1, 2, all)]
homologousUCE <- homologousUCE[, !colnames(homologousUCE) %in% columns_to_remove]

write.csv(homologousUCE, file = "homologous_UCEs_extracted.csv", row.names = TRUE)
```

</details>

## Step 4: Principal components analysis of landmark disparity

<details>
  <summary>Click to expand content!</summary>

Do we see evidence of inversions?

</details>
