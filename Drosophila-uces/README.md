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

cat UCEcount_*.csv > Total_UCE_counts.txt
awk '(NR == 1) || (FNR == 1)' bwa_mem_align_UCEs_*.fasta.sam > Chromosomes_lengths.txt

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
  colnames(df)[colnames(df) == "V1"] <- "chromosomes"
  colnames(df)[colnames(df) == "V4"] <- df_name
  homologousUCE <- left_join(homologousUCE, df, by = "chromosomes")
}

columns_to_remove <- colnames(homologousUCE)[apply(homologousUCE == 1, 2, all)]
homologousUCE <- homologousUCE[, !colnames(homologousUCE) %in% columns_to_remove]

write.csv(homologousUCE, file = "homologous_UCEs_extracted.csv", row.names = TRUE)
```

After this step you will have a CSV file called ```homologous_UCEs_extracted.csv``` which contains two columns for each species: (1) the 'V2' direction that landmarks were mapped on the chromosomes (0 [forward] or 16 [reverse]) and (2) the mapping location (in base pairs) of each landmark. Sometimes chromosomes are assembled with opposite complementarities. For example, if we look at the first 9 landmark positions in the ```homologous_UCEs_extracted.csv``` file we should see this:   


V2.x |	mapped_D_americana_c5.fasta | V2.y | mapped_D_flavomontana_c5.fasta | V2.x.x | mapped_D_montana_c5.fasta | V2.y.y | mapped_D_novamexicana_c5.fasta | V2.x.x.x | mapped_D_virilis_c5_A.fasta | V2.y.y.y | mapped_D_virilis_c5_B.fasta
------------ | -------------  | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
16 | 9468162 | 0 | 17552991 | 0 | 15039180 | 16 | 9451004 | 16 | 9109883 | 16 | 9207235
16 | 9468140 | 0 | 17553031 | 0 | 15039220 | 16 | 9450982 | 16 | 9109843 | 16 | 9207195
16 | 22024688 | 0 | 8535617 | 16 | 22114936 | 16 | 21807240 | 16 | 21487388 | 16 | 21571438
16 | 22024688 | 0 | 8535614 | 16 | 22114936 | 16 | 21807240 | 16 | 21487388 | 16 | 21571438
16 | 19346268 | 0 | 11231717 | 16 | 19370570 | 16 | 19130702 | 16 | 18910683 | 16 | 18976268
16 | 19346266 | 0 | 11231665 | 16 | 19370568 | 16 | 19130700 | 16 | 18910681 | 16 | 18976266
0 | 25602768 | 16 | 5075950 | 0 | 25564769 | 0 | 25229190 | 0 | 24875739 | 0 | 24933228
0 | 25602808 | 16 | 5075916 | 0 | 25564809 | 0 |  25229230 | 0 | 24875779 | 0 | 24933268
16 | 25543223 | 0 | 5138695 | 16 | 25503819 | 16 | 25169374 | 16 | 24816581 | 16 | 24874057

We can see that for almsot all of these landmarks (7 out of 9) _Drosophila flavomontana_ has the opposite landmark direction as the other taxa. This pattern continues in the rest of the ```homologous_UCEs_extracted.csv``` file, so we conclude that this species chromosome sequence was uploaded with an opposite complementarity when compared to the other _Drosophila_. Determining which taxa to 'flip' is arbitrary for disparity analysis, but for the tutorial dataset we will 'flip' _Drosophila flavomontana_ so that the landmark positions are more consistent with the other taxa. To 'flip' we need to know the total lengths of each chromosome. We can find the total lengths of of the chromosomes in the ```Chromosome_lengths.tsv``` files we generated during **Step 2**.

Species | GenBank/Reference | Chromosome 5 Length 
------------ | -------------  | -------------
_Drosophila americana_	| CM061086.1 | 27587546  
_Drosophila flavomontana_ | Poikela et al. (2024) | 27217941
_Drosophila montana_	| Poikela et al. (2024) | 26508887
_Drosophila novamexicana_ | CM061080.1 | 26715699
_Drosophila virilia_	| CM017608.2 | 27902728
_Drosophila virilia_	| CM061075.1 | 27785111

The 'flipping' step can be conducted in Excel or similar spreadsheet editor. We can also use R to 'flip' the _Drosophila flavomontana_ chromosome using a function command: 

```
landmarkflip <- read.csv("homologous_UCEs_extracted.csv", row.names = 1)
fun1 <- function(x) {27217941-x}
Drfl <-lapply(landmarkflip$mapped_D_flavomontana_c5.fasta, fun1)
landmarkflip$mapped_D_flavomontana_c5.fasta <- Drfl
A <-as.numeric(landmarkflip$mapped_D_flavomontana_c5.fasta)
landmarkflip$mapped_D_flavomontana_c5.fasta <-A
```
>Note: In chromosomes with relatively conserved landmark placements, it should be obvious which taxa need to be 'flipped'. However, when landmarks are more evolutionarily labile it may be diffcult to justify a 'flipping' operation, so we encourage users to think about this operation carefully. 

Now all of the chromosomes are positioned correctly and we can remove the mapping direction information and export our final file:

```
data2 <- landmarkflip[-c(2,4,6,8,10,12)]

write.csv(data2, file = "homologous_UCEs_extracted_flipped.csv")
```

After the last step we have a file called ```homologous_UCEs_extracted_flipped.csv``` that is ready for the final preparation steps. This includes accounting for a specific caveat of using the UCE probe set. The UCE probe set was developed to capture UCEs across diverse taxa, as such some UCEs are targeted by multiple probes, so to control for the variation this creates in mapping, we average the probe placement across landmarks targeting the same UCE. We will remove the probe numbers (e.g. p1) in order to merge information from the probes targeting multiple parts of the same UCE landmark. After we average the UCE positions, we transpose the matrix to prepare it for the PCA: 

```
data3 <- read.csv("homologous_UCEs_extracted_flipped.csv")
data3$chromosomes <- sub("_.*", "", data3$chromosomes)
data3$X <-NULL

data4 <-aggregate(data3, by = list(data3$chromosomes), mean)
data4$chromosomes<-NULL
colnames(data4)[colnames(data4) == "Group.1"] <- "chromosomes"

transposed_data <- data4 %>% t() %>% as.data.frame()

write.csv(transposed_data, file = "homologous_UCEs_chromosome_5_PCA.csv")
```
  
After averaging the positions, there should be **191 UCE landmarks**. This final prep will produce the file ```homologous_UCEs_chromosome_5_PCA.csv``` which is the input file for Genomic Disparity Analysis as outlined in the final step of the tutorial. 

</details>

## Step 4: Principal components analysis of landmark disparity

<details>
  <summary>Click to expand content!</summary>

>In this final step we will visualize the disparity in landmark placement using Principal Components Analysis (PCA). 

At the end of **Step 3** we will have the file ```homologous_UCEs_chromosome_5_PCA.csv```, which we will now load into R for the PCA:  

```
library(dplyr)
library(ggplot2)
library(cowplot)

C5_PCA <- read.csv("homologous_UCEs_chromosome_5_PCA.csv", row.names = 1)
names(C5_PCA) <- C5_PCA[1,]
C5_PCA <- C5_PCA[-1,]

PCprep <- C5_PCA %>% mutate_at(1:191, as.numeric)

PCAC5 <-prcomp(PCprep)
```

After the PCA has completed there will be five lists of results including "sdev", "rotation", "center", "scale", and "x". The "x" variable contains the PC scores that we will use for Genomic Disparity analysis. We will now extract them and visualize the genomic disparity results in ggplot2.

```
PCAC5_scores <-as.data.frame(PCAC5$x)

Species <-c("americana", "flavomontana", "montana", "novamexicana", "virilis", "virilis")

PCAC5_plots <-cbind(PCAC5_scores, Species)

P1 <-ggplot(PCAC5_plots, aes(x = PC1, y = PC2, color = Species)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("americana", "flavomontana", "montana", "novamexicana", "virilis"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

P2 <-ggplot(PCAC5_plots, aes(x = PC3, y = PC4, color = Species)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("americana", "flavomontana", "montana", "novamexicana", "virilis"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

plot_grid(P1, P2, ncol = 1)
```

</details>
