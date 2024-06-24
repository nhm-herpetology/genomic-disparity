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

We need the Diptera 2K probe set.

```
wget https://raw.githubusercontent.com/nhm-herpetology/genomic-disparity/main/Drosophila-uces/Tetrapods-UCE-5Kv1.fasta
```

</details>

## Step 3: Preparing landmarked chromosomes for downstream analysis

<details>
  <summary>Click to expand content!</summary>

We need to do a few prep steps on the mapped chromsome 5 data before they are ready for Genomic Disparity Analysis...

</details>

## Step 4: Principal components analysis of landmark disparity

<details>
  <summary>Click to expand content!</summary>

Do we see evidence of inversions?

</details>
