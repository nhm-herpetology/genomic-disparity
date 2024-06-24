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

First, we need to obtain the data for the chromosome 5 FASTA data which we will download from this GitHub repository and NCBI:

```
./chromosome_download.sh
```

</details>

## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>
  
>Landmarks can be any single-copy, conserved sequence that can be aligned to chromosomes in your dataset, but we have used ultraconserved elements (UCEs) in this tutorial as an example. More information about dipteran UCEs can be found [here](https://www.ultraconserved.org/)

We also need the Diptera 2K probe set.

```
wget 
```

</details>
