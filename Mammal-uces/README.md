# Placental Mammal UCE Tutorial

>This protocol was developed as part of the ULTRAMOD project https://github.com/AshwiniVM/ULTRAMOD

**Dependencies**

* NCBI esearch
* BWA
* R statistical software

## Step 1: Download chromosome sequences from genome assemblies
<details>
  <summary>Click to expand content!</summary>

We will need to download chromosome-level genome assemblies from NCBI or other repository. 

</details>


## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>

  Landmarks can be any conserved sequence that can be aligned to genomes in your dataset, but we have used ultraconserved elements (UCEs) as an example. BWA is used to map the landmarks to different chromosomes.  

</details>

## Step 3: Cluster chromosomes based on landmark similarity

<details>
  <summary>Click to expand content!</summary>

  MMDS in R statistical software is used to identify which chromosomes likley contain homologous blocks of genomes (i.e. supergenes, Marian fragments etc.). 

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
