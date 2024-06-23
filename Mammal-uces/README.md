# Genomic Disparity Analysis: Placental Mammal UCE Tutorial

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
Taxon_1 Bos_indicus "CM003021.1" "CM003022.1" "CM003023.1" "CM003024.1" "CM003025.1" "CM003026.1" "CM003027.1" "CM003028.1" "CM003029.1" "CM003030.1" "CM003031.1" "CM003032.1" "CM003033.1" "CM003034.1" "CM003035.1" "CM003036.1" "CM003037.1" "CM003038.1" "CM003039.1" "CM003040.1" "CM003041.1" "CM003042.1" "CM003043.1" "CM003044.1" "CM003045.1" "CM003046.1" "CM003047.1" "CM003048.1" "CM003049.1"
Taxon_2 Bos_taurus "CM008168.2" "CM008169.2" "CM008170.2" "CM008171.2" "CM008172.2" "CM008173.2" "CM008174.2" "CM008175.2" "CM008176.2" "CM008177.2" "CM008178.2" "CM008179.2" "CM008180.2" "CM008181.2" "CM008182.2" "CM008183.2" "CM008184.2" "CM008185.2" "CM008186.2" "CM008187.2" "CM008188.2" "CM008189.2" "CM008190.2" "CM008191.2" "CM008192.2" "CM008193.2" "CM008194.2" "CM008195.2" "CM008196.2"
Taxon_3 Bubalus_bubalis "CM034271.1" "CM034272.1" "CM034273.1" "CM034274.1" "CM034275.1" "CM034276.1" "CM034277.1" "CM034278.1" "CM034279.1" "CM034280.1" "CM034281.1" "CM034282.1" "CM034283.1" "CM034284.1" "CM034285.1" "CM034286.1" "CM034287.1" "CM034288.1" "CM034289.1" "CM034290.1" "CM034291.1" "CM034292.1" "CM034293.1" "CM034294.1"
Taxon_4 Capra_aegagrus "CM003214.1" "CM003215.1" "CM003216.1" "CM003217.1" "CM003218.1" "CM003219.1" "CM003220.1" "CM003221.1" "CM003222.1" "CM003223.1" "CM003225.1" "CM003226.1" "CM003227.1" "CM003228.1" "CM003229.1" "CM003230.1" "CM003231.1" "CM003232.1" "CM003233.1" "CM003234.1" "CM003235.1" "CM003236.1" "CM003237.1" "CM003238.1" "CM003239.1" "CM003240.1" "CM003241.1" "CM003242.1" "CM003243.1"
Taxon_5 Capra_hircus "CM004562.1" "CM004563.1" "CM004564.1" "CM004565.1" "CM004566.1" "CM004567.1" "CM004568.1" "CM004569.1" "CM004570.1" "CM004571.1" "CM004572.1" "CM004573.1" "CM004574.1" "CM004575.1" "CM004576.1" "CM004577.1" "CM004578.1" "CM004579.1" "CM004580.1" "CM004581.1" "CM004582.1" "CM004583.1" "CM004584.1" "CM004585.1" "CM004586.1" "CM004587.1" "CM004588.1" "CM004589.1" "CM004590.1"
Taxon_6 Giraffa_tippelskirchi "CM018103.1" "CM018104.1" "CM018105.1" "CM018106.1" "CM018107.1" "CM018108.1" "CM018109.1" "CM018110.1" "CM018111.1" "CM018112.1" "CM018113.1" "CM018114.1" "CM018115.1" "CM018116.1"
Taxon_7 Ovis_aries "CM028704.1" "CM028705.1" "CM028706.1" "CM028707.1" "CM028708.1" "CM028709.1" "CM028710.1" "CM028711.1" "CM028712.1" "CM028713.1" "CM028714.1" "CM028715.1" "CM028716.1" "CM028717.1" "CM028718.1" "CM028719.1" "CM028720.1" "CM028721.1" "CM028722.1" "CM028723.1" "CM028724.1" "CM028725.1" "CM028726.1" "CM028727.1" "CM028728.1" "CM028729.1"
Taxon_8 Ceratotherium_simum "CM043809.1" "CM043810.1" "CM043811.1" "CM043812.1" "CM043813.1" "CM043814.1" "CM043815.1" "CM043816.1" "CM043817.1" "CM043818.1" "CM043819.1" "CM043820.1" "CM043821.1" "CM043822.1" "CM043823.1" "CM043824.1" "CM043825.1" "CM043826.1" "CM043827.1" "CM043828.1" "CM043829.1" "CM043830.1" "CM043831.1" "CM043832.1" "CM043833.1" "CM043834.1" "CM043835.1" "CM043836.1" "CM043837.1" "CM043838.1" "CM043839.1" "CM043840.1" "CM043841.1" "CM043842.1" "CM043843.1" "CM043844.1" "CM043845.1" "CM043846.1" "CM043847.1" "CM043848.1"
Taxon_9 Equus_asinus "CM027690.2" "CM027691.2" "CM027692.2" "CM027693.2" "CM027694.2" "CM027695.2" "CM027696.2" "CM027697.2" "CM027698.2" "CM027699.2" "CM027700.2" "CM027701.2" "CM027702.2" "CM027703.2" "CM027704.2" "CM027705.2" "CM027706.2" "CM027707.2" "CM027708.2" "CM027709.2" "CM027710.2" "CM027711.2" "CM027712.2" "CM027713.2" "CM027714.2" "CM027715.2" "CM027716.2" "CM027717.2" "CM027718.2" "CM027719.2"
Taxon_10 Equus_caballus "CM027690.2" "CM027691.2" "CM027692.2" "CM027693.2" "CM027694.2" "CM027695.2" "CM027696.2" "CM027697.2" "CM027698.2" "CM027699.2" "CM027700.2" "CM027701.2" "CM027702.2" "CM027703.2" "CM027704.2" "CM027705.2" "CM027706.2" "CM027707.2" "CM027708.2" "CM027709.2" "CM027710.2" "CM027711.2" "CM027712.2" "CM027713.2" "CM027714.2" "CM027715.2" "CM027716.2" "CM027717.2" "CM027718.2" "CM027719.2"
Taxon_11 Felis_catus "CM031412.1" "CM031413.1" "CM031414.1" "CM031415.1" "CM031416.1" "CM031417.1" "CM031418.1" "CM031419.1" "CM031420.1" "CM031421.1" "CM031422.1" "CM031423.1" "CM031424.1" "CM031425.1" "CM031426.1" "CM031427.1" "CM031428.1" "CM031429.1"
Taxon_12 Neomonachus_schauinslandi "CM035896.1" "CM035898.1" "CM035899.1" "CM035900.1" "CM035901.1" "CM035902.1" "CM035903.1" "CM035904.1" "CM035905.1" "CM035906.1" "CM035907.1" "CM035908.1" "CM035909.1" "CM035910.1" "CM035911.1" "CM035912.1"
Taxon_13 Panthera_tigris "CM031431.1" "CM031432.1" "CM031433.1" "CM031434.1" "CM031435.1" "CM031436.1" "CM031437.1" "CM031438.1" "CM031439.1" "CM031440.1" "CM031441.1" "CM031442.1" "CM031443.1" "CM031444.1" "CM031445.1" "CM031446.1" "CM031447.1" "CM031448.1"
Taxon_14 Cricetulus_griseus "CM023436.1" "CM023437.1" "CM023438.1" "CM023439.1" "CM023440.1" "CM023441.1" "CM023442.1" "CM023443.1" "CM023444.1"
Taxon_15 Mus_caroli "LT608242.1" "LT608244.1" "LT608232.1" "LT608246.1" "LT608240.1" "LT608245.1" "LT608243.1" "LT608237.1" "LT608231.1" "LT608233.1" "LT608241.1" "LT608234.1" "LT608247.1" "LT608238.1" "LT608239.1" "LT608236.1" "LT608229.1" "LT608235.1" "LT608248.1"
Taxon_16 Mus_musculus "CM000994.3" "CM000995.3" "CM000996.3" "CM000997.3" "CM000998.3" "CM000999.3" "CM001000.3" "CM001001.3" "CM001002.3" "CM001003.3" "CM001004.3" "CM001005.3" "CM001006.3" "CM001007.3" "CM001008.3" "CM001009.3" "CM001010.3" "CM001011.3" "CM001012.3"
Taxon_17 Mus_pahari "LT608296.1" "LT608286.1" "LT608290.1" "LT608287.1" "LT608292.1" "LT608307.1" "LT608288.1" "LT608301.1" "LT608291.1" "LT608289.1" "LT608304.1" "LT608305.1" "LT608299.1" "LT608308.1" "LT608295.1" "LT608294.1" "LT608298.1" "LT608302.1" "LT608297.1" "LT608293.1" "LT608306.1" "LT608303.1" "LT608309.1"
Taxon_18 Mus_spretus "OW971678.1" "OW971679.1" "OW971680.1" "OW971682.1" "OW971684.1" "OW971683.1" "OW971685.1" "OW971687.1" "OW971688.1" "OW971686.1" "OW971689.1" "OW971692.1" "OW971691.1" "OW971690.1" "OW971693.1" "OW971694.1" "OW971695.1" "OW971696.1" "OW971697.1"
Taxon_19 Rattus_norvegicus "CM070391.1" "CM070392.1" "CM070393.1" "CM070394.1" "CM070395.1" "CM070396.1" "CM070397.1" "CM070398.1" "CM070399.1" "CM070400.1" "CM070401.1" "CM070402.1" "CM070403.1" "CM070404.1" "CM070405.1" "CM070406.1" "CM070407.1" "CM070408.1" "CM070409.1" "CM070410.1"
Taxon_20 Peromyscus_maniculatus "CM010879.2" "CM010880.2" "CM010881.2" "CM010882.2" "CM010883.2" "CM010884.2" "CM010885.2" "CM010886.1" "CM010887.2" "CM010888.2" "CM010889.2" "CM010890.2" "CM010891.2" "CM010892.2" "CM010893.2" "CM010894.2" "CM010895.2" "CM010896.2" "CM010897.2" "CM010898.1" "CM010899.2" "CM010900.2" "CM010901.2"
Taxon_21 Gorilla_gorilla "CM055446.2" "CM068950.1" "CM055449.2" "CM055450.2" "CM055451.2" "CM055452.2" "CM055453.2" "CM055454.2" "CM055455.2" "CM055456.2" "CM055457.2" "CM055458.2" "CM055459.2" "CM055460.2" "CM055461.2" "CM055462.2" "CM055463.2" "CM055464.2" "CM055465.2" "CM055466.2" "CM055467.2" "CM055468.2" "CM068951.1"
Taxon_22 Macaca_fascicularis "BLPH02000001.1" "BLPH02000002.1" "BLPH02000003.1" "BLPH02000004.1" "BLPH02000005.1" "BLPH02000006.1" "BLPH02000007.1" "BLPH02000008.1" "BLPH02000009.1" "BLPH02000010.1" "BLPH02000011.1" "BLPH02000012.1" "BLPH02000013.1" "BLPH02000014.1" "BLPH02000015.1" "BLPH02000016.1" "BLPH02000017.1" "BLPH02000018.1" "BLPH02000019.1" "BLPH02000020.1"
Taxon_23 Macaca_mulatta "CM014336.1" "CM014337.1" "CM014338.1" "CM014339.1" "CM014340.1" "CM014341.1" "CM014342.1" "CM014343.1" "CM014344.1" "CM014345.1" "CM014346.1" "CM014347.1" "CM014348.1" "CM014349.1" "CM014350.1" "CM014351.1" "CM014352.1" "CM014353.1" "CM014354.1" "CM014355.1"
Taxon_24 Pan_troglodytes "CM054434.2" "CM068905.1" "CM054437.2" "CM054438.2" "CM054439.2" "CM054440.2" "CM054441.2" "CM054442.2" "CM054443.2" "CM054444.2" "CM054445.2" "CM054446.2" "CM054447.2" "CM054448.2" "CM054449.2" "CM054450.2" "CM054451.2" "CM054452.2" "CM054453.2" "CM054454.2" "CM054455.2" "CM054456.2" "CM068906.1"
Taxon_25 Papio_anubis "CM018180.1" "CM018181.1" "CM018182.1" "CM018183.2" "CM018184.2" "CM018185.2" "CM018186.2" "CM018187.2" "CM018188.1" "CM018189.1" "CM018190.2" "CM018191.2" "CM018192.2" "CM018193.1" "CM018194.1" "CM018195.2" "CM018196.1" "CM018197.1" "CM018198.1"
Taxon_26 Piliocolobus_tephrosceles "CM019240.1" "CM019241.1" "CM019242.1" "CM019243.1" "CM019244.1" "CM019245.1" "CM019246.1" "CM019247.1" "CM019248.1" "CM019249.1" "CM019250.1" "CM019251.1" "CM019252.1" "CM019253.1" "CM019254.1" "CM019255.1" "CM019256.1" "CM019257.1" "CM019258.1" "CM019259.1" "CM019260.1"
```

After the configuration file is ready we make the downloading shell script executable and then run it. Note: for the script to work the Entrez Direct UNIX E-utilities needs to be in your $PATH.

```
chmod +x chromosome_download.sh
```

```
./chromosome_download.sh
```

Depending on the number of taxa you are using, this may download a substantial amount of data. It may take up a while to complete but progress updates will be sent from the script as each taxon is processed for you to track progress. The tutorial mammal dataset (26 species) takes **~6.25 hours** to download using 8 CPUS and 40G of RAM.   

</details>


## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>
  
>Landmarks can be any single-copy, conserved sequence that can be aligned to chromosomes in your dataset, but we have used ultraconserved elements (UCEs) in this tutorial as an example. More information about tetrapod UCEs can be found [here](https://www.ultraconserved.org/)

We will use the mapping shell script ```landmark_mapping.sh``` to identfy the location of different landmarks on the various chromosomes. We will make the script executable and then run it. Note: BWA and samtools need to be in your $PATH for the script to work.  

```
chmod +x landmark_mapping.sh
```

```
./landmark_mapping.sh
```

This script will take a while to execute, but like the download script, it will report progress along the way to keep you updated. Allowing 16 CPUS and 40GB of RAM the tutorial dataset takes **~17 hours** to be processed. When the script is done running each Taxon directory should have summary files called ```Total_UCE_counts.txt``` and ```Chromosome_lengths.txt``` which summarize the number of UCEs mapped to each chromosome and the length (in bp) of each chromosome, respectively. There should also be new directories for archiving the FASTA and SAMTOOLS files. Finally, the parent directory should have an R input file called ```sample_input_pres_abs.csv``` which will be used in the next step.

</details>

## Step 3: Cluster chromosomes based on landmark similarity

<details>
  <summary>Click to expand content!</summary>

>MMDS in R statistical software is used to identify which chromosomes likely contain homologous blocks of genomes (i.e. supergenes, Marian fragments etc.). 

After the last step we should now have a 4.7 MB file called ```sample_input_pres_abs.csv``` which is used in R to identify which chromosomes have many landmarks in common. The first few lines of this file should look like: 

```
chromosomes,uces
Bos_indicus_CM003021.1,uce-95_p1
Bos_indicus_CM003021.1,uce-110_p1
Bos_indicus_CM003021.1,uce-117_p2
Bos_indicus_CM003021.1,uce-127_p1
Bos_indicus_CM003021.1,uce-153_p1
Bos_indicus_CM003021.1,uce-197_p1
Bos_indicus_CM003021.1,uce-232_p2
Bos_indicus_CM003021.1,uce-264_p4
Bos_indicus_CM003021.1,uce-279_p1
Bos_indicus_CM003021.1,uce-286_p1
Bos_indicus_CM003021.1,uce-288_p1
Bos_indicus_CM003021.1,uce-319_p1
Bos_indicus_CM003021.1,uce-322_p1
```

You can execute the R script, ```chromosome_cluster.R``` to run all at once, but for the purpose of the tutorial we will walk through the key steps here. Start with loading the following R libraries:

```
library(ggplot2)
library(reshape2)
```

Next, we format the input file for presence-absence analysis of landmarks:

```
Input_pres_abs <- as.data.frame(read.csv("sample_input_pres_abs.csv", stringsAsFactors = F))
matrix <-dcast(Input_pres_abs, chromosomes ~ uces, length)
write.csv(matrix, file = 'landmark_pres_abs.csv')

```

These commands will produce a 5.6 MB CSV file ```landmark_pres_abs.csv``` that contains the presence-absence of each landmark. Then we extract the relevant columns for the MMDS:

```
X <- matrix[3:5113]
```

>NOTE: The tutorial dataset has a total of 5,110 mapped landmarks, this value will differ when using other input data and the user will need to calculate the number of non-string columns to proceed.  

```
d <- dist(X)
fit <- cmdscale(d,eig=TRUE, k=2)
temp <-cbind(matrix[1],fit$points)
colnames(temp)[2]  <- "Axis1"
colnames(temp)[3]  <- "Axis2"
write.csv(temp, file = 'MMDS_sample.csv')
```

We then plot the results using: 

```
data <- read.csv("MMDS_sample.csv", header =T, row.names = 1)
ggplot(data, aes(x=Axis1, y=Axis2)) + geom_point() + geom_text(size=4,label=rownames(data),check_overlap = F) + xlab("Chromosome Landmark Similarity Axis 1") + ylab("Chromosome Landmark Similarity Axis 2") + theme_classic()
```

Using the tutorial data this should produce a plot that looks like this: 

![Landmarks-1](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/Landmarks-1.jpg)

Data points in this plot represent 567 different chromosomes from the species we used in the pipeline. Because we used presence/absence of landmarks as the input data for this analysis, the placement of different data points should largely correspond to how many landmarks they share.We use the MMDS scores to identify those chromosomes that should have many UCE landmarks in common, indicating they likely contain homologous genomic regions. 

**Identifying chromosome sets for Genomic Disparity Analysis**

We need to set thresholds for landmark similarity in order to identify chromosomes that will be analyzed together. For example, we can see in the plot above that on Axis 1 there is a cluster of chromosomes with MMDS scores greater than 10. We can use this threshold to extract the names of these chromosomes from the ```MMDS_sample.csv file```. We can open the CSV file (using excel or similar program) and sort the Axis 1 scores from largest to smallest to identify the following chromosomes as belonging to this Chromosome Set, which we will call **Chromosome Set 1**.  

Chromosome ID | Species + Chromosome GenBank  | Axis 1 Score  
------------ | -------------  | -------------
250	| Felis_catus_CM031419.1	| 18.45603982
479	| Panthera_tigris_CM031438.1	| 18.45352839
353	| Mus_musculus_CM000995.3	| 16.76595889
181	| Cricetulus_griseus_CM023440.1	| 16.73988689
393	| Mus_spretus_OW971679.1	| 16.73905454
510	| Peromyscus_maniculatus_CM010882.2	| 16.64799578
348	| Mus_caroli_LT608244.1	| 16.64046105
374	| Mus_pahari_LT608290.1	| 16.62978242
426	| Ovis_aries_CM028705.1	| 16.62075814
461	| Pan_troglodytes_CM054447.2	| 16.50669384
551	| Rattus_norvegicus_CM070393.1	| 16.49902926
282	| Gorilla_gorilla_CM055457.2	| 16.45284339
326	| Macaca_mulatta_CM014347.1	| 16.40570999
498	| Papio_anubis_CM018189.1	| 16.39394662
539	| Piliocolobus_tephrosceles_CM019250.1	| 16.38174391
306	| Macaca_fascicularis_NW_025540829.1	| 16.36260442
262	| Giraffa_tippelskirchi_CM018105.1	| 16.26180669
31	| Bos_taurus_CM008169.2	| 16.20595656
110	| Capra_hircus_CM004563.1	| 16.20479132
412	| Neomonachus_schauinslandi_CM035899.1	| 16.16641741
59	| Bubalus_bubalis_CM034272.1	| 16.12303766
2	| Bos_indicus_CM003022.1	| 15.69361074
188	| Equus_asinus_CM027693.2	| 15.42305316
217	| Equus_caballus_CM027693.2	| 15.42305316
154	| Ceratotherium_simum_CM043826.1	| 15.34565128
82	| Capra_aegagrus_CM003215.1	| 14.51796749

The next closest score on Axis 1 is 2.65, so we will call these 26 chromosomes (one for each species in the dataset) **Chromosome Set 1**. Based on the MMDS result, it is clear that the chromosomes contained in **Chromosome Set 1** share many landmarks suggesting these chromosomes contain homologous genomic regions. However, identifying an MMDS score threshold is not always this clear. For example, althought it also includes one chromosome for each of the 26 mammal species, **Chromosome Set 2** described in Mohan et al. (2024) has a much narrower gap based on an Axis 2 MMDS score threshold (scores less than -8). When using other datasets, users are encouraged to experiment with chromosome set thresholds to determine how robust downstream results are. The plot below shows the two chromosome sets and the threshold values used:

![Landmarks-2](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/Landmarks-2.jpg)

For the rest of the tutorial we will focus on processing **Chromosome Set 1**. Following the procedure described above, you should have identified the following chromosomes as belonging to **Chromosome Set 1**:

Species | Order  | Chromosome (in MMDS set 1) 
------------ | -------------  | ------------- 
_Bos indicus_	| Artiodactyla | 2 (CM003022.1)
_Bos taurus_	| Artiodactyla | 2 (CM008169.2)
_Bubalus bubalis_	| Artiodactyla | 2 (CM034272.1)
_Capra aegagrus_	| Artiodactyla | 2 (CM003215.1)
_Capra hircus_	| Artiodactyla | 2 (CM004563.1)
_Giraffa tippelskirchi_	| Artiodactyla | 3 (CM018105.1)
_Ovis aries_	| Artiodactyla | 2 (CM028705.1)
_Ceratotherium simum_	| Perissodactyla | 9 (CM043826.1)
_Equus asinus_ | Perissodactyla | 4 (CM027693.2)
_Equus caballus_	| Perissodactyla | 18 (CM027693.2)
_Felis catus_ | Carnivora | C1 (CM031419.1)
_Neomonachus schauinslandi_	| Carnivora | 3 (CM035899.1)
_Panthera tigris_ | Carnivora | C1 (CM031438.1)
_Cricetulus griseus_ | Rodentia | 6 (CM023440.1)
_Mus caroli_ | Rodentia | 2 (LT608244.1)
_Mus musculus_ | Rodentia | 2 (CM000995.3)
_Mus pahari_ | Rodentia | 3 (LT608290.1)
_Mus spretus_	| Rodentia | 2 (OW971679.1)
_Rattus norvegicus_	| Rodentia | 3 (CM070393.1)
_Peromyscus maniculatus_ | Rodentia | 4 (CM010882.2)
_Gorilla gorilla_	| Primates | 3 (CM055457.2)
_Macaca fascicularis_ | Primates | 12 (NW_025540829.1)
_Macaca mulatta_	| Primates | 12 (CM014347.1)
_Pan troglodytes_ | Primates | 2B (CM054447.2)
_Papio anubis_ | Primates | 10 (CM018189.1)
_Piliocolobus tephrosceles_ | Primates | 11 (CM019250.1)

**Now that we have identified the 26 chromosomes we will analyze, they will be further prepped in the next step.** 

</details>

## Step 4: Preparing landmarked chromosomes for downstream analysis

<details>
  <summary>Click to expand content!</summary>

>We need to remove any landmarks that the chromosomes do not share, check the directionality/orientation of the chromosomes, and merge UCE landmarks before Genomic Disparity Analysis 
  
Now that we have identified the 26 chromosomes that comprise **Chromosome Set 1**, we will extract them for further analysis using R. We will do this using the ```landmark_pres_abs.csv``` file generated during the previous step. 

```
library(reshape2)
data1 <- read.csv("landmark_pres_abs.csv", header =T, row.names = 1)
data2 <- t(data1)
write.csv(data2, file = 'landmark_exclusion_matrix.csv')
data3 <- read.csv("landmark_exclusion_matrix.csv", skip = 1)
```

Now we will subset the columns that correspond to the UCE landmarks + the 26 chromosomes of interest:

```
data4 <-data3[c("chromosomes", "Felis_catus_CM031419.1", "Panthera_tigris_CM031438.1", "Mus_musculus_CM000995.3", "Cricetulus_griseus_CM023440.1", "Mus_spretus_OW971679.1", "Peromyscus_maniculatus_CM010882.2", "Mus_caroli_LT608244.1", "Mus_pahari_LT608290.1", "Ovis_aries_CM028705.1", "Pan_troglodytes_CM054447.2", "Rattus_norvegicus_CM070393.1", "Gorilla_gorilla_CM055457.2", "Macaca_mulatta_CM014347.1", "Papio_anubis_CM018189.1", "Piliocolobus_tephrosceles_CM019250.1", "Macaca_fascicularis_NW_025540829.1", "Giraffa_tippelskirchi_CM018105.1", "Bos_taurus_CM008169.2", "Capra_hircus_CM004563.1", "Neomonachus_schauinslandi_CM035899.1", "Bubalus_bubalis_CM034272.1", "Bos_indicus_CM003022.1", "Equus_asinus_CM027693.2", "Equus_caballus_CM027693.2", "Ceratotherium_simum_CM043826.1", "Capra_aegagrus_CM003215.1")]

write.csv(data4, file = 'cluster1_exclusion_matrix.csv')
```

Next we will remove all of the landmarks (UCE probes) that are not shares across all of the mammal species: 

```
chr <- data4
chr_matchedUCEs <- apply(chr, 1, function(row) all(row != 0))
chr_clean <- chr[chr_matchedUCEs,]
write.csv(chr_clean, file = "present_landmarks.csv")
```

>This should result in an output file that has 220 UCE probes that were present (=1) in all 16 species. 

Now we will use these to extract position for each probe on each chromosome from the SAM files generated during **Step 2**. We will use the shell script ```chromosome_retriever.sh``` to collect the necessary SAM files (now in TSV format). We will make the script executable and then run it in the same directory that the previous shell scripts were executed in. 

```
chmod +x chromosome_retriever.sh
```

```
./chromosome_retriever.sh
```

This script will create a directory called ```chromosome_set``` that contains all 26 chromosomes from **Chromosome Set 1**. We will access it in the following steps using R. Now that we have the relevant chromosome files isolated, we need to extract the landmark, direction and position information from each file. 

```
library(matrixStats)
library(dplyr)

homologousUCE <- read.csv("present_landmarks.csv", header =T, row.names = 1)

folder_path <- "./chromosome_set"
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

After this step you will have a CSV file called ```homologous_UCEs_extracted.csv``` which contains two columns for each species: (1) the 'V2' direction that landmarks were mapped on the chromosomes (0 [forward] or 16 [reverse]) and (2) the mapping location (in base pairs) of each landmark. In the tutorial dataset there should be several examples of species with chromsomes that have been assembled with opposite complementarities. For example, if we look at the first 20 landmark positions for _Capra aegragrus_ and _Capra hircus_ in the ```homologous_UCEs_extracted.csv``` file we should see this:   

V2.y.y | Capra_aegagrus_CM003215.1.fasta	  | V2.x.x.x | Capra_hircus_CM004563.1.fasta
------------ | -------------  | ------------- | -------------	        	
16	 | 50568525	   |                             0	 |       84091246
16	 |   50568465	 |                             0	 |      84091306
16	 |   50568405	 |                            0	   |     84091366
16	 |   36102876  |                            	0	 |       99207598
16	 |   38905700	 |                             0	 |       96044026
16	 |   50308872	 |                             0	 |       84336628
0	   |   18829966	 |                             16	 |     114936983
0	   |   18830026	 |                             16	 |     114936923
16	 |   50555754	 |                             0	 |       84103333
16	 |   51538534	 |                             0   |     	83172024
16	 |   50991589	 |                             0   |     	83702175
16	 |   51789652	 |                             0	 |       82921983
0	   |   30815279	 |                             16	 |     104073711
16	 |   38900786	 |                             0	 |       96048941
0	   |   31036783	 |                             16	 |     103852106
16	 |   49929220	 |                             0	 |       84716751
16   | 	49778672	 |                             0	 |       84869062
0	   |  29898324	 |                             16	 |     104364483
16	 |   50992659	 |                             0	 |       83701105
16	 |   50992599	 |                             0	 |       83701165

We can see that the direction of the landmark placement is opposite in all cases which suggests that these chromosomes have been submitted to NCBI with opposite complementarities. Another way to see evidence of this is to visualize the placement of the landmarks on the chromosomes where we can easily see the need to 'flip' the chromsomes for some taxa so that the landmark positions are homologous in an evolutionary sense. In the image below UCE landmarks are indicated in orange on the grey chromosomes, homologous UCE landmarks are connected with a light blue line. 

![Capra](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/Capra.jpg)

If we take the reverse complement of all the landmark positions for _Capra aegragrus_ we see that the landmark positions are now oriented in the same direction. 

![Capra_flipped](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/Capra_flipped.jpg)

Determining which taxa to 'flip' is arbitrary for disparity analysis, but for the tutorial dataset we will 'flip' those taxa with the same landmark directionalities as _Capra aegragrus_. That includes seven other species: _Bos indicus_, _Bos taurus_, _Bubalus bubalus_, _Cricetulus griseus_, _Equus asinus_, _Equus caballus_, and _Ovis aries_. In order to 'flip' taxa we can subtract the total length of the chromosome from the existing BWA-inferred landmark positions. We can find the total lengths of of the chromosomes in the ```Chromosome_lengths.tsv``` files we generated during **Step 2**.

Species | Accession	  | Length (bp)
------------ | -------------  | -------------     	
_Bos indicus_ | CM003022.1 | 140680885
_Bos taurus_ | CM008169.2 | 136231102
_Bubalus bubalus_ | CM034272.1 | 188164321
_Capra aegragrus_	 | CM003215.1	| 128023632
_Cricetulus griseus_ | CM023440.1 | 155611870
_Equus asinus_ | CM027693.2 | 92916174
 _Equus caballus_ | CM027693.2 | 82641348
 _Ovis aries_ | CM028705.1 | 250202058

This 'flipping' step can be conducted in Excel or similar spreadsheet editor. We can also use R to 'flip' these chromosomes using function commands: 

```
landmarkflip <- read.csv("homologous_UCEs_extracted.csv", header =T, row.names = 1)

fun1 <- function(x) {140680885-x}
fun2 <- function(x) {136231102-x}
fun3 <- function(x) {188164321-x}
fun4 <- function(x) {128023632-x}
fun5 <- function(x) {155611870-x}
fun6 <- function(x) {92916174-x}
fun7 <- function(x) {82641348-x}
fun8 <- function(x) {250202058-x}

Boin <-lapply(landmarkflip$Bos_indicus_CM003022.1.fasta, fun1)
Bota <-lapply(landmarkflip$Bos_taurus_CM008169.2.fasta, fun2)
Bubu <-lapply(landmarkflip$Bubalus_bubalis_CM034272.1.fasta, fun3)
Caae <-lapply(landmarkflip$Capra_aegagrus_CM003215.1.fasta, fun4)
Crgr <-lapply(landmarkflip$Cricetulus_griseus_CM023440.1.fasta, fun5)
Eqas <-lapply(landmarkflip$Equus_asinus_CM027693.2.fasta, fun6)
Eqca <-lapply(landmarkflip$Equus_caballus_CM027693.2.fasta, fun7)
Ovar <-lapply(landmarkflip$Ovis_aries_CM028705.1.fasta, fun8)

landmarkflip$Bos_indicus_CM003022.1.fasta <- Boin
landmarkflip$Bos_taurus_CM008169.2.fasta <- Bota
landmarkflip$Bubalus_bubalis_CM034272.1.fasta <- Bubu
landmarkflip$Capra_aegagrus_CM003215.1.fasta <- Caae
landmarkflip$Cricetulus_griseus_CM023440.1.fasta <- Crgr
landmarkflip$Equus_asinus_CM027693.2.fasta <- Eqas
landmarkflip$Equus_caballus_CM027693.2.fasta <- Eqca
landmarkflip$Ovis_aries_CM028705.1.fasta <- Ovar


A <-as.numeric(landmarkflip$Bos_indicus_CM003022.1.fasta)
B <-as.numeric(landmarkflip$Bos_taurus_CM008169.2.fasta)
C <-as.numeric(landmarkflip$Bubalus_bubalis_CM034272.1.fasta)
D <-as.numeric(landmarkflip$Capra_aegagrus_CM003215.1.fasta)
E <-as.numeric(landmarkflip$Cricetulus_griseus_CM023440.1.fasta)
F <-as.numeric(landmarkflip$Equus_asinus_CM027693.2.fasta)
G <-as.numeric(landmarkflip$Equus_caballus_CM027693.2.fasta)
H <-as.numeric(landmarkflip$Ovis_aries_CM028705.1.fasta)

landmarkflip$Bos_indicus_CM003022.1.fasta <- A
landmarkflip$Bos_taurus_CM008169.2.fasta <- B
landmarkflip$Bubalus_bubalis_CM034272.1.fasta <- C
landmarkflip$Capra_aegagrus_CM003215.1.fasta <- D
landmarkflip$Cricetulus_griseus_CM023440.1.fasta <- E
landmarkflip$Equus_asinus_CM027693.2.fasta <- F
landmarkflip$Equus_caballus_CM027693.2.fasta <- G
landmarkflip$Ovis_aries_CM028705.1.fasta <- H
```
 
>Note: In chromosomes with relatively conserved landmark placements, it should be obvious which taxa need to be 'flipped'. However, when landmarks are more evolutionarily labile it may be diffcult to justify a 'flipping' operation, so we encourage users to think about this operation carefully. 

Now all of the chromosomes are positioned correctly and we can remove the mapping direction information and export our final file:

```
data5 <- landmarkflip[-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)]

write.csv(data5, file = "homologous_UCEs_extracted_flipped.csv")
```

After the last step we have a file called ```homologous_UCEs_extracted_flipped.csv``` that is ready for the final preparation steps. This includes accounting for a specific caveat of using the UCE probe set. The UCE probe set was developed to capture UCEs across diverse taxa, as such some UCEs are targeted by multiple probes, so to control for the variation this creates in mapping, we average the probe placement across landmarks targeting the same UCE. We will remove the probe numbers (e.g. p1) in order to merge information from the probes targeting multiple parts of the same UCE landmark. After we average the UCE positions, we transpose the matrix to prepare it for the PCA: 

```
data6 <- read.csv("homologous_UCEs_extracted_flipped.csv")
data6$landmarks <- sub("_.*", "", data6$landmarks)
data6$X <-NULL

data7 <-aggregate(data6, by = list(data6$landmarks), mean)
data7$landmarks<-NULL
colnames(data7)[colnames(data7) == "Group.1"] <- "chromosomes"

transposed_data <- data7 %>% t() %>% as.data.frame()

write.csv(transposed_data, file = "homologous_UCEs_Set1_PCA.csv")
```
  
After averaging the positions, there should be **186 UCE landmarks**. This final prep will produce the file ```homologous_UCEs_cluster1_PCA.csv``` which is the input file for Genomic Disparity Analysis as outlined in the final step of the tutorial. 

</details>

## Step 5: Principal components analysis of landmark disparity

<details>
  <summary>Click to expand content!</summary>

>In this final step we will visualize the disparity in landmark placement using Principal Components Analysis (PCA). We will specifically focus on PCs 1, 2 and 3 from this analysis. 

At the end of **Step 4** we will have the file ```homologous_UCEs_Set1_PCA.csv```, which we will now load into R for the PCA:  

```
library(dplyr)
library(ggplot2)

Set1_PCA <- read.csv("homologous_UCEs_Set1_PCA.csv", row.names = 1)
names(Set1_PCA) <- Set1_PCA[1,]
Set1_PCA <- Set1_PCA[-1,]

PCprep <- Set1_PCA %>% mutate_at(1:186, as.numeric)

PCA1 <-prcomp(PCprep)

```

After the PCA has completed there will be five lists of results including "sdev", "rotation", "center", "scale", and "x". The "x" variable contains the PC scores that we will use for Genomic Disparity analysis. We will now extract them and prepare for visualizing the results in ggplot2.

```
PCA1_scores <-as.data.frame(PCA1$x)

gplot(PCA1_scores, aes(x = PC1, y = PC2, label = row.names(PCA1_scores))) + geom_point(size = 3) + geom_text(size = 3)  + theme_classic()
gplot(PCA1_scores, aes(x = PC2, y = PC3, label = row.names(PCA1_scores))) + geom_point(size = 3) + geom_text(size = 3)  + theme_classic()

```

![PC1_PC2_dirty](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/PC1_PC2_dirty.jpg)

![PC2_PC3_dirty](https://github.com/nhm-herpetology/genomic-disparity/blob/main/Mammal-uces/PC2_PC3_dirty.jpg)

Using these two plots we can see the variation in UCE landmark placement across the three most explanatory axes. For clearer interpretation, we can tidy up the plots by adding information about the mammalian orders, and making the plots again. 

```
Orders <-c("Artiodactlya", "Artiodactlya", "Artiodactlya", "Artiodactlya", "Artiodactlya", "Perissodactyla", "Rodentia", "Perissodactyla", "Perissodactyla", "Carnivora", "Artiodactlya", "Primates", "Primates", "Rodentia", "Rodentia", "Rodentia", "Rodentia", "Carnivora", "Artiodactlya", "Primates", "Carnivora", "Primates", "Rodentia", "Primates", "Rodentia")

PCA1_plots <-cbind(PCA1_scores, Orders)

ggplot(PCA1_plots, aes(x = PC1, y = PC2, color = Orders)) + geom_point(size = 3) + theme_classic()

ggplot(PCA1_plots, aes(x = PC2, y = PC3, color = Orders)) + geom_point(size = 3) + theme_classic()

```

</details>
