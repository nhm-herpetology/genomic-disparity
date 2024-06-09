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

Depending on the number of taxa you are using, this may download a substantial amount of data. It may take up a while to complete but progress updates will be sent from the script as each taxon is processed for you to track progress. The tutorial mammal dataset (26 species) takes XX hours to download using 8 CPUS and 40G of RAM.   

</details>


## Step 2: Reference mapping of landmark sequences
<details>
  <summary>Click to expand content!</summary>
  
>Landmarks can be any single-copy, conserved sequence that can be aligned to chromosomes in your dataset, but we have used ultraconserved elements (UCEs) in this tutorial as an example. More information about tetrapod UCEs can be found [here](https://www.ultraconserved.org/)

We will use the mapping shell script ```landmark_mapping.sh``` to identfy the location of different landmarks on the various chromosomes. We will make the script executable and then run it. Note: BWA and samtools need to be in your $PATH for the script to work.  

```
chmod +x chromosome_download.sh
```

```
./chromosome_download.sh
```


The UCE probe set was developed to capture UCEs across diverse taxa, as such some UCEs are targeted by multiple probes, so to control for the variation this creates in mapping, we average the probe placement across landmarks targeting the same UCE. 


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
