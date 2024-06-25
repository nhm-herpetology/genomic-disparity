#R commands for Genomic Disparity Analysis Mammal-uce tutorial 
#See GitHub for more details

#Chromosome clustering commands

library(ggplot2)
library(reshape2)
library(matrixStats)
library(dplyr)
library(cowplot)

Input_pres_abs <- as.data.frame(read.csv("sample_input_pres_abs.csv", stringsAsFactors = F))
matrix <-dcast(Input_pres_abs, chromosomes ~ uces, length)
write.csv(matrix, file = 'landmark_pres_abs.csv')

X <- matrix[3:5113]

d <- dist(X)
fit <- cmdscale(d,eig=TRUE, k=2)
temp <-cbind(matrix[1],fit$points)
colnames(temp)[2]  <- "Axis1"
colnames(temp)[3]  <- "Axis2"
write.csv(temp, file = 'MMDS_sample.csv')

data <- read.csv("MMDS_sample.csv", header =T, row.names = 1)
ggplot(data, aes(x=Axis1, y=Axis2)) + geom_point() + geom_text(size=4,label=rownames(data),check_overlap = F) + xlab("Chromosome Landmark Similarity Axis 1") + ylab("Chromosome Landmark Similarity Axis 2") + theme_classic()

#Chromosome preparation steps

library(reshape2)
data1 <- read.csv("landmark_pres_abs.csv", header =T, row.names = 1)
data2 <- t(data1)
write.csv(data2, file = 'landmark_exclusion_matrix.csv')
data3 <- read.csv("landmark_exclusion_matrix.csv", skip = 1)

data4 <-data3[c("chromosomes", "Felis_catus_CM031419.1", "Panthera_tigris_CM031438.1", "Mus_musculus_CM000995.3", "Cricetulus_griseus_CM023440.1", "Mus_spretus_OW971679.1", "Peromyscus_maniculatus_CM010882.2", "Mus_caroli_LT608244.1", "Mus_pahari_LT608290.1", "Ovis_aries_CM028705.1", "Pan_troglodytes_CM054447.2", "Rattus_norvegicus_CM070393.1", "Gorilla_gorilla_CM055457.2", "Macaca_mulatta_CM014347.1", "Papio_anubis_CM018189.1", "Piliocolobus_tephrosceles_CM019250.1", "Macaca_fascicularis_NW_025540829.1", "Giraffa_tippelskirchi_CM018105.1", "Bos_taurus_CM008169.2", "Capra_hircus_CM004563.1", "Neomonachus_schauinslandi_CM035899.1", "Bubalus_bubalis_CM034272.1", "Bos_indicus_CM003022.1", "Equus_asinus_CM027693.2", "Equus_caballus_CM027693.2", "Ceratotherium_simum_CM043826.1", "Capra_aegagrus_CM003215.1")]

write.csv(data4, file = 'cluster1_exclusion_matrix.csv')

chr <- data4
chr_matchedUCEs <- apply(chr, 1, function(row) all(row != 0))
chr_clean <- chr[chr_matchedUCEs,]
write.csv(chr_clean, file = "present_landmarks.csv")

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

data5 <- landmarkflip[-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)]

write.csv(data5, file = "homologous_UCEs_extracted_flipped.csv")

data6 <- read.csv("homologous_UCEs_extracted_flipped.csv")
data6$landmarks <- sub("_.*", "", data6$landmarks)
data6$X <-NULL

data7 <-aggregate(data6, by = list(data6$landmarks), mean)
data7$landmarks<-NULL
colnames(data7)[colnames(data7) == "Group.1"] <- "chromosomes"

transposed_data <- data7 %>% t() %>% as.data.frame()

write.csv(transposed_data, file = "homologous_UCEs_Set1_PCA.csv")

#PCA and plotting commands

Set1_PCA <- read.csv("homologous_UCEs_Set1_PCA.csv", row.names = 1)
names(Set1_PCA) <- Set1_PCA[1,]
Set1_PCA <- Set1_PCA[-1,]

PCprep <- Set1_PCA %>% mutate_at(1:186, as.numeric)

PCA1 <-prcomp(PCprep)

PCA1_scores <-as.data.frame(PCA1$x)

Order <-c("Artiodactyla", "Artiodactyla", "Artiodactyla", "Artiodactyla", "Artiodactyla", "Perissodactyla", "Rodentia", "Perissodactyla", "Perissodactyla", "Carnivora", "Artiodactyla", "Primates", "Primates", "Rodentia", "Rodentia", "Rodentia", "Rodentia", "Carnivora", "Artiodactyla", "Primates", "Carnivora", "Primates", "Rodentia", "Primates", "Rodentia")

PCA1_plots <-cbind(PCA1_scores, Order)

P1 <-ggplot(PCA1_plots, aes(x = PC1, y = PC2, color = Order)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("Artiodactyla","Carnivora","Perissodactyla","Primates","Rodentia"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

P2 <-ggplot(PCA1_plots, aes(x = PC3, y = PC4, color = Order)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("Artiodactyla","Carnivora","Perissodactyla","Primates","Rodentia"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

plot_grid(P1, P2, ncol = 1)

write.csv(PCA1$rotation, file = "component_loadings_PCA.csv")

cent1 <- aggregate(cbind(PC1, PC2) ~Order, data = PCA1_plots, FUN = mean)
cent2 <- aggregate(cbind(PC3, PC4) ~Order, data = PCA1_plots, FUN = mean)
segs1 <- merge(PCA1_plots, setNames(cent1, c('Order','oPC1','oPC2')), by = 'Order', sort = FALSE)
segs2 <- merge(PCA1_plots, setNames(cent2, c('Order','oPC3','oPC4')), by = 'Order', sort = FALSE)

out1 <-ggplot(PCA1_plots, aes(x = PC1, y = PC2, color = Order)) + geom_segment(data = segs1, mapping = aes(xend = oPC1, yend = oPC2)) + geom_point(data = cent1, size = 5) + geom_point() + coord_fixed() + xlab('PC1 (98.9%)') + ylab('PC2 (1.0%)') + scale_color_manual(breaks = c("Artiodactyla", "Carnivora","Perissodactyla","Primates","Rodentia"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

out2 <-ggplot(PCA1_plots, aes(x = PC3, y = PC4, color = Order)) + geom_segment(data = segs2, mapping = aes(xend = oPC3, yend = oPC4)) + geom_point(data = cent2, size = 5) + geom_point() + coord_fixed() + xlab('PC3 (0.3%)') + ylab('PC4 (0.03%)') + scale_color_manual(breaks = c("Artiodactyla", "Carnivora","Perissodactyla","Primates","Rodentia"), values=c("orange", "red","brown","blue","purple")) + theme_classic()

plot_grid(out1, out2, ncol = 1)