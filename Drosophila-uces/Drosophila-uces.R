#R commands for Genomic Disparity Analysis Drosophila-uce tutorial
#See GitHub page for more details

library(reshape2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(cowplot)

Input_pres_abs <- as.data.frame(read.csv("sample_input_pres_abs.csv", stringsAsFactors = F))
matrix <-dcast(Input_pres_abs, chromosomes ~ uces, length)
write.csv(matrix, file = 'landmark_pres_abs.csv')

data1 <- read.csv("landmark_pres_abs.csv", row.names = 1)
chr <- data1
chr_matchedUCEs <- apply(chr, 1, function(row) all(row != 0))
chr_clean <- chr[chr_matchedUCEs,]
write.csv(chr_clean, file = "present_landmarks.csv")

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

landmarkflip <- read.csv("homologous_UCEs_extracted.csv", row.names = 1)
fun1 <- function(x) {27217941-x}
Drfl <-lapply(landmarkflip$mapped_D_flavomontana_c5.fasta, fun1)
landmarkflip$mapped_D_flavomontana_c5.fasta <- Drfl
A <-as.numeric(landmarkflip$mapped_D_flavomontana_c5.fasta)
landmarkflip$mapped_D_flavomontana_c5.fasta <-A

data2 <- landmarkflip[-c(2,4,6,8,10,12)]

write.csv(data2, file = "homologous_UCEs_extracted_flipped.csv")

data3 <- read.csv("homologous_UCEs_extracted_flipped.csv")
data3$chromosomes <- sub("_.*", "", data3$chromosomes)
data3$X <-NULL

data4 <-aggregate(data3, by = list(data3$chromosomes), mean)
data4$chromosomes<-NULL
colnames(data4)[colnames(data4) == "Group.1"] <- "chromosomes"

transposed_data <- data4 %>% t() %>% as.data.frame()

write.csv(transposed_data, file = "homologous_UCEs_chromosome_5_PCA.csv")

#PCA commands and plot commands

C5_PCA <- read.csv("homologous_UCEs_chromosome_5_PCA.csv", row.names = 1)
names(C5_PCA) <- C5_PCA[1,]
C5_PCA <- C5_PCA[-1,]

PCprep <- C5_PCA %>% mutate_at(1:191, as.numeric)

PCAC5 <-prcomp(PCprep)

PCAC5_scores <-as.data.frame(PCAC5$x)

Species <-c("americana", "flavomontana", "montana", "novamexicana", "virilis", "virilis")

PCAC5_plots <-cbind(PCAC5_scores, Species)

P1 <-ggplot(PCAC5_plots, aes(x = PC1, y = PC2, color = Species)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("americana", "flavomontana", "montana", "novamexicana", "virilis"), values=c("orange", "pink2","darkred","darkgreen","purple")) + theme_classic()

P2 <-ggplot(PCAC5_plots, aes(x = PC3, y = PC4, color = Species)) + geom_point(size = 4, alpha=0.9) + scale_color_manual(breaks = c("americana", "flavomontana", "montana", "novamexicana", "virilis"), values=c("orange", "pink2","darkred","darkgreen","purple")) + theme_classic()

plot_grid(P1, P2, ncol = 1)

