#!/bin/bash

for dir in */; do
	(
	echo "Starting mapping on $dir"
	cd $dir

wget https://raw.githubusercontent.com/nhm-herpetology/museum-NGS-training/main/Unit_03/Bioinformatics_Lab/Tetrapods-UCE-5Kv1.fasta

for i in $(ls *.fasta)
do
	bwa index $i
	bwa mem $i Tetrapods-UCE-5Kv1.fasta -t 6 > bwa_mem_align_UCEs_$i.sam 
	samtools view -S -b bwa_mem_align_UCEs_$i.sam > UCE_$i.bam 
	samtools sort UCE_$i.bam  -o UCE_$i.sorted.bam 
	samtools index UCE_$i.sorted.bam 
	samtools view -F 4 UCE_$i.bam > mapped_$i.sam 
	wc -l mapped_$i.sam  | tr "," "." > UCEcount_$i.csv
done

echo "Finished with mapping... tidying file structure"

rm mapped_Tetrapods-UCE-5Kv1.fasta.sam
rm UCE_Tetrapods-UCE-5Kv1.fasta.bam
rm UCE_Tetrapods-UCE-5Kv1.fasta.sorted.bam
rm UCE_Tetrapods-UCE-5Kv1.fasta.sorted.bam.bai
rm UCEcount_Tetrapods-UCE-5Kv1.fasta.csv

cat UCEcount_*.csv > Total_UCE_counts.txt

echo "Total UCE landmarks mapped to chromosomes in this assembly:"

echo "$(<Total_UCE_counts.txt )"

awk '(NR == 1) || (FNR == 1)' bwa_mem_align_UCEs_*.fasta.sam > Chromosomes_lengths.tsv
 
mkdir fasta_files
mkdir samtools_files
mv *.fasta fasta_files/
mv *.bam samtools_files/
mv *.bai samtools_files/
mv *.amb samtools_files/
mv *.ann samtools_files/
mv *.bwt samtools_files/
mv *.pac samtools_files/
mv *.sa samtools_files/
mv bwa_mem_align_*.sam samtools_files/
rm UCEcount_*.csv

for f in *.sam; do
    mv "$f" "$(basename "$f" .sam).tsv"
done

echo "Merging data for comparison with other taxa... " 

cut -f1,3 *.tsv > Landmarks_merged.tsv
sed 's/\t/,/g' Landmarks_merged.tsv > Landmarks_merged.csv
awk -F, -v OFS=, '{ print $2,$1 }' Landmarks_merged.csv > Landmarks_merged-temp.csv

VAR2="$(awk '{print $2}' Taxon_*.txt)"
sed -e 's/^/'$VAR2'_/' Landmarks_merged-temp.csv > Landmarks_merged-done.csv
mv -- Landmarks_merged-done.csv "${VAR2}_Landmarks.land" 

rm Landmarks_merged-temp.csv 
rm Landmarks_merged.csv

mv *.land ../

)done

cat *.land > All_taxa.land

sed ' 1 s/.*/chromosomes,uces/' All_taxa.land > sample_input_pres_abs.csv
