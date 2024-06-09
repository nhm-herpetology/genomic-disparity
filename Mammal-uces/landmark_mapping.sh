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

)done
