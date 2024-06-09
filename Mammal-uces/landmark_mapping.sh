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

cat UCEcount_*.csv > Total_UCE_counts.txt

)done
