#!/usr/bin/env bash

read -p "Enter your Organism: " ORGANISM_NAME

if [ ${ORGANISM_NAME} == "human" ]
    then
    ORGANISM="Homo_sapiens"
    GENOME="GRCh38"
    FASTQ_FILE=human
elif [ ${ORGANISM_NAME} == "mouse" ]
    then
    ORGANISM="mus_musculus"
    GENOME="GRCm38"
    FASTQ_FILE=mouse
elif [ ${ORGANISM_NAME} == "plant" ]
    then
    ORGANISM="Arabidopsis_thaliana"
    GENOME="TAIR10"
    FASTQ_FILE=plant
elif [ ${ORGANISM_NAME} == "bacterium" ]
    then
    ORGANISM="Escherichia_coli_str_k_12_substr_mg1655"
    GENOME="ASM584v2"
    FASTQ_FILE=bacterium
elif [ ${ORGANISM_NAME} == "mushrooms" ]
    then
    ORGANISM="Saccharomyces_cerevisiae"
    GENOME="R64-1-1"
    FASTQ_FILE=mushrooms
else
    echo "You have entered incorrectly. Human, Mouse, Plant, Bacterium, Mushrooms"
    read
fi

SAMPLE=FASTQ_FILE # Change based on your file base name.
SAMPLE=CPD3 # Change based on your file base name.

# fastq-dump --stdout SRR1976056 >${SAMPLE}.fastq && fastq-dump --stdout SRR1976057 >>${SAMPLE}.fastq # CPD example. Comment out this line when you run the pipeline for your local file. 
# fastq-dump --stdout SRR3062635 >${SAMPLE}.fastq && fastq-dump --stdout SRR3062636 >>${SAMPLE}.fastq # (6-4)PP example. Comment out this line when you run the pipeline for your local file. 
fastq-dump --stdout SRR3062593 >${SAMPLE}.fastq && fastq-dump --stdout SRR3062594 >>${SAMPLE}.fastq  && fastq-dump --stdout SRR3062595 >>${SAMPLE}.fastq # (6-4)PP example. Comment out this line when you run the pipeline for your local file. 

GENOME_DIR=/data/genomes/${GENOME}
BOWTIE2_IND=${GENOME_DIR}/Bowtie2/genome

echo "Cut adapter"
cutadapt -a TGGAATTCTCGGGTGCCAAGG AACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq

echo "Align with the reference genome"
bowtie2 -p 4 -x $BOWTIE2_IND -U ${SAMPLE}_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam

echo "Convert to bam"
samtools view -q 20 -b -o ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.sam

echo "Convert to bed"
bedtools bamtobed -i ${SAMPLE}_cutadapt.bam >${SAMPLE}_cutadapt.bed

echo "Sort bed file. Use -u to remove duplicates"
sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE}_cutadapt.bed >${SAMPLE}_cutadapt_sorted.bed

echo "Count number of mapped reads"
grep -c "^" ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_readCount.txt 

echo "Get the read length distribution of the aligned and deduplicated reads"
awk '{print $3-$2}' ${SAMPLE}_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}'

echo "Get the certain-sized reads (eg 26)"
awk '{ if ($3-$2 == 26) { print } }' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_26.bed

echo "Retrieve sequences in fasta format"
bedtools getfasta -fi ${GENOME_DIR}/genome.fa -bed ${SAMPLE}_cutadapt_sorted_26.bed -fo ${SAMPLE}_cutadapt_sorted_26.fa -s

echo "Get the dinucleotide content of the reads"
fa2kmerAbundanceTable.py -i ${SAMPLE}_cutadapt_sorted_26.fa -k 2 --percentage -o ${SAMPLE}_cutadapt_sorted_26_dinucleotideTable.txt

echo "Split into strands"
awk '{if($6=="+"){print}}' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_plus.bed
awk '{if($6=="-"){print}}' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_minus.bed

echo "Convert bed to bedgraph"
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted_plus.bed -g ${GENOME_DIR}/genome.fa.fai -bg -scale $(cat ${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print 1000000/$1}') >${SAMPLE}_cutadapt_sorted_plus.bdg
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted_minus.bed -g ${GENOME_DIR}/genome.fa.fai -bg -scale $(cat ${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print -1000000/$1}') >${SAMPLE}_cutadapt_sorted_minus.bdg

echo "Convert bedgraph to bigwig"
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted_plus.bdg ${GENOME_DIR}/genome.fa.fai ${SAMPLE}_cutadapt_sorted_plus.bw
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted_minus.bdg ${GENOME_DIR}/genome.fa.fai ${SAMPLE}_cutadapt_sorted_minus.bw

echo "Count read values for transcribed (TS) and nontranscribed (NTS) strands"
bedtools intersect -sorted -a ${GENOME_DIR}/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -S -F 0.5 > ${SAMPLE}_cutadapt_sorted_TScount.txt
bedtools intersect -sorted -a ${GENOME_DIR}/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -s -F 0.5 > ${SAMPLE}_cutadapt_sorted_NTScount.txt

# chr5:168,478,195-168,593,555
# chr17:7,400,388-7,597,844