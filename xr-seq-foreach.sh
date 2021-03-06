#!/usr/bin/env bash


	jq -c '.data[] | select(.protocol == "XR-seq") | .title' sampleList.json | while read -r item; do

	echo $item

	genome=$(jq --raw-output '.data[] | select(.title == '$item') | .reference_genome' sampleList.json)
	srr_list=$(jq --raw-output '.data[] | select(.title == '$item') | "fastq-dump --stdout " + .SRA.runs[] + " >${SAMPLE}.fastq"' sampleList.json)
	echo $srr_list
	echo $genome
	
	SAMPLE=$item

	$srr_list #fastq-dump arragment
	GENOME_DIR=/data/genomes/$genome
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

done
