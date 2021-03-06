#! /bin/bash
# inspired from https://github.com/milospjanic/bam2bigwig/blob/master/bam2bigwig.sh

module load bedtools
module load ucsc


bam=$1

bamToBed -i $bam -split > ${bam%.*}.accepted_hits.bed

cat ${bam%.*}.accepted_hits.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=Hg38_shRNA_hybrid.sizes stdin | sort -k1,1 -k2,2n > ${bam%.*}.accepted_hits.bedGraph

bedGraphToBigWig ${bam%.*}.accepted_hits.bedGraph Hg38_shRNA_hybrid.sizes ${bam%.*}.bw

rm -f ${bam%.*}.accepted_hits.bedGraph ${bam%.*}.accepted_hits.bed
