#!/bin/bash

function usage() {
	echo "Usage: SurVirus Biowulf wrapper for hg38+shRNA"
	echo "bash SurVirus.sh <outdir> <R1.fastq.gz> <R2.fastq.gz> <threads>"
	echo "<outdir> is deleted if it already exists"
}

if [ $# -lt 4 ];then
	usage && exit 1
fi

module load python/2.7
module load bwa/0.7.17
module load samtools/1.10

R1=$2
R2=$3
outdir=$1
threads=$4
extraparameters=$5

if [ -d $outdir ]; then rm -rf $outdir;fi
mkdir -p $outdir
python /data/kopardevn/SandBox/SurVirus/surveyor.py \
	$R1,$R2 \
	$outdir \
	/data/CCBR/projects/ccbr1060/resources/SurVirus/host/hg38.fa \
	/data/CCBR/projects/ccbr1060/resources/SurVirus/shRNA/shRNA.fa \
	/data/CCBR/projects/ccbr1060/resources/SurVirus/host_shRNA/hg38_shRNA.fa \
	--threads $threads \
	--bwa `which bwa` \
	--samtools `which samtools` \
	--fq \
	$5
