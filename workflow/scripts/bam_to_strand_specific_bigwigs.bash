#!/bin/bash
module load python

ARGPARSE_DESCRIPTION="split bam into strand specific bigwigs"      # this is optional
source <(curl -s https://raw.githubusercontent.com/CCBR/Tools/master/scripts/argparse.bash) || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--sample',required=True, help='sample name')
parser.add_argument('--bam',required=True, help='input BAM file')
parser.add_argument('--outdir',required=False, help='output DIR')
parser.add_argument('--region',required=False, help='region of BAM')
EOF

module load samtools
module load bedtools
module load ucsc
module load parallel
set -e -x -o pipefail

#set outdir
if [ ! $OUTDIR ];then
OUTDIR=$(pwd)
fi

#set tmpdir
TMPDIR=/lscratch/$SLURM_JOBID
if [ ! -d /lscratch/$SLURM_JOBID ];then
TMPDIR=/dev/shm
fi

#out bigwigs
if [ $REGION ];then
FWD_BW="${OUTDIR}/${SAMPLE}.${REGION}.fwd.bw"
REV_BW="${OUTDIR}/${SAMPLE}.${REGION}.rev.bw"
else
FWD_BW="${OUTDIR}/${SAMPLE}.fwd.bw"
REV_BW="${OUTDIR}/${SAMPLE}.rev.bw"
fi

threads=4

samtools view -H $BAM | grep "^@SQ" | cut -f2,3 | sed "s/SN://g" | sed "s/LN://g" > ${TMPDIR}/${SAMPLE}.genome
echo "samtools view -@ $threads -b -f 128 -F 16 $BAM $REGION | bedtools genomecov -bg -split -ibam stdin -g ${TMPDIR}/${SAMPLE}.genome > ${TMPDIR}/${SAMPLE}.fwd1.bg && bedSort ${TMPDIR}/${SAMPLE}.fwd1.bg ${TMPDIR}/${SAMPLE}.fwd1.bg" > ${TMPDIR}/${SAMPLE}.dobg
echo "samtools view -@ $threads -b -f 80 $BAM $REGION | bedtools genomecov -bg -split -ibam stdin -g ${TMPDIR}/${SAMPLE}.genome > ${TMPDIR}/${SAMPLE}.fwd2.bg && bedSort ${TMPDIR}/${SAMPLE}.fwd2.bg ${TMPDIR}/${SAMPLE}.fwd2.bg" >> ${TMPDIR}/${SAMPLE}.dobg
echo "samtools view -@ $threads -b -f 144 $BAM $REGION | bedtools genomecov -bg -split -ibam stdin -g ${TMPDIR}/${SAMPLE}.genome > ${TMPDIR}/${SAMPLE}.rev1.bg && bedSort ${TMPDIR}/${SAMPLE}.rev1.bg ${TMPDIR}/${SAMPLE}.rev1.bg" >> ${TMPDIR}/${SAMPLE}.dobg
echo "samtools view -@ $threads -b -f 64 -F 16 $BAM $REGION | bedtools genomecov -bg -split -ibam stdin -g ${TMPDIR}/${SAMPLE}.genome > ${TMPDIR}/${SAMPLE}.rev2.bg && bedSort ${TMPDIR}/${SAMPLE}.rev2.bg ${TMPDIR}/${SAMPLE}.rev2.bg" >> ${TMPDIR}/${SAMPLE}.dobg
parallel -j 4 < ${TMPDIR}/${SAMPLE}.dobg
bedtools unionbedg -i ${TMPDIR}/${SAMPLE}.fwd1.bg ${TMPDIR}/${SAMPLE}.fwd2.bg|awk -F"\t" -v OFS="\t" '{sum=0;for (i=4;i<=NF;i+=1) {sum+=$i};print $1,$2,$3,sum}' - > ${TMPDIR}/${SAMPLE}.fwd.bg
bedtools unionbedg -i ${TMPDIR}/${SAMPLE}.rev1.bg ${TMPDIR}/${SAMPLE}.rev2.bg|awk -F"\t" -v OFS="\t" '{sum=0;for (i=4;i<=NF;i+=1) {sum+=$i};print $1,$2,$3,sum}' - > ${TMPDIR}/${SAMPLE}.rev.bg 
bedSort ${TMPDIR}/${SAMPLE}.fwd.bg ${TMPDIR}/${SAMPLE}.fwd.bg 
bedGraphToBigWig ${TMPDIR}/${SAMPLE}.fwd.bg ${TMPDIR}/${SAMPLE}.genome $FWD_BW
bedSort ${TMPDIR}/${SAMPLE}.rev.bg ${TMPDIR}/${SAMPLE}.rev.bg 
bedGraphToBigWig ${TMPDIR}/${SAMPLE}.rev.bg ${TMPDIR}/${SAMPLE}.genome $REV_BW
rm -f ${TMPDIR}/${SAMPLE}.fwd*.bg ${TMPDIR}/${SAMPLE}.rev*.bg ${TMPDIR}/${SAMPLE}.dobg ${TMPDIR}/${SAMPLE}.genome
