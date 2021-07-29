bam=$1
region=$2
fwdbam=$3
revbam=$4
samtools view -b -f 128 -F 16 $bam $region > ${bam%.*}.fwd1.bam
samtools view -b -f 80 $bam $region > ${bam%.*}.fwd2.bam
samtools merge -f -o $fwdbam ${bam%.*}.fwd1.bam ${bam%.*}.fwd2.bam
samtools view -b -f 144 -F 16 $bam $region > ${bam%.*}.rev1.bam
samtools view -b -f 64 $bam $region > ${bam%.*}.rev2.bam
samtools merge -f -o $revbam ${bam%.*}.rev1.bam ${bam%.*}.rev2.bam
rm -rf ${bam%.*}.fwd?.bam ${bam%.*}.rev?.bam
