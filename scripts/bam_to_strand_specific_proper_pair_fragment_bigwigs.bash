sample=$1
bam="./${sample}/STAR1p/${sample}_p1.Aligned.sortedByCoord.out.bam"
sh bam_split_by_strand.bash $bam chr_shRNA ${sample}.fwd.bam ${sample}.rev.bam
for dirn in "fwd" "rev";do
samtools sort -n -o ${sample}.${dirn}.qsorted.bam ${sample}.${dirn}.bam
bedtools bamtobed -i ${sample}.${dirn}.qsorted.bam -bedpe > ${sample}.${dirn}.qsorted.bam.bedpe
awk -F"\t" -v OFS="\t" '{if ($2>$4){s=$4}else{s=$2}; if ($3>$6){e=$3}else{e=$6};print $1,s,e}' ${sample}.${dirn}.qsorted.bam.bedpe | \
  bedItemOverlapCount hg19 -chromSize=Hg38_shRNA_hybrid.sizes stdin |sort -k1,1 -k2,2n > ${sample}.${dirn}.bg
bedGraphToBigWig ${sample}.${dirn}.bg Hg38_shRNA_hybrid.sizes ${sample}.shRNA_proper_pair_fragments.${dirn}.bw
rm -f ${sample}.${dirn}.qsorted.bam.bedpe ${sample}.${dirn}.bg ${sample}.${dirn}.qsorted.bam ${sample}.${dirn}.bam
done
rm -f ${sample}.bam

