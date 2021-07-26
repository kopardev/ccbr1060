#!/bin/bash
for f in `find . -name "*_p1.Chimeric.out.junction"`;do
sample=$(basename $f|awk -F"_p1" '{print $1}')
grep $sample filtered.chimeric.junctions|awk -v OFS="\t" '{print $2,$3,$3,".",$1,"."}'|sort -k1,1 -k2,2n > ${sample}.filtered.chimeric.junctions.bed
bedtools intersect -wa -wb -a ${sample}.filtered.chimeric.junctions.bed -b integration_windows.bed |awk -F"\t" '{a[$(NF-2)":"$(NF-1)"-"$NF]+=$5}END{for (i in a) { print i,a[i] }}' > ${sample}.integration_windows.counts
done
