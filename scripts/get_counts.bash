#!/bin/bash
for f in `find . -name "*_p1.Chimeric.out.junction"`;do
sample=$(basename $f|awk -F"_p1" '{print $1}')
less ./${sample}/${sample}_p1.Chimeric.out.junction|grep chr_shRNA|awk -F"\t" '{if ($1!=$4){read[$10]++;if (read[$10]==1){print}}}' > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA
awk -F"\t" '{if ($1!=$4){if ($1=="chr_shRNA"){print $4"|"$5"|"$6"|"$10}else{print $1"|"$2"|"$3"|"$10}}}' ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index
awk -F"|" -v OFS="\t" '{print $1,$2,$2,$4,".",$3}' ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index.bed
bedtools intersect -wa -wb -a integration_windows.bed -b ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index.bed |awk '{print $4"|"$5"|"$NF"|"$7}' > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index.filtered_to_integration_windows
paste ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.tsv
cat ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index.filtered_to_integration_windows ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.tsv |awk '{seen[$1]++;if (seen[$1]==2){print}}' > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.tsv
bedtools intersect -wa -wb -b integration_windows.annotated.bed -a ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.index.bed |awk -F"\t" -v OFS="\t" '{print $1"|"$2"|"$6"|"$4,$10}' > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA_to_integration_window.lookup.txt
paste ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA_to_integration_window.lookup.txt ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.tsv |cut -f1-2,4- > ./${sample}/${sample}_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.annotated.tsv
python get_detailed_chimera_counts.py ${sample}/${sample}_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.annotated.tsv ${sample}/${sample}_shRNA_integration_window_counts.tsv ${sample}/${sample}_shRNA_integration_window_coordinates.tsv
grep True ${sample}/${sample}_shRNA_integration_window_coordinates.tsv|awk '{print $NF}' > ${sample}/${sample}_shRNA.OST.readids
grep False ${sample}/${sample}_shRNA_integration_window_coordinates.tsv|awk '{print $NF}' > ${sample}/${sample}_shRNA.not_OST.readids
done

for i in {1..8};do awk -F"\t" -v t="test${i}" '{print t,$0}' test${i}/test${i}_shRNA_integration_window_counts.tsv;done|grep -v integration|sort -k2,2|grep True|awk '{seen[$2] = seen[$2]","$1}END{for (i in seen){print i,seen[i]}}' |sed "s/ ,/\t/g"|sed "s/|/\t/g"|awk -F"\t" -v OFS="\t" '{print $1,$NF}' > integration_window.samples_w_OST.lookup.txt
