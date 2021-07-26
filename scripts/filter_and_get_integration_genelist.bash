if [ "$#" == "0" ];then
threshold_nreads=50
else
threshold_nreads=$1
fi
for f in `find . -name "*_p1.Chimeric.out.junction"`;do 
  sample=$(basename $f|awk -F"_p1" '{print $1}')
  grep chr_shRNA $f | awk -F"\t" '{if ($1!=$4){read[$10]++;if (read[$10]==1){print}}}' > ${f}.shRNA
  awk -F"\t" -v OFS="\t" '{if ($1=="chr_shRNA"){print $4,$5}else{print $1,$2}}' ${f}.shRNA | \
    sort | uniq -c | awk -v OFS="\t" -v s=$sample -v t=$threshold_nreads '{if ($1>t){print $1,$2,$3,s}}'
done | sort -k1,1 -k2,2n -k3,3 > filtered.chimeric.junctions
cat filtered.chimeric.junctions |awk '{array[$2"_"$3]+=$1}END{for (i in array) print i"\t"array[i]}'|sed "s/_/\t/g" > filtered.chimeric.junctions.aggregate_across_samples
python create_integration_windows.py | grep SiteList | awk -v OFS="\t" '{print $4,$6-100,$6+99}' | sed "s/,//g" | sort -k1,1 -k2,2n > integration_windows.bed
bedtools merge -i integration_windows.bed > integration_windows.bed.tmp && mv integration_windows.bed.tmp integration_windows.bed
bedtools intersect -loj -wa -wb -a integration_windows.bed -b /data/CCBR/projects/ccbr1060/resources/hg38/genes.bed6 | awk -F"\t" '{seen[$1":"$2"-"$3]++;if (seen[$1":"$2"-"$3]==1){print}}' > integration_windows.bed.annotated.txt
awk -F"\t" -v OFS="\t" '{print $1":"$2"-"$3,$7"|"$NF}' integration_windows.bed.annotated.txt > integration_windows.annotation_lookup.txt
awk '{print $1}' integration_windows.annotation_lookup.txt | sed "s/:/\t/g" | sed "s/-/\t/g"  > integration_windows.annotation_lookup.txt.a
awk '{print $2}' integration_windows.annotation_lookup.txt >  integration_windows.annotation_lookup.txt.b
paste  integration_windows.annotation_lookup.txt.a  integration_windows.annotation_lookup.txt.b | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3"|"$4,".","."}' > integration_windows.annotated.bed
rm -f  integration_windows.annotation_lookup.txt.a  integration_windows.annotation_lookup.txt.b
bedtools closest -d -a integration_windows.bed -b lab_verified_sites.annotated.bed |awk -F"\t" -v OFS="\t" '{print $1":"$2"-"$3,$7"|"$10}' |awk '{seen[$1]++;if (seen[$1]==1) {print}}' > integration_windows.nearest_lab_verified_site_lookup.txt
cut -f7 integration_windows.bed.annotated.txt | sort | uniq

