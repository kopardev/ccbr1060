#!/bin/bash
module load bedtools
inBed=$1
genesBed6=$2
lookuptable=$3
ncol=$(awk '{print NF}' $inBed|sort|uniq)
if [ "$ncol" == "3" ];then
	awk -v OFS="\t" '{print $1,$2,$3,".",".","."}' $inBed > ${inBed}.bed6
	to_annotate="${inBed}.bed6"
fi
if [ "$ncol" == "6" ];then
	to_annotate="$inBed"
fi
bedtools intersect -loj -wa -wb -a $to_annotate -b $genesBed6 | \
 awk -F"\t" -v OFS="\t" '{
   if ($8==-1)
	   {print $1":"$2"-"$3"|"$5"|"$6,"NA"}
   else
	   {print $1":"$2"-"$3"|"$5"|"$6,$7":"$8"-"$9"|"$10"|"$12}
   }' > $lookuptable

# cleanup
if [ -f ${inBed}.bed6 ];then rm -f ${inBed}.bed6 ;fi
