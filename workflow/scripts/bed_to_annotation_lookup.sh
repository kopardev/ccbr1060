#!/bin/bash
# @Description:
# @param1: bed6 file to be annotated (bedfile with only 3 columns then dots will be added to columns 4 through 6 to make this work)
# @param2: bed6 to extract annotations from "genes.bed6" already created for hg38
# @param3: tab-delimited lookup table for @param1
#
# output (@params3) format example:
# column1: region to be annotated: <chr>:<start>-<end>|<score>|<strand> format
# column2: annotation from genes.bed: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand> format
# chr13:95131648-95131847|45|.	chr13:95019835-95301475|ENSG00000125257_ABCC4|-
# chr13:95133345-95133544|52|.	chr13:95019835-95301475|ENSG00000125257_ABCC4|-
# chr13:95149166-95149365|20|.	chr13:95019835-95301475|ENSG00000125257_ABCC4|-
# chr1:71072449-71072648|16|.	chr1:71063291-71081289|ENSG00000132485_ZRANB2|-
# chr1:71076694-71076893|6|.	chr1:71063291-71081289|ENSG00000132485_ZRANB2|-
# 
# genes.bed6 (@params2) format:
# 1. chromosome
# 2. start
# 3. end
# 4. <ensemblID>_<HUGOname>
# 5. "."  for score
# 6. strand
# example lines:
# chr1	11869	14409	ENSG00000223972_DDX11L1	.	+
# chr1	14404	29570	ENSG00000227232_WASH7P	.	-
# chr1	17369	17436	ENSG00000278267_MIR6859-1	.	-
# chr1	29554	31109	ENSG00000243485_MIR1302-2HG	.	+
# chr1	30366	30503	ENSG00000284332_MIR1302-2	.	+

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
