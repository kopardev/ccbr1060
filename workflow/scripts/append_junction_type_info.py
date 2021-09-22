import pandas as pd
import sys
pd.set_option('display.max_columns', None)

# adds columns regarding junction info status to aggregate OST calls
# arguments:
# 1. INPUT : <annotated junctions files from annotate_junctions.py>
# 2. INPUT: <aggregate OST calls file>
# 3. OUTPUT: <aggregate OST calls file with junction type columns added>
#
# create <aggregate OST calls file>
#
# create header first
#
# echo -ne "nreads\t" > <aggregate OST calls file>
# head -n1 <annotated junctions files from annotate_junctions.py> |cut -f22- >> <aggregate OST calls file>
#
# now add data by searching for is_OST == "yes"
#
# grep yes <annotated junctions files from annotate_junctions.py> |cut -f22-|sort|uniq -c |awk -v OFS="\t" '{print $1,$2,$3,$4,$5}'|sort -k1,1nr >> <aggregate OST calls file>

d=pd.read_csv(sys.argv[1],sep="\t",header=0)
l=pd.read_csv(sys.argv[2],sep="\t",header=0)
l['junction_type_encompassing_junction']=0
l['junction_type_other']=0
l['junction_type_GT/AG']=0
l['junction_type_CT/AC']=0
jtcode=dict()
jtcode['junction_type_encompassing_junction']=-1
jtcode['junction_type_other']=0
jtcode['junction_type_GT/AG']=1
jtcode['junction_type_CT/AC']=2
for index,row in l.iterrows():
	d2=d[(d['host_integration_window']==row['host_integration_window']) & (d['shRNA_integration_window']==row['shRNA_integration_window']) & (d['is_OST']=="yes")]
	for junction_type in jtcode.keys():
		#print(d2['junction_type'])
		#print(jtcode[junction_type])
		#exit()
		d3=d2[d2['junction_type']==jtcode[junction_type]]
		l.loc[index,junction_type]=len(d3)

l.to_csv(sys.argv[3],sep="\t",header=True,index=False)
