import pandas as pd
import sys

pd.set_option('display.max_columns', None)

d=pd.read_csv(sys.argv[1],sep="\t",header=0)
sn=sys.argv[2]

donor_is_shRNA = d['chr_donorA']=="chr_shRNA"
acceptor_is_shRNA = d['chr_acceptorB']=="chr_shRNA"
either_is_shRNA = donor_is_shRNA | acceptor_is_shRNA
both_is_shRNA = donor_is_shRNA & acceptor_is_shRNA
only_either_is_shRNA = either_is_shRNA & (~both_is_shRNA )

d2 = d[only_either_is_shRNA]

d3 = d2[d2['chr_donorA']=="chr_shRNA"]
d4 = d2[d2['chr_acceptorB']=="chr_shRNA"]

print("\t".join(["sample","host_integration_window","host_integration_window_annotation","gene_name","shRNA_coordinate","host_integration_window_strand","gene_strand","shRNA_strand","shRNA_integration_window","is_OST"]))

for index, row in d3.iterrows():
	shRNA_coord=str(row['brkpt_donorA'])
	shRNA_strand=row['strand_donorA']
	if shRNA_strand=="+": shRNA_strand="-"
	if shRNA_strand=="-": shRNA_strand="+"
	hiw=row['host_integration_window']
	hiw_annotation=row['host_integration_window_annotation']
	hiw_strand="."
	if hiw!="no_overlap": 
		hiw_strand=hiw[::-1][0]	
	else:
		continue
	genecoords,genename,genestrand,ensemblID,hugoname=".",".",".",".","."
	if hiw_annotation!="unknown":
		a1=hiw_annotation.split("#")[1]
		if a1!="NA":
			genecoords,genename,genestrand=a1.split("|")
			ensemblID,hugoname=genename.split("_")
	shRNA_iw=row['shRNA_integration_window']
	ost=row['is_OST']
	print("\t".join([sn,hiw,hiw_annotation,hugoname,shRNA_coord,hiw_strand,genestrand,shRNA_strand,shRNA_iw,ost]))
for index, row in d4.iterrows():
	shRNA_coord=str(row['brkpt_acceptorB'])
	shRNA_strand=row['strand_acceptorB']
	if shRNA_strand=="+": shRNA_strand="-"
	if shRNA_strand=="-": shRNA_strand="+"
	hiw=row['host_integration_window']
	hiw_strand="."
	if hiw!="no_overlap": 
		hiw_strand=hiw[::-1][0]	
	else:
		continue
	genecoords,genename,genestrand,ensemblID,hugoname=".",".",".",".","."
	if hiw_annotation!="unknown":
		a1=hiw_annotation.split("#")[1]
		if a1!="NA":
			genecoords,genename,genestrand=a1.split("|")
			ensemblID,hugoname=genename.split("_")
	hiw_annotation=row['host_integration_window_annotation']
	shRNA_iw=row['shRNA_integration_window']
	ost=row['is_OST']
	print("\t".join([sn,hiw,hiw_annotation,hugoname,shRNA_coord,hiw_strand,genestrand,shRNA_strand,shRNA_iw,ost]))
