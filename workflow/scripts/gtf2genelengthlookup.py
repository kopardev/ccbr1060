import gtfparse
import sys
debug=0

if len(sys.argv)!=2:
	print("1 argument required! Need path to GTF file!")
	exit()

df = gtfparse.read_gtf(sys.argv[1])
bigdict = dict()
if debug==1: print("list of genes=",list(df.gene_id.unique()))
for gene in list(df.gene_id.unique()):
	bigdict[gene]=list()
	x=df[(df["gene_id"]==gene) & (df["feature"]!="gene")]
	if debug==1: print("gene=",gene)
	if debug==1: print("list of transcripts=",list(x.transcript_id.unique()))
	for transcript in list(x.transcript_id.unique()):
		y=x[(x["transcript_id"]==transcript) & (x["feature"]=="exon")]
		l=0
		for i,exonrow in y.iterrows():
			l+=(exonrow["end"]-(exonrow["start"]-1))
		if debug==1: print("transcript=",transcript)
		if debug==1: print("l=",l)
		bigdict[gene].append(l)
	if debug==1: print(bigdict)

for gene in list(df.gene_id.unique()):
	x=df[(df["gene_id"]==gene) & (df["feature"]=="gene")]
	gene_name=x.gene_name.unique()[0]
	maxlen=max(bigdict[gene])
	r=[gene,gene_name,str(maxlen)]
	print("\t".join(r))