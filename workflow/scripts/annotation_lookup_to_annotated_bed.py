# convert annotation lookup to bed format
# input is TSV lookuptable with 2 columns:
# 	1: region to be annotated: <chr>:<start>-<end>|<score>|<strand> format
# 	2: annotation from genes.bed: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand> format
# output is BED6 fromat with columns:
# 	1: chr ... chromosome of the region to be annotated
# 	2: start ... start coordinate of the region to be annotated
# 	3: end ... end coordinate of the region to be annotated
# 	4: name ... <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2>
# 	5: score ... score of the region to be annotated
# 	6: strand ... strand of the region to be annotated 
import sys
lookuptable=sys.argv[1]
d=open(lookuptable).readlines()
bigdict=dict()
for l in d:
	l=l.strip().split("\t")
	if not l[0] in bigdict:
		bigdict[l[0]]=list()
	bigdict[l[0]].append(l[1])
for k,v in bigdict.items():
	a=k.split(":")
	chrom=a[0]
	b=a[1].split("-")
	start=b[0]
	# minus strand correction
	if len(b)==3 and b[2]=="" :
		b[1]=b[1]+"-"
	c=b[1].split("|")
	end=c[0]
	score=c[1]
	strand=c[2]
	x=list()
	x.append(k)
	x.extend(v)
	name="#".join(x)
	outline=[chrom,start,end,name,score,strand]
	print("\t".join(outline))

