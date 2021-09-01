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

