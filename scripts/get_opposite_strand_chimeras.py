import sys

for line in sys.stdin:
	l=line.strip().split("\t")
	s1=l[1][-1]
	if s1==".":
		continue
	if l[2]=="chr_shRNA":
		c=l[3]
		s2=l[7]
		s3=l[4]
	else:
		c=l[6]
		s2=l[4]
		s3=l[7]
	if s1==s2:
		print("{0}\t{1}\t{2}".format(l[1],c,s3))
