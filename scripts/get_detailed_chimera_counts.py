import sys
infile=sys.argv[1] # eg. test5/test5_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.annotated.tsv
# % head test5/test5_p1.Chimeric.out.junction.shRNA.filtered_to_integration_windows.annotated.tsv
# chr5|73643317|+|NB501223:102:HMC2KBGXF:1:11212:14310:16860	chr5:73643062-73643416|ENSG00000214944.10|ARHGEF28|+	chr5	73643317	+	chr_shRNA	8407	-	0	0	1	NB501223:102:HMC2KBGXF:1:11212:14310:16860	73643270	47M104S	8107	149M47p104M47S	1	300	256	296	296	0
# chr5|113546953|+|NB501223:102:HMC2KBGXF:1:11212:23853:17770	chr5:113546853-113547052|ENSG00000047188.16|YTHDC2|+	chr5	113546953	+	chr_shRNA	4554	-	2	0	1	NB501223:102:HMC2KBGXF:1:11212:23853:17770	113546719	150M-12p96M48S	4506	48M96S	1	294	245	291	291	0
# chr5|73642754|+|NB501223:102:HMC2KBGXF:1:11301:7615:2612	chr5:73642533-73642853|ENSG00000214944.10|ARHGEF28|+	chr5	73642754	+	chr_shRNA	4554	-	2	0	1	NB501223:102:HMC2KBGXF:1:11301:7615:2612	73642482	151M75p46M104S	4450	104M46S	1	301	196	297	297	0
out_counts_file=sys.argv[2]
out_shRNA_coordinates_file=sys.argv[3]

def readtsvline(line):
    l=line.strip().split("\t")
    integration_window=l[1]
    gene_strand=integration_window.split("|")[-1]
    if l[2]=="chr_shRNA":
        shRNA_strand=l[4]
        shRNA_coordinate=l[3]
        genomic_strand=l[7]
        genomic_coordinate=l[6]
    else:
        shRNA_strand=l[7]
        shRNA_coordinate=l[6]
        genomic_strand=l[4]
        genomic_coordinate=l[3]
    return integration_window,gene_strand,genomic_strand,genomic_coordinate,shRNA_strand,shRNA_coordinate

def get_sums(counts):
    gene_strand=counts['gene_strand']
    same_strand_sum=0
    opposite_strand_sum=0
    sumall=counts['++']+counts['+-']+counts['-+']+counts['--']
    if gene_strand=="+":
        same_strand_sum=counts['++']+counts['+-']
        opposite_strand_sum=counts['-+']+counts['--']
    elif gene_strand=="-":
        same_strand_sum=counts['-+']+counts['--']
        opposite_strand_sum=counts['++']+counts['+-']
    sumall=same_strand_sum+opposite_strand_sum
    # if same_strand_sum/sumall>0.2:
    if same_strand_sum>=5:
        opposite_strand_transcription=True
    else:
        opposite_strand_transcription=False
    # print(gene_strand,same_strand_sum,opposite_strand_sum,sumall,opposite_strand_transcription)
    return same_strand_sum,opposite_strand_sum,sumall,opposite_strand_transcription


counts=dict()
outcoord=open(out_shRNA_coordinates_file,'w')
outcoord.write("\t".join(["integration_window","shRNA_strand","shRNA_coordinate"])+'\n')
for line in open(infile).readlines():
    iw,gene_strand,gs,gc,ss,sc=readtsvline(line)
    outcoord.write("{0}\t{1}\t{2}\n".format(iw,ss,sc))
    if not iw in counts:
        counts[iw]=dict()
        counts[iw]["++"]=0
        counts[iw]["+-"]=0
        counts[iw]["-+"]=0
        counts[iw]["--"]=0
        counts[iw]['gene_strand']=gene_strand
    counts[iw][gs+ss]+=1
outcoord.close()

outcounts=open(out_counts_file,'w')
outcounts.write("\t".join(["integration_window","gene_strand","counts_++","counts_+-","counts_-+","counts_--","same_strand_counts","opposite_strand_counts","total_counts","opposite_strand_transcription"])+'\n')
for iw in counts.keys():
    # print(counts[iw])
    sss,oss,s,ost=get_sums(counts[iw])
    outcounts.write('\t'.join([iw,counts[iw]['gene_strand'],str(counts[iw]['++']),str(counts[iw]['+-']),str(counts[iw]['-+']),str(counts[iw]['--']),str(sss),str(oss),str(s),str(ost)])+'\n')
outcounts.close()