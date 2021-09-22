
def create_bedtool(chrom,start,end,strand):
    if strand=="+":
        st="-"
    if strand=="-":
        st="+"
    bt=pybedtools.BedTool("{0}\t{1}\t{2}\t.\t.\t{3}".format(chrom,start,end,st),from_string=True)
    return bt

def find_overlapping(bt1,bt2):
    # print(bt1,bt2)
    for c in bt1.closest(bt2,d=True,s=True):
        # print(c)
        if int(c[-1])==0:
            return c[9]
    else:
        return "no_overlap"


import pybedtools
import sys,os

bigDict=dict()


# load shRNA integration sites
shRNA_integration_sites=pybedtools.BedTool("integration_windows.shRNA.bed")

# load host integration sites
host_integration_sites=pybedtools.BedTool("integration_windows.annotated.bed")
tmp=open("integration_windows.annotated.tmp.bed",'w')
for i in host_integration_sites:
    newname1=i.name+"|-"
    newname2=i.name+"|+"
    tmp.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(i.chrom,i.start,i.end,newname1,"-"))
    tmp.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(i.chrom,i.start,i.end,newname2,"+"))
tmp.close()
host_integration_sites=pybedtools.BedTool("integration_windows.annotated.tmp.bed")

for i in host_integration_sites:
    bigDict[i.name]=dict()
    for j in shRNA_integration_sites:
        bigDict[i.name][j.name]=0


junctions_data=open(sys.argv[1]).readlines()
for jd in junctions_data:
    jd=jd.strip().split("\t")
    # is shRNA chimera?
    is_shRNA_chimera=False
    if jd[0] == "chr_shRNA" or jd[3] == "chr_shRNA":
        is_shRNA_chimera=True
    if is_shRNA_chimera==False:
        continue
    if jd[0]=="chr_shRNA":
        bt1=create_bedtool(jd[0],jd[1],jd[1],jd[2])
        bt2=create_bedtool(jd[3],jd[4],jd[4],jd[5])
    else:
        bt2=create_bedtool(jd[0],jd[1],jd[1],jd[2])
        bt1=create_bedtool(jd[3],jd[4],jd[4],jd[5])
    # bt1 is shRNA and bt2 is host
    # print(bt1)
    # print(bt2)
    bt1_closest=find_overlapping(bt1,shRNA_integration_sites)
    # print(bt1_closest)
    if bt1_closest=="no_overlap":
        continue
    bt2_closest=find_overlapping(bt2,host_integration_sites)
    # print(bt2_closest)
    if bt2_closest=="no_overlap":
        continue
    # print(jd)
    # print(bt1)
    # print(bt2)
    bigDict[bt2_closest][bt1_closest]+=1

print("",end="\t")
for j in shRNA_integration_sites:
    print(j.name,end="\t")
print("")

for i in host_integration_sites:
    print(i.name,end="\t")
    for j in shRNA_integration_sites:
        print(str(bigDict[i.name][j.name]),end="\t")
    print("")
os.remove("integration_windows.annotated.tmp.bed")

