import pandas as pd
import sys
import pybedtools
import uuid
import os

# arguments
# 1. INPUT sample level junctions file
# 2. INPUT group level integration windows bed file ... annotated
# 3. INPUT shRNA integration windows bed file ... annotated
# 4. OUTPUT output junctions file ... annotated
# 5. OUTPUT detailed counts file

random_str=str(uuid.uuid4())

pd.set_option('display.max_columns', None)

# load junctions file
d=pd.read_csv(sys.argv[1],sep="\t",header=0)

bigdict=dict()

# load shRNA integration sites
shRNA_integration_sites=pybedtools.BedTool(sys.argv[3])

# load host integration sites
host_integration_sites=pybedtools.BedTool(sys.argv[2])
tmp=open(random_str+".bed",'w')
iw_annot_strands=dict()
iw_oldname=dict()
for i in host_integration_sites:
    oldname=i.name
    bigdict[oldname+"##+"]=dict()
    bigdict[oldname+"##-"]=dict()
    iwid=oldname.split("|")[0]
    iw_oldname[iwid]=oldname
    if len(oldname.split("#"))==2:
        if oldname.split("#")[1]=="NA":
            annot_strand="unknown"
        else:
            annot_strand=oldname[-1]
    else:
        annot_strand="undetermined"
    iw_annot_strands[iwid]=annot_strand
    newname1="{0}:{1}-{2}|{3}".format(i.chrom,i.start,i.end,"-")
    newname2="{0}:{1}-{2}|{3}".format(i.chrom,i.start,i.end,"+")
    tmp.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(i.chrom,i.start,i.end,newname1,"-"))
    tmp.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(i.chrom,i.start,i.end,newname2,"+"))
tmp.close()
host_integration_sites=pybedtools.BedTool(random_str+".bed")


for i in shRNA_integration_sites:
    for k in bigdict.keys():
        bigdict[k][i.name]=0


def create_bedtool(chrom,start,end,strand):
    if strand=="+":
        st="-"
    if strand=="-":
        st="+"
    bt=pybedtools.BedTool("{0}\t{1}\t{2}\t.\t.\t{3}".format(chrom,start,end,st),from_string=True)
    return bt

def find_overlapping(bt1,bt2,strand=False):
    for c in bt1.closest(bt2,d=True,s=strand):
        if int(c[-1])==0:
            annot=c[9]
            iwid=annot.split("#")[0]
            return iwid
        else:
            return "no_overlap"

#d['sample_name']=sys.argv[2]
#d['host_integration_window']='unknown'
#d['shRNA_integration_window']='unknown'
#>>> d.columns
#Index(['chr_donorA', 'brkpt_donorA', 'strand_donorA', 'chr_acceptorB',
#	       'brkpt_acceptorB', 'strand_acceptorB', 'junction_type',
#	              'repeat_left_lenA', 'repeat_right_lenB', 'read_name', 'start_alnA',
#		             'cigar_alnA', 'start_alnB', 'cigar_alnB', 'num_chim_aln',
#			            'max_poss_aln_score', 'non_chim_aln_score', 'this_chim_aln_score',
#				           'bestall_chim_aln_score', 'PEmerged_bool', 'readgrp', 'sampleName',
#					          'host_integration_window', 'shRNA_integration_window'],
#d=d[(d['chr_donorA']=="chr_shRNA" or d['chr_acceptorB']=="chr_shRNA") and d['chr_donorA']!=d['chr_acceptorB']]

# filter out non-chimeric reads and chimeric reads that do not involve shRNA
d=d[((d['chr_donorA']=="chr_shRNA") | (d['chr_acceptorB']=="chr_shRNA")) & (d['chr_donorA']!=d['chr_acceptorB'])]
d['host_integration_window']='unknown'
d['host_integration_window_annotation']='unknown'
d['shRNA_integration_window']='unknown'
d['is_OST']='unknown'

def is_OST(alignment_strand,annot_strand):
    if alignment_strand=="+" and annot_strand=="-":
        return "yes"
    elif alignment_strand=="-" and annot_strand=="+":
        return "yes"
    elif alignment_strand=="-" and annot_strand=="-":
        return "no"
    elif alignment_strand=="+" and annot_strand=="+":
        return "no"
    else:
        return "undetermined"



for index, row in d.iterrows():
    if row['chr_donorA']!="chr_shRNA":
        bt_host=create_bedtool(row['chr_donorA'],row['brkpt_donorA'],row['brkpt_donorA'],row['strand_donorA'])
        bt_shRNA=create_bedtool(row['chr_acceptorB'],row['brkpt_acceptorB'],row['brkpt_acceptorB'],row['strand_acceptorB'])
    else:
        bt_shRNA=create_bedtool(row['chr_donorA'],row['brkpt_donorA'],row['brkpt_donorA'],row['strand_donorA'])
        bt_host=create_bedtool(row['chr_acceptorB'],row['brkpt_acceptorB'],row['brkpt_acceptorB'],row['strand_acceptorB'])
    closest_host_iw=find_overlapping(bt_host,host_integration_sites,True)
    closest_shRNA_iw=find_overlapping(bt_shRNA,shRNA_integration_sites,True)
    if closest_host_iw=="no_overlap":
        iw_alignment_strand="unknown"
        iw_annot_strand="unknown"
    else:
        iwid,iw_alignment_strand=closest_host_iw.split("|")
        iw_annot_strand=iw_annot_strands[iwid]
        oldname=iw_oldname[iwid]
        d.loc[index,'host_integration_window_annotation']=oldname
        if closest_shRNA_iw!="no_overlap":
            bigdict[oldname+"##"+iw_alignment_strand][closest_shRNA_iw]+=1
    # print(iwid,iw_annot_strand,iw_alignment_strand)
    d.loc[index,'host_integration_window']=closest_host_iw
    d.loc[index,'shRNA_integration_window']=closest_shRNA_iw
    d.loc[index,'is_OST']=is_OST(iw_alignment_strand,iw_annot_strand)

        # print(d[[index]])
        # exit()

# open output file for writing
bigdictoutputfile=open(sys.argv[5],'w')

bigdictoutputfile.write("host_integration_window#annotation##alignment_strand")
for j in shRNA_integration_sites:
    bigdictoutputfile.write("\t"+j.name)
bigdictoutputfile.write("\n")

for i in bigdict.keys():
    bigdictoutputfile.write(i)
    for j in shRNA_integration_sites:
        bigdictoutputfile.write("\t"+str(bigdict[i][j.name]))
    bigdictoutputfile.write("\n")
bigdictoutputfile.close()



d.to_csv(sys.argv[4],sep="\t",header=True,index=False)
os.remove(random_str+".bed")
# print(random_str+".bed")