
import pandas
from functools import reduce
import glob

def readdf(filename,sn):
	data=pandas.read_csv(filename,sep=" ",header=None,names=["integration_window",sn])
	return data

countsfiles=glob.glob("*.counts")
countsfiles.sort()
list_of_dfs=[]

anno=pandas.read_csv("integration_windows.annotation_lookup.txt",sep="\t",header=None,names=["integration_window","gene"])
list_of_dfs.append(anno)

anno2=pandas.read_csv("integration_windows.nearest_lab_verified_site_lookup.txt",sep="\t",header=None,names=["integration_window","nearest_verified_site"])
list_of_dfs.append(anno2)


for f in countsfiles:
	samplename=f.split(".")[0]
	list_of_dfs.append(readdf(f,samplename))

mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="integration_window"),list_of_dfs)
mergeddf.fillna(0,inplace=True)
mergeddf.drop_duplicates(inplace=True)
mergeddf=mergeddf.infer_objects()
mergeddf.to_csv("integration_windows.counts_matrix.tsv",sep="\t",index=False)
