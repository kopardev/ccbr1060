
import pandas
from functools import reduce
import glob

def readdf(filename,sn):
    data=pandas.read_csv(filename,sep="\t",header=0,usecols=[0,8],names=["integration_window_plus",sn])
    data["integration_window"]=list(map(lambda x:x.split("|")[0],data["integration_window_plus"]))
    data=data.drop(columns=["integration_window_plus"])
    data=data.astype({sn: 'int32'})
    # print(sn)
    # print(data)
    return data

list_of_dfs=[]

anno=pandas.read_csv("integration_windows.annotation_lookup.txt",sep="\t",header=None,names=["integration_window","gene"])
list_of_dfs.append(anno)

# print(anno)
anno2=pandas.read_csv("integration_windows.nearest_lab_verified_site_lookup.txt",sep="\t",header=None,names=["integration_window","nearest_verified_site"])
list_of_dfs.append(anno2)

anno3=pandas.read_csv("integration_window.samples_w_OST.lookup.txt",sep="\t",header=None,names=["integration_window","samples_w_OST"])
list_of_dfs.append(anno3)

# print(anno2)
# exit()

countsfiles=glob.glob("**/*_shRNA_integration_window_counts.tsv")
countsfiles.sort()

for f in countsfiles:
    samplename=f.split("/")[0]
    print(samplename,f)
    # exit()
    list_of_dfs.append(readdf(f,samplename))

mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="integration_window"),list_of_dfs)
mergeddf.fillna(0,inplace=True)
mergeddf.drop_duplicates(inplace=True)
mergeddf=mergeddf.infer_objects()
mergeddf.to_csv("integration_windows.counts_matrix.v2.tsv",sep="\t",index=False)
