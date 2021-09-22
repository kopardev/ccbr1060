import sys,os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def t_test(x,y,alternative='both-sided'):
        _, double_p = ttest_ind(x,y,equal_var = False)
        if alternative == 'both-sided':
            pval = double_p
        elif alternative == 'greater':
            if np.mean(x) > np.mean(y):
                pval = double_p/2.
            else:
                pval = 1.0 - double_p/2.
        elif alternative == 'less':
            if np.mean(x) < np.mean(y):
                pval = double_p/2.
            else:
                pval = 1.0 - double_p/2.
        return pval


annotated_lvnearest=sys.argv[1]
bigdict=dict()
for x in open(annotated_lvnearest).readlines():
	x=x.strip().split("\t")
	y=x[0].split("#")
	iwid=y[0]
	bigdict[iwid]=dict()
	bigdict[iwid]["lvnearestdist"]=x[2]
	bigdict[iwid]["iwanno"]=y[1]
	if not bigdict[iwid]["lvnearestdist"]=="-1":
		z=x[1].split("#")
		bigdict[iwid]["lvnearest"]=z[0]
		bigdict[iwid]["lvnearestanno"]=z[1]
	else:
		bigdict[iwid]["lvnearest"]="."
		bigdict[iwid]["lvnearestanno"]="."

nreadsfile=sys.argv[2]
nreads=dict()
for x in open(nreadsfile).readlines():
	x=x.strip().split("\t")
	nreads[x[0]]=x[1]

samplesfile=sys.argv[3]
SAMPLESDF = pd.read_csv(samplesfile,sep="\t",header=0,index_col="sampleName")

samples=list()
for x in sys.argv[4:]:
	fn=os.path.basename(x)
	s=fn.split(".counts.")[0]
	samples.append(s)
	for z in open(x).readlines():
		z=z.strip().split("\t")
		bigdict[z[0]][s]=z[1]


header=list()
header.append("integration_window")
header.append("integration_window_annotation")
header.append("nearest_lab_verified_integration_site")
header.append("nearest_lab_verified_integration_site_annotation")
header.append("distance_to_nearest_lab_verified_integration_site")
for s in samples:
	header.append(s+"_raw_counts")
for s in samples:
	header.append(s+"_normalized_counts")
header.append("t_test_p_value")
print("\t".join(header))
for i in bigdict.keys():
	row=list()
	row.append(i)
	row.append(bigdict[i]["iwanno"])
	row.append(bigdict[i]["lvnearest"])
	row.append(bigdict[i]["lvnearestanno"])
	row.append(bigdict[i]["lvnearestdist"])
	normalized_counts=list()
	no_dox=[]
	dox=[]
	for s in samples:
		try:
			row.append(bigdict[i][s])
			nc=float(bigdict[i][s])*1000000/float(nreads[s])
			if SAMPLESDF.loc[s,'dox']==0:
				no_dox.append(nc)
			if SAMPLESDF.loc[s,'dox']==1:
				dox.append(nc)
			normalized_counts.append(str(nc))
		except KeyError:
			row.append('0')
			normalized_counts.append('0')
			if SAMPLESDF.loc[s,'dox']==0:
				no_dox.append(0)
			if SAMPLESDF.loc[s,'dox']==1:
				dox.append(0)
	row.extend(normalized_counts)
	pval=str(t_test(no_dox,dox,'less'))
	row.append(pval)
	print("\t".join(row))
