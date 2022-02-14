# ex. usage in Snakefile:
# python {params.script1} {input.iwnearestlvlookup} {input.nreads} {params.samplemanifest} $counts_files > {output.count_matrix}
# @Inputs:
# @param: iwnearestlvlookup <tsv_file>
# TSV with columns:
# 	1. iw in format: <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
# 	2. lv in format: <chr>:<start>-<end>|<score>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <region to-be annotated>#<annotation1>#<annotation2> 
# 	3. distance in bp. Note of col3 (distance in bp) = "-1", then no lab-verified integration site found on the chromosome of interest
# @param: nreads <tsv_file>
# TSV with columns:
# 	1. samplename or replicate name
# 	2. sum of aligned reads per sample according to STAR
# @param: samplemanifest <tsv_file> with header
# TSV with columns:
# 1. sampleName --> this is actually replicate name
# 2. group --> replicates that belong to the same group ...its name ... This is actually samplename
# 3. dox --> Dox added Y or N
# 4. path_to_R1_fastq
# 5. path_to_R2_fastq
# @params: list of counts files <space-separated list> 1 for each replicate
# each counts file is a TSV with columns:
# 1. iw --> chr:start-end|score|strand ... score = number of integration sites in this window for all replicates of this sample
# 2. count --> count chimeric_shRNA_donor_acceptor_sites in this window for this replicate
# @Outputs:
# @param: count_matrix <tsv_file>
# TSV with columns:
# 1. integration_window
# format: <chr>:<start>-<end>|<score>|<strand> score=number of integration sites in each integration window aggregated accross sample replicates
# 2. integration_window_annotation
# format: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <annotation1>#<annotation2> 
# 3. nearest_lab_verified_integration_site
# format: <chr>:<start>-<end>|<score>|<strand> score=number of reads lab-verified
# 4. nearest_lab_verified_integration_site_annotation
# format: <chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>#<chr>:<start>-<end>|<ensembl_id>_<genename>|<strand>... ie <annotation1>#<annotation2>
# 5. distance_to_nearest_lab_verified_integration_site
# format: distance in bp. Note if = "-1", then no lab-verified integration site found on the chromosome of interest	
# 6. Res1_raw_counts	Res2_raw_counts	Res3_raw_counts
# raw counts for each replicate in tab-delimited columns
# 7. Res1_normalized_counts	Res2_normalized_counts	Res3_normalized_counts
# normalized counts for each replicate in tab-delimited columns	
# 8. t_test_p_value

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
