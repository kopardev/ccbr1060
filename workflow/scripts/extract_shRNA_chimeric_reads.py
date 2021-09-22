import pysam
import sys
import pandas as pd


# Arguments
# 1. INPUT: annotated junctions file
# 2. INPUT: STAR BAM file
# 3. OUTPUT: output BAM file name

# load annotated junctions file
junctions=pd.read_csv(sys.argv[1],sep="\t",header=0)

readids=set(junctions['read_name'])

inBAM = pysam.AlignmentFile(sys.argv[2], "rb")
outBAM = pysam.AlignmentFile(sys.argv[3], "wb", template=inBAM)
for read in inBAM.fetch():
    if read.query_name in readids:
        outBAM.write(read)
inBAM.close()
outBAM.close()
