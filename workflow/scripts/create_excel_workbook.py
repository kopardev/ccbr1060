import pandas as pd
import sys
import os

# create Excel workbook
# 1. OUTPUT excel file name
# 2. INPUT space separated list of files that need to be loaded in as sheets in the excel file
# string before the first dot in the file basename is used as the name of the sheet.

writer = pd.ExcelWriter(sys.argv[1], engine='xlsxwriter')
for f in sys.argv[2:]:
	bn=os.path.basename(f)
	sn=bn.split(".")[0]
	df=pd.read_csv(f,sep="\t",header=0)
	df.to_excel(writer, sheet_name=sn, index=False)

writer.close()
