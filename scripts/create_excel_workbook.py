import pandas as pd
import sys
import os

# create Excel workbook
# 1. OUTPUT excel file name
# 2. INPUT comma separated list of files that need to be loaded in as sheets in the excel file
# 3. INPUT comma separated list of tab names ... one for each file

writer = pd.ExcelWriter(sys.argv[1], engine='xlsxwriter')
tabnames = sys.argv[3].split(",")
for i,f in enumerate(sys.argv[2].split(",")):
	sn=tabnames[i]
	df=pd.read_csv(f,sep="\t",header=0)
	df.to_excel(writer, sheet_name=sn, index=False)

writer.close()
