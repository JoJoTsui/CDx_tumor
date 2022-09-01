#!/usr/bin/env python
#coding=utf-8
import pandas as pd
import sys
import os
import re

dfiles = sys.argv[1:]
outfile = dfiles.pop(-1)
maxcoln = 0

writer = pd.ExcelWriter(outfile)
for dfile in dfiles:
	size = os.path.getsize(dfile)
	basen = os.path.basename(dfile)
	basen = re.sub(".txt","",basen)
	if size == 0:
		df = pd.DataFrame()
	else:
		cmd = "awk -F'\t' 'BEGIN{max=0}{if(NF>max){max=NF}}END{print max}' " + dfile
		maxcoln = os.popen(cmd).readline().strip()
		df = pd.read_csv(dfile,sep="\t",names=range(1,int(maxcoln)+1),low_memory=False)
		#df = pd.read_csv(dfile,sep="\t",names=range(1,100),low_memory=False)
	df.to_excel(writer,sheet_name=basen[0:80],index=False,header=False)
writer.save()
writer.close()
