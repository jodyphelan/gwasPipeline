import sys
from collections import defaultdict
import subprocess
import json


bfiles = sys.argv[1].split(",")
files = [x.split("/")[-1] for x in bfiles]

var_dict = defaultdict(int)
pos_dict = {}
set_pos = set()

for i in bfiles:
	with open(i+".bim") as f:
		for line in f.readlines():
			rid,map,pos = line.split()[1:4]
			if rid not in set_pos:
				set_pos.add(rid)
				var_dict[rid] += 1
				pos_dict[rid] = pos
			else:
				if pos == pos_dict[rid]:
					var_dict[rid] += 1

with open("extract.txt","w") as o:
	for rid in var_dict:
		if var_dict[rid]==len(files):
			o.write(rid+"\n")

for i in range(len(files)):
	plink_cmd = "plink --bfile %s --extract extract.txt --make-bed --out %s" % (bfiles[i],files[i])
	subprocess.call(plink_cmd,shell=True)

temp_name = None
temp_file = files
i=0
while len(temp_file)>0:
	if temp_name==None:
		f1 = temp_file.pop(0)
		temp_name = "temp_"+str(i)
	else:
		f1 = temp_name
		temp_name = "temp_"+str(i)
	f2 = temp_file.pop(0)
	if len(temp_file)==0:
		temp_name = "merged"
	plink_cmd = "plink --bfile %s --bmerge %s --make-bed --out %s" % (f1,f2,temp_name)
	subprocess.call(plink_cmd,shell=True)
	i += 1

open("merge_log.json","w").write(json.dumps(var_dict))
