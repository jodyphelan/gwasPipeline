#! /usr/bin/python
import sys
from collections import defaultdict
import subprocess
import json
import os.path

script_dir = os.path.dirname(os.path.realpath(__file__))


bfiles = sys.argv[1].split(",")
prefix = sys.argv[2]
files = [x.split("/")[-1] for x in bfiles]

if os.path.isdir("temp_merging"):
	print "ERROR: Temp directory already present, please clean up"
	quit()
else:
	subprocess.call("mkdir temp_merging",shell=True)

var_dict = defaultdict(int)
pos_dict = {}
set_pos = set()
bad_set = set()

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
				else:
					bad_set.add(rid)


with open("exclude.txt","w") as o:
	for rid in var_dict:
		if rid in bad_set:
			o.write(rid+"\n")

for i in range(len(files)):
	print "Removing ambiguous loci from %s" % (files[i])
	plink_cmd = "plink --bfile %s --exclude exclude.txt --make-bed --out temp_merging/%s > log 2> err" % (bfiles[i],files[i])
	subprocess.call(plink_cmd,shell=True)

temp_name = None
temp_file = files
i=0
while len(temp_file)>0:
	print "Performing round %s of merging" % i
	if temp_name==None:
		f1 = temp_file.pop(0)
		temp_name = "temp_"+str(i)
	else:
		f1 = temp_name
		temp_name = "temp_"+str(i)
	f2 = temp_file.pop(0)
	if len(temp_file)==0:
		temp_name = prefix 

	plink_cmd = "plink --bfile temp_merging/%s --bmerge temp_merging/%s --make-bed --out temp_merging/%s --allow-no-sex > log 2>err" % (f1,f2,temp_name)
	print plink_cmd
	subprocess.call(plink_cmd,shell=True)
	i += 1

subprocess.call("mv temp_merging/%s.fam preimpute/; mv temp_merging/%s.bim preimpute/; mv temp_merging/%s.bed preimpute/" % (prefix,prefix,prefix),shell=True)
subprocess.call("rm -r temp_merging exclude.txt",shell=True)

if "--pca" in sys.argv:
	print "Performing PCA"
	
	idxfile = "plots/"+prefix+".index.txt"
	pngfile = "plots/"+prefix+".pca.png"
	pca_prefix = "plots/"+prefix
	eigenval = pca_prefix+".eigenval"
	eigenvec = pca_prefix+".eigenvec"
	with open(idxfile,"w") as o:
		for f in bfiles:
			for l in open(f+".fam"):
				arr = l.rstrip().split()
				o.write("%s\t%s\n" % (arr[1],f))
	subprocess.call("plink --bfile preimpute/%s --pca --out %s >log 2>err" % (prefix,pca_prefix),shell=True)
	subprocess.call("Rscript %s/plot_gemma_pca.r %s %s %s %s" % (script_dir,eigenvec,eigenval,idxfile,pngfile),shell=True)


open("merge_log.json","w").write(json.dumps(var_dict))
