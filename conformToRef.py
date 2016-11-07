import sys
from tqdm import tqdm
import subprocess
import gzip 
import string

chrom = sys.argv[1]
sample = sys.argv[2]

raw_sample_vcf ="%s.%s.raw.vcf"%(sample,chrom)
raw_sample_bed ="%s.%s.bed"%(sample,chrom)
tabix = "/home/jody/software/bcftools-1.3.1/htslib-1.3.1/tabix"
plink_cmd = "plink --bfile %s.flipped --chr %s --recode vcf --out %s.%s.raw" % (sample,chrom,sample,chrom)
subprocess.call(plink_cmd,shell=True)
idSet = set()
ref_var_dict = {}
pos2id_dict = {}
posSet = set()
print "Processing chr"+chrom+" VCF"
refVCF = "refVCF/chr"+chrom+".1kg.phase3.v5a.vcf.gz"
#with gzip.open(refVCF,"rb") as f:
#	for i,line in enumerate(tqdm(f)):
#		if line[0] == "#": 
#			continue
#		chr,pos,rid,ref,alt = line.rstrip().split()[:5]
#		idSet.add(rid)
#		ref_var_dict[rid] = {"chr":chr,"pos":pos,"ref":ref,"alt":alt}


with open(raw_sample_bed,"w") as o:
	with open(raw_sample_vcf) as f:
		for i,line in enumerate(tqdm(f)):
			if line[0] == "#":
				continue
			chr,pos,rid,ref,alt = line.rstrip().split()[:5]
			o.write("%s\t%s\t%s\n"%(chr,pos,pos))

proc = subprocess.Popen([tabix,refVCF,"-T",raw_sample_bed],stdout=subprocess.PIPE)
for line in proc.stdout:
	chr,pos,rid,ref,alt = line.rstrip().split()[:5]
	idSet.add(rid)
	ref_var_dict[rid] = {"chr":chr,"pos":pos,"ref":ref,"alt":alt}
	pos2id_dict[pos] = rid 

logFh = open("%s.log" % sample,"w")
targetFh = open("%s.%s.preimpute.vcf" % (sample,chrom),"w")
with open("%s.%s.raw.vcf"%(sample,chrom)) as f:
	for i,line in enumerate(tqdm(f)):
		if line[0] == "#":
			targetFh.write(line)
			continue
		arr = line.rstrip().split()
		chr,pos,rid,ref,alt = arr[:5]
		if rid in idSet:
			if (ref == ref_var_dict[rid]["alt"]) and (alt == ref_var_dict[rid]["ref"]):
				logFh.write(rid+"\n")
				arr[3],arr[4] = arr[4],arr[3]	
				for i in range(9,len(arr)):
					arr[i] = arr[i].translate(string.maketrans("01","10"))
				targetFh.write("\t".join(arr)+"\n")	
			elif(alt == ref_var_dict[rid]["alt"]) and (ref == ref_var_dict[rid]["ref"]):
				targetFh.write("\t".join(arr)+"\n")
			else:
				pass
		else:
			pass
				  


