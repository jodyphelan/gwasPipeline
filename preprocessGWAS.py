import sys
from tqdm import tqdm
import subprocess
import gzip

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def gz_file_len(fname):
	with gzip.open(fname,"rb") as f:
		for i,l in enumerate(f):
			pass
	return i+1

def fa2dict(filename):
	fa_dict = {}
	seq_name = ""
	with open(filename) as f:
		for i in tqdm(range(file_len(filename))):
			line = f.readline().rstrip()
			if line[0] == ">":
				seq_name = line[1:].split()[0]
				fa_dict[seq_name] = []
			else:
				fa_dict[seq_name].append(line)
	result = {}
	for seq in fa_dict:
		result[seq] = "".join(fa_dict[seq])
	return result

sample_name = sys.argv[1]
bim_file =  "genotypes/"+sample_name+".bim"
ref_file = sys.argv[2]

def flip_snp():
	print "Loading reference"
	fa_dict = fa2dict(ref_file)
	setAT = set(["A","T"])
	setCG = set(["C","G"])
	setACGT = set(["A","C","G","T"])
	ambFile = open(sample_name+".ambiguous.txt","w")
	revFile = open(sample_name+".reverse.txt","w")
	errFile = open(sample_name+".notFound.txt","w")
	okFile = open(sample_name+".ok.txt","w")
	print "Analysing bim file"
	with open(bim_file) as f:
		for i in tqdm(range(file_len(bim_file))):
			line = f.readline()
			chr,rid,temp,pos,ref,alt = line.rstrip().split()
			if chr not in fa_dict:
				errFile.write(rid+"\n")
				continue
			if (ref in setAT and alt in setAT) or (ref in setCG and alt in setCG):
				ambFile.write(rid+"\n")
				continue
			ref_nuc = fa_dict[chr][int(pos)-1]
			if ref_nuc not in setACGT:
				ambFile.write(rid+"\n")
				continue
			if ref not in setACGT:
				#for removing D and I
				ambFile.write(rid+"\n")
				continue
			if (ref_nuc == ref) or (ref_nuc == alt):
				okFile.write(rid+"\n")
			elif (ref_nuc != ref) and (ref_nuc != ref):
				revFile.write(rid+"\n")
	subprocess.call("cat %s.notFound.txt %s.ambiguous.txt > %s.exclude.txt" % (sample_name,sample_name,sample_name), shell=True)
	plink_cmd = "plink --no-sex --bfile %s --exclude %s.exclude.txt --flip %s.reverse.txt --make-bed --out %s" % ("genotypes/"+sample_name,sample_name,sample_name,sample_name+".flipped")
	subprocess.call(plink_cmd,shell=True)

flip_snp()

