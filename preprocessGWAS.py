import sys
from tqdm import tqdm
import subprocess
import gzip
import os
import argparse
import json

script_dir = os.path.dirname(os.path.realpath(__file__))

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


def flip_snp(ref_file,sample,threads):
	print "Loading reference"
	fa_dict = fa2dict(ref_file)
	bim_file =  "genotypes/"+sample+".bim"

	setAT = set(["A","T"])
	setCG = set(["C","G"])
	setACGT = set(["A","C","G","T"])
	ambFile = open(sample+".ambiguous.txt","w")
	revFile = open(sample+".reverse.txt","w")
	errFile = open(sample+".notFound.txt","w")
	okFile = open(sample+".ok.txt","w")
	log_dict = {"ok":[],"amb":[],"rev":[],"err":[]}
	print "Analysing bim file"
	with open(bim_file) as f:
		for i in tqdm(range(file_len(bim_file))):
			line = f.readline()
			chr,rid,temp,pos,ref,alt = line.rstrip().split()
			if chr not in fa_dict:
				errFile.write(rid+"\n")
				log_dict["err"].append(rid)
				continue
			if (ref in setAT and alt in setAT) or (ref in setCG and alt in setCG):
				ambFile.write(rid+"\n")
				log_dict["amb"].append(rid)
				continue
			ref_nuc = fa_dict[chr][int(pos)-1]
			if ref_nuc not in setACGT:
				ambFile.write(rid+"\n")
				log_dict["amb"].append(rid)
				continue
			if ref not in setACGT:
				#for removing D and I
				ambFile.write(rid+"\n")
				log_dict["amb"].append(rid)
				continue
			if (ref_nuc == ref) or (ref_nuc == alt):
				log_dict["ok"].append(rid)
				okFile.write(rid+"\n")
			elif (ref_nuc != ref) and (ref_nuc != ref):
				revFile.write(rid+"\n")
				log_dict["rev"].append(rid)

	ambFile.close()
	revFile.close()
	errFile.close()
	okFile.close()
	subprocess.call("cat %s.notFound.txt %s.ambiguous.txt > %s.exclude.txt" % (sample,sample,sample), shell=True)
	plink_cmd = "plink --no-sex --bfile %s --exclude %s.exclude.txt --flip %s.reverse.txt --make-bed --out %s" % ("genotypes/"+sample,sample,sample,sample+".flipped")
	subprocess.call(plink_cmd,shell=True)
	open(sample+".flipFilter.json","w").write(json.dumps(log_dict))
	xargs_cmd = "cat chromosomes.txt | xargs -i -P%s sh -c \"python %s/conformToRef.py {} %s\"" % (threads,script_dir,sample)
	subprocess.call(xargs_cmd,shell=True)


def cleanup(base_dir,sample):
	subprocess.call("rm %s.*.bed" % sample, shell=True)
	subprocess.call("rm %s.nosex" % sample, shell=True)
	if not os.path.isdir(base_dir+"/preimpute"):
		subprocess.call("mkdir "+base_dir+"/preimpute",shell=True)
	subprocess.call("mv %s* preimpute" % sample, shell=True)



def preprocess(args):
	flip_snp(args.ref,args.sample,args.threads)
	cleanup(args.base_dir,args.sample)

def impute(args):
	sample = args.sample
	chr = args.chr
	beagle_cmd = "java8 -Xmx50g -jar ~/software/beagle.03May16.862.jar gt=preimpute/%s.%s.preimpute.vcf ref=refVCF/chr%s.1kg.phase3.v5a.vcf.gz map=map/plink.chr%s.GRCh37.map out=%s.%s.imputed.vcf nthreads=%s" % (sample,chr,chr,chr,sample,chr,args.threads)
	subprocess.call(beagle_cmd,shell=True)

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_raw = subparsers.add_parser('preprocess', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('ref',help='RefFile')
parser_raw.add_argument('base_dir',help='RefFile')
parser_raw.add_argument('threads',help='RefFile')
parser_raw.set_defaults(func=preprocess)

parser_raw = subparsers.add_parser('impute', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('chr',help='RefFile')
parser_raw.add_argument('threads',help='RefFile')
parser_raw.set_defaults(func=impute)


args = parser.parse_args()
args.func(args)
