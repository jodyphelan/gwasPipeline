#! /home/jody/software/anaconda2/bin/python
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
#	subprocess.call("rm %s.*.bed" % sample, shell=True)
	if not os.path.isdir(base_dir+"/preimpute"):
		subprocess.call("mkdir "+base_dir+"/preimpute",shell=True)
	subprocess.call("mv %s* preimpute" % sample, shell=True)



def preprocess(args):
	flip_snp("ref_fasta/ref_fasta.fa",args.sample,args.threads)
	cleanup(args.base_dir,args.sample)

def impute(args):
	if not os.path.isdir("imputed_vcf"):
		subprocess.call("mkdir imputed_vcf",shell=True)
	
	sample = args.sample
	chr = args.chr
	beagle_cmd = "java8 -Xmx50g -jar ~/software/beagle.03May16.862.jar gt=preimpute/%s.%s.preimpute.vcf ref=ref_vcf/ref_vcf_%s.vcf.gz map=ref_map/ref_map_%s.map out=imputed_vcf/%s.%s.imputed nthreads=%s" % (sample,chr,chr,chr,sample,chr,args.threads)
	subprocess.call(beagle_cmd,shell=True)

def init(args):
	config = [x.rstrip() for x in open(args.config).readlines()]
	data_dict = {}
	for line in config:
		arr = line.split()
		data_dict[arr[0]] = arr[1] 

	ref_check = "ref_fasta" in data_dict.keys()
	print "Checing for reference Fasta: %s" % ref_check
	vcf_test = sorted([int(x[8:]) for x in filter(lambda x: x[4:7]=="vcf", data_dict.keys())]) == range(1,23)
	print "Checing for reference VCF files: %s" % vcf_test
	map_test = sorted([int(x[8:]) for x in filter(lambda x: x[4:7]=="map", data_dict.keys())]) == range(1,23)
	print "Checing for reference map files: %s" % map_test
	
	print "Creating directory structure"
	if not os.path.isdir("ref_vcf"):
		subprocess.call("mkdir ref_vcf",shell=True)
	if not os.path.isdir("ref_map"):
		subprocess.call("mkdir ref_map",shell=True)
	if not os.path.isdir("ref_fasta"):
		subprocess.call("mkdir ref_fasta",shell=True)
	if not os.path.isdir("genotypes"):
		subprocess.call("mkdir genotypes",shell=True)

	fa_cmd = "ln -s %s ref_fasta/%s" % (data_dict["ref_fasta"],"ref_fasta.fa")
	subprocess.call(fa_cmd,shell=True)	
	for i in range(1,23):
		ref_vcf = "ref_vcf_%s" % i
		vcf_cmd = "ln -s %s ref_vcf/%s" % (data_dict[ref_vcf], ref_vcf+".vcf.gz")
		subprocess.call(vcf_cmd,shell=True)
		ref_map = "ref_map_%s" % i
		map_cmd = "ln -s %s ref_map/%s" % (data_dict[ref_map], ref_map+".map")
		subprocess.call(map_cmd,shell=True)	
	with open("chromosomes.txt","w") as o:
		if "chromosomes" in data_dict.keys():
			for i in data_dict["chromosomes"].split(","):
				o.write(i+"\n")
		else:
			for i in range(1,23):
				o.write(str(i)+"\n")
	print "Indexing VCFs"
	index_cmd = "cat chromosomes.txt | xargs -i -P20 sh -c \"~/software/bcftools-1.3.1/htslib-1.3.1/tabix ref_vcf/ref_vcf_{}.vcf.gz\""
	subprocess.call(index_cmd,shell=True)	


def fa2json(args):
	fa_dict = fa2dict(args.fasta)
	open(args.out,"w").write(json.dumps(fa_dict))




def mergeImputed(sample):
	chromosomes = [x.rstrip() for x in open("chromosomes.txt").readlines()]		
	with gzip.open(sample+".imputed.vcf.gz","wb") as o:
		test = True
		with gzip.open("imputed_vcf/"+sample+"."+chromosomes[0]+".imputed.vcf.gz","rb") as f:
			while test==True:
				line = f.readline()
				if line[0]=="#":
					o.write(line)
				else:
					test = False
		for i in tqdm(range(len(chromosomes))):
			chrom = chromosomes[i]
			for line in gzip.open("imputed_vcf/"+sample+"."+chrom+".imputed.vcf.gz","rb"):
				if line[0]!="#":
					o.write(line)
				
def imputationQC(args):
	vcf_file =  "imputed_vcf/"+args.sample+".all.vcf.gz" 
	preQC_prefix = "imputation_QC/"+args.sample+".preQC"
	postQC_prefix = "imputation_QC/"+args.sample+".postQC"
	#mergeImputed(args.sample)
	if not os.path.isdir("imputation_QC"):
		subprocess.call("mkdir imputation_QC",shell=True)
#	subprocess.call("plink --vcf "+vcf_file+" --make-bed --out "+preQC_prefix+" --const-fid",shell=True)
	subprocess.call("plink --bfile "+preQC_prefix+" --maf "+args.maf+" --geno "+args.geno+" --hwe "+args.hwe+" --make-bed --out "+postQC_prefix,shell=True)

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_raw = subparsers.add_parser('preprocess', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('base_dir',help='RefFile')
parser_raw.add_argument('threads',help='RefFile')
parser_raw.set_defaults(func=preprocess)

parser_raw = subparsers.add_parser('impute', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('chr',help='RefFile')
parser_raw.add_argument('threads',help='RefFile')
parser_raw.set_defaults(func=impute)

parser_raw = subparsers.add_parser('init', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('config',help='RefFile')
parser_raw.set_defaults(func=init)

parser_raw = subparsers.add_parser('fa2json', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('fasta',help='RefFile')
parser_raw.add_argument('out',help='RefFile')
parser_raw.set_defaults(func=fa2json)

parser_raw = subparsers.add_parser('imputeQC', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('maf',help='RefFile')
parser_raw.add_argument('geno',help='RefFile')
parser_raw.add_argument('hwe',help='RefFile')
parser_raw.set_defaults(func=imputationQC)

args = parser.parse_args()
args.func(args)
