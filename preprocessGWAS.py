#! /usr/bin/python
import sys
import subprocess
import gzip
import os
import argparse
import json

script_dir = os.path.dirname(os.path.realpath(__file__))

beagle = "%s/bin/beagle.08Jun17.d8b.jar" % script_dir
tabix = "%s/bin/tabix" % script_dir
plink = "%s/bin/plink" % script_dir

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
	with gzip.open(filename) if filename[-3:]==".gz" else open(filename) as f:
		for i in (range(file_len(filename))):
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


def flip_snp(ref_file,sample):
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
		for i in (range(file_len(bim_file))):
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
			if chr not in fa_dict or int(pos)-1>len(fa_dict[chr]):
				print len(fa_dict[chr])
				print "%s\t%s" % (chr,pos)
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
	plink_cmd = "%s --allow-no-sex --bfile %s --exclude %s.exclude.txt --flip %s.reverse.txt --make-bed --out %s" % (plink,"genotypes/"+sample,sample,sample,sample+".flipped")
	subprocess.call(plink_cmd,shell=True)
	open(sample+".flipFilter.json","w").write(json.dumps(log_dict))

def split_by_chr(base_dir,sample,threads):
	xargs_cmd = "cat chromosomes.txt | xargs -i -P%s sh -c \"python %s/conformToRef.py {} %s\"" % (threads,script_dir,sample)
	subprocess.call(xargs_cmd,shell=True)

def cleanup(base_dir,sample):
#	subprocess.call("rm %s.*.bed" % sample, shell=True)
	if not os.path.isdir(base_dir+"/preimpute"):
		subprocess.call("mkdir "+base_dir+"/preimpute",shell=True)
	subprocess.call("mv %s* preimpute" % sample, shell=True)



def preprocess(args):
	flip_snp("ref_fasta/ref_fasta.fa",args.sample)
	cleanup(args.base_dir,args.sample)

def impute(args):
	if not os.path.isdir("imputed_vcf"):
		subprocess.call("mkdir imputed_vcf",shell=True)
	
	sample = args.sample
	chrom = args.chr
	print "python %s/conformToRef.py %s preimpute/%s" % (script_dir,chrom,sample)
	subprocess.call("python %s/conformToRef.py %s preimpute/%s" % (script_dir,chrom,sample),shell=True)
	beagle_cmd = "java -Xmx50g -jar %s gt=preimpute/%s.%s.preimpute.vcf ref=ref_vcf/ref_vcf_%s.vcf.gz map=ref_map/ref_map_%s.map out=imputed_vcf/%s.%s.imputed nthreads=%s" % (beagle,sample,chrom,chrom,chrom,sample,chrom,args.threads)
	print beagle_cmd
	subprocess.call(beagle_cmd,shell=True)

def init(args):
	data_dict = {}
	for line in [x.rstrip() for x in open(args.config).readlines()]:
		arr = line.split()
		data_dict[arr[0]] = arr[1] 

	
	num_chromosomes = len([d for d in data_dict.keys() if "ref_vcf" in d])+1


	ref_check = "ref_fasta" in data_dict.keys()
	print "Checking for reference Fasta: %s" % ref_check
	vcf_test = sorted([int(x[8:]) for x in filter(lambda x: x[4:7]=="vcf", data_dict.keys())]) == range(1,num_chromosomes)
	print "Checking for reference VCF files: %s" % vcf_test
	map_test = sorted([int(x[8:]) for x in filter(lambda x: x[4:7]=="map", data_dict.keys())]) == range(1,num_chromosomes)
	print "Checking for reference map files: %s" % map_test
	
	print "Creating directory structure"
	for x in ["ref_vcf","ref_map","ref_fasta","genotypes","plots","logs"]:
		if not os.path.isdir(x):
			subprocess.call("mkdir %s"%x,shell=True)
	fa_cmd = "ln -s %s ref_fasta/%s" % (data_dict["ref_fasta"],"ref_fasta.fa")
	subprocess.call(fa_cmd,shell=True)	
	for i in range(1,num_chromosomes):
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
			for i in range(1,num_chromosomes):
				o.write(str(i)+"\n")
	print "Indexing VCFs"
	index_cmd = "cat chromosomes.txt | xargs -i -P20 sh -c \"%s ref_vcf/ref_vcf_{}.vcf.gz\"" % tabix
	subprocess.call(index_cmd,shell=True)	

	file_paths = data_dict["genotypes"].split(",")
	file_prefix = [x.split("/")[-1] for x in file_paths]
	print "Linking Genotypes"
	for i,x in enumerate(file_paths):
		y = file_prefix[i]
		subprocess.call("ln -s %s.bed genotypes/%s.bed" % (x,y),shell=True)
		subprocess.call("ln -s %s.bim genotypes/%s.bim" % (x,y),shell=True)
		subprocess.call("ln -s %s.fam genotypes/%s.fam" % (x,y),shell=True)

	with open("runAnalysis.sh","w") as o:
		for i,x in enumerate(file_paths):
			y = file_prefix[i]	
			o.write("%s preprocess %s .\n" % (sys.argv[0],y))
		temp = ",".join(["preimpute/"+x+".flipped" for x in file_prefix])
		if len(file_paths)>1:
			o.write("%s/relaxed_merge.py %s final.preimpute --pca\n" % (script_dir,temp))
		else:	
			for x in ["fam","bim","bed"]:
				o.write("ln -s preimpute/%s.flipped.%s preimpute/final.preimpute.%s\n" % (file_prefix[0],x,x))
		subprocess.call("bash runAnalysis.sh",shell=True)
						

	
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
		for i in (range(len(chromosomes))):
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
	subprocess.call("%s --bfile "+preQC_prefix+" --maf "+args.maf+" --geno "+args.geno+" --hwe "+args.hwe+" --make-bed --out %s" % (plink,postQC_prefix),shell=True)

def download_1kg():
	os.chdir(script_dir)
	os.mkdir("1kg")
	os.chdir("1kg")
	subprocess.call("seq 1 22 | xargs -i -P 20 sh -c \"wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz.tbi\" && seq 1 22 | xargs -i -P 20 sh -c \"wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz\" && wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip && unzip plink.GRCh37.map.zip && wget ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz",shell=True)

def generate_1kg_config(args):
	paths = {"d":script_dir,"f":os.path.realpath(args.prefix)}
	config = """genotypes	%(f)s	
ref_fasta	%(d)s/1kg/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz
ref_vcf_10	%(d)s/1kg/chr10.1kg.phase3.v5a.vcf.gz
ref_vcf_11      %(d)s/1kg/chr11.1kg.phase3.v5a.vcf.gz
ref_vcf_12      %(d)s/1kg/chr12.1kg.phase3.v5a.vcf.gz
ref_vcf_13      %(d)s/1kg/chr13.1kg.phase3.v5a.vcf.gz
ref_vcf_14      %(d)s/1kg/chr14.1kg.phase3.v5a.vcf.gz
ref_vcf_15      %(d)s/1kg/chr15.1kg.phase3.v5a.vcf.gz
ref_vcf_16      %(d)s/1kg/chr16.1kg.phase3.v5a.vcf.gz
ref_vcf_17      %(d)s/1kg/chr17.1kg.phase3.v5a.vcf.gz
ref_vcf_18      %(d)s/1kg/chr18.1kg.phase3.v5a.vcf.gz
ref_vcf_19      %(d)s/1kg/chr19.1kg.phase3.v5a.vcf.gz
ref_vcf_20      %(d)s/1kg/chr20.1kg.phase3.v5a.vcf.gz
ref_vcf_21      %(d)s/1kg/chr21.1kg.phase3.v5a.vcf.gz
ref_vcf_22      %(d)s/1kg/chr22.1kg.phase3.v5a.vcf.gz
ref_vcf_1      %(d)s/1kg/chr1.1kg.phase3.v5a.vcf.gz
ref_vcf_2      %(d)s/1kg/chr2.1kg.phase3.v5a.vcf.gz
ref_vcf_3      %(d)s/1kg/chr3.1kg.phase3.v5a.vcf.gz
ref_vcf_4      %(d)s/1kg/chr4.1kg.phase3.v5a.vcf.gz
ref_vcf_5      %(d)s/1kg/chr5.1kg.phase3.v5a.vcf.gz
ref_vcf_6      %(d)s/1kg/chr6.1kg.phase3.v5a.vcf.gz
ref_vcf_7      %(d)s/1kg/chr7.1kg.phase3.v5a.vcf.gz
ref_vcf_8      %(d)s/1kg/chr8.1kg.phase3.v5a.vcf.gz
ref_vcf_9      %(d)s/1kg/chr9.1kg.phase3.v5a.vcf.gz
ref_map_10	%(d)s/1kg/plink.chr10.GRCh37.map
ref_map_11	%(d)s/1kg/plink.chr11.GRCh37.map
ref_map_12	%(d)s/1kg/plink.chr12.GRCh37.map
ref_map_13	%(d)s/1kg/plink.chr13.GRCh37.map
ref_map_14	%(d)s/1kg/plink.chr14.GRCh37.map
ref_map_15	%(d)s/1kg/plink.chr15.GRCh37.map
ref_map_16	%(d)s/1kg/plink.chr16.GRCh37.map
ref_map_17	%(d)s/1kg/plink.chr17.GRCh37.map
ref_map_18	%(d)s/1kg/plink.chr18.GRCh37.map
ref_map_19	%(d)s/1kg/plink.chr19.GRCh37.map
ref_map_1	%(d)s/1kg/plink.chr1.GRCh37.map
ref_map_20	%(d)s/1kg/plink.chr20.GRCh37.map
ref_map_21	%(d)s/1kg/plink.chr21.GRCh37.map
ref_map_22	%(d)s/1kg/plink.chr22.GRCh37.map
ref_map_2	%(d)s/1kg/plink.chr2.GRCh37.map
ref_map_3	%(d)s/1kg/plink.chr3.GRCh37.map
ref_map_4	%(d)s/1kg/plink.chr4.GRCh37.map
ref_map_5	%(d)s/1kg/plink.chr5.GRCh37.map
ref_map_6	%(d)s/1kg/plink.chr6.GRCh37.map
ref_map_7	%(d)s/1kg/plink.chr7.GRCh37.map
ref_map_8	%(d)s/1kg/plink.chr8.GRCh37.map
ref_map_9	%(d)s/1kg/plink.chr9.GRCh37.map
""" % paths	
	open("config.txt","w").write(config)
	
parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_raw = subparsers.add_parser('preprocess', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('sample',help='RefFile')
parser_raw.add_argument('base_dir',help='RefFile')
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

parser_raw = subparsers.add_parser('config', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_raw.add_argument('prefix',help='RefFile')
parser_raw.set_defaults(func=generate_1kg_config)


args = parser.parse_args()
args.func(args)
