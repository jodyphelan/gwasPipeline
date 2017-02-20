import sys


script,famfile,metafile,metaname,outfile = sys.argv

pheno = {}
p_NA,p_1,p_0 = 0,0,0

for l in open(metafile):
	arr = l.rstrip().split()
	if arr[1]=="NA":
		p = "-9"
		p_NA+=1
	elif arr[1]==metaname:
		p = "1"
		p_1+=1
	else:
		p = "0"
		p_0+=1
	pheno[arr[0]] = p

with open(outfile,"w") as o:
	for l in open(famfile):
		arr = l.rstrip().split()
		if arr[1] not in pheno:
			print "ERROR: %s not found" % arr[1]
			quit()
		arr[5] = pheno[arr[1]]
		o.write("%s\n" % (" ".join(arr)))

if p_1==0:
	print "Warning: 0 subjects with positive phenotype"
print "Number of subjects:%s\nMissing:%s\nNegative:%s\nPositive:%s" % (sum([p_NA,p_1,p_0]),p_NA,p_0,p_1)
