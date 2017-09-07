mkdir 1kg
cd 1kg
seq 1 22 | xargs -i -P 20 sh -c "wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz.tbi"
seq 1 22 | xargs -i -P 20 sh -c "wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz"
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
unzip plink.GRCh37.map.zip
wget ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz
cd ../

echo "genotypes	$PWD/$1	
ref_fasta	$PWD/1kg/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz
ref_vcf_10	$PWD/1kg/chr10.1kg.phase3.v5a.vcf.gz
ref_vcf_11      $PWD/1kg/chr11.1kg.phase3.v5a.vcf.gz
ref_vcf_12      $PWD/1kg/chr12.1kg.phase3.v5a.vcf.gz
ref_vcf_13      $PWD/1kg/chr13.1kg.phase3.v5a.vcf.gz
ref_vcf_14      $PWD/1kg/chr14.1kg.phase3.v5a.vcf.gz
ref_vcf_15      $PWD/1kg/chr15.1kg.phase3.v5a.vcf.gz
ref_vcf_16      $PWD/1kg/chr16.1kg.phase3.v5a.vcf.gz
ref_vcf_17      $PWD/1kg/chr17.1kg.phase3.v5a.vcf.gz
ref_vcf_18      $PWD/1kg/chr18.1kg.phase3.v5a.vcf.gz
ref_vcf_19      $PWD/1kg/chr19.1kg.phase3.v5a.vcf.gz
ref_vcf_20      $PWD/1kg/chr20.1kg.phase3.v5a.vcf.gz
ref_vcf_21      $PWD/1kg/chr21.1kg.phase3.v5a.vcf.gz
ref_vcf_22      $PWD/1kg/chr22.1kg.phase3.v5a.vcf.gz
ref_vcf_1      $PWD/1kg/chr1.1kg.phase3.v5a.vcf.gz
ref_vcf_2      $PWD/1kg/chr2.1kg.phase3.v5a.vcf.gz
ref_vcf_3      $PWD/1kg/chr3.1kg.phase3.v5a.vcf.gz
ref_vcf_4      $PWD/1kg/chr4.1kg.phase3.v5a.vcf.gz
ref_vcf_5      $PWD/1kg/chr5.1kg.phase3.v5a.vcf.gz
ref_vcf_6      $PWD/1kg/chr6.1kg.phase3.v5a.vcf.gz
ref_vcf_7      $PWD/1kg/chr7.1kg.phase3.v5a.vcf.gz
ref_vcf_8      $PWD/1kg/chr8.1kg.phase3.v5a.vcf.gz
ref_vcf_9      $PWD/1kg/chr9.1kg.phase3.v5a.vcf.gz
ref_map_10	$PWD/1kg/plink.chr10.GRCh37.map
ref_map_11	$PWD/1kg/plink.chr11.GRCh37.map
ref_map_12	$PWD/1kg/plink.chr12.GRCh37.map
ref_map_13	$PWD/1kg/plink.chr13.GRCh37.map
ref_map_14	$PWD/1kg/plink.chr14.GRCh37.map
ref_map_15	$PWD/1kg/plink.chr15.GRCh37.map
ref_map_16	$PWD/1kg/plink.chr16.GRCh37.map
ref_map_17	$PWD/1kg/plink.chr17.GRCh37.map
ref_map_18	$PWD/1kg/plink.chr18.GRCh37.map
ref_map_19	$PWD/1kg/plink.chr19.GRCh37.map
ref_map_1	$PWD/1kg/plink.chr1.GRCh37.map
ref_map_20	$PWD/1kg/plink.chr20.GRCh37.map
ref_map_21	$PWD/1kg/plink.chr21.GRCh37.map
ref_map_22	$PWD/1kg/plink.chr22.GRCh37.map
ref_map_2	$PWD/1kg/plink.chr2.GRCh37.map
ref_map_3	$PWD/1kg/plink.chr3.GRCh37.map
ref_map_4	$PWD/1kg/plink.chr4.GRCh37.map
ref_map_5	$PWD/1kg/plink.chr5.GRCh37.map
ref_map_6	$PWD/1kg/plink.chr6.GRCh37.map
ref_map_7	$PWD/1kg/plink.chr7.GRCh37.map
ref_map_8	$PWD/1kg/plink.chr8.GRCh37.map
ref_map_9	$PWD/1kg/plink.chr9.GRCh37.map
" > config
