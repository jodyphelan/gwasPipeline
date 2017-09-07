mkdir 1kg
cd 1kg

seq 1 22 | xargs -i -P 20 sh -c "wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz.tbi" 

seq 1 22 | xargs -i -P 20 sh -c "wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr{}.1kg.phase3.v5a.vcf.gz"

wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip 
unzip plink.GRCh37.map.zip 

wget ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz

cd ../
