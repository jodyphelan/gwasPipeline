# gwasPipeline

Scripts to process genotype data:
* Preprocessing
* Imputation


### Installation
Basic:
```
git clone --recursive https://github.com/jodyphelan/gwasPipeline.git
cd gwasPipeline
bash install_prerequisites.sh
```
If you want to use 1000 Genomes reference panels for imputation:
```
bash download_1kg.sh 4
```
Where 4 can be replaced with the number of threads available.

### Usage
If you have a plink formatted bim/bed/fam dataset e.g. `data.bim data.bed data.fam`:
```
/path/to/preprocessGWAS.py config /path/to/data
/path/to/preprocessGWAS.py init config.txt
```
This will use the 1KG data downloaded earlier and perform all the steps required before imputation. This can then be followed by imputation:
```
/path/to/preprocessGWAS.py impute final.preimpute <chromosome> <threads>
```

