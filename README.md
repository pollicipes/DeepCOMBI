# DeepCOMBI
DeepCOMBI is a NN algorithm that exploits SNP relationships in GWAS datasets. This is a repository used to test the SNPs that DeepCOMBI produced. 

Paper under review in NAR Genomics & Bioinformatics.

DeepCOMBI: Explainable artificial intelligence for the analysis and discovery in genome-wide association studies https://www.biorxiv.org/content/10.1101/2020.11.06.371542v1 

The main script, `combi_check_SNPs_1.sh` takes a list of snps (output from DeepCOMBI) and checks whether they are on the GWAS Catalog, or if a LD proxy is found in there.

The `snps_to_check_27_09_2019` folder contains the output list for DeepCOMBI

Download the 1000 Genomes files, v5a.
