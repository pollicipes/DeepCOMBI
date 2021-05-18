# DeepCOMBI
DeepCOMBI is a NN algorithm that exploits SNP relationships in GWAS datasets. This is a repository used to test the SNPs that DeepCOMBI produced. 

Paper under review in NAR Genomics & Bioinformatics.

DeepCOMBI: Explainable artificial intelligence for the analysis and discovery in genome-wide association studies https://www.biorxiv.org/content/10.1101/2020.11.06.371542v1 

The main script, `combi_check_SNPs_1.sh` takes a list of snps (output from DeepCOMBI) and checks whether they are on the GWAS Catalog, or if a LD proxy is found in there.

The `snps_to_check_27_09_2019` folder contains the output list for DeepCOMBI

Download the 1000 Genomes files, v5a.

---

This scripts validate snps found by combi https://www.nature.com/articles/srep36671 new neural network approach.

Download last release of the GWAS Catalog, both the associations file gwas_catalog_v1.0.2-associations_e100_r2020-06-30.tsv and the ancestry file gwas_catalog-ancestry_r2020-06-30.tsv. By July 2020 these two were the last ones.

Make a working folder and put there: 1- The script ./combi_check_SNPs.sh 2- Both GWAS Catalog files. Check and open the script to set the name of the correct Catalog files you are using. 3- Another folder containing the COMBI output SNPs as a list.

The script takes as a first argument the trait as reported in the EFO terms field of the catalog. For type 1 diabetes, you would write: "type I diabetes mellitus"

As a second argument add the name as you have it written in the outputs of the COMBI. In this case it would be "Type1Diabetes". Run it like this:

./combi_check_SNPs.sh "type I diabetes mellitus" "Type1Diabetes"

And the script runs and performs the search.
