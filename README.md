# DeepCOMBI
DeepCOMBI is a NN algorithm that exploits SNP relationships in GWAS datasets. This is a repository with the script used to test the SNPs that DeepCOMBI produced. 

Paper under review in NAR Genomics & Bioinformatics.

DeepCOMBI: Explainable artificial intelligence for the analysis and discovery in genome-wide association studies https://www.biorxiv.org/content/10.1101/2020.11.06.371542v1 

The main script, `combi_check_SNPs_1.sh` takes a list of snps (output from DeepCOMBI) and checks whether they are on the GWAS Catalog, or if a LD proxy is found in there.

### Some prerequisites

You'll also need the binaries for `tabix`, `liftOver`, `vcftools` and `PLINK`. Bioinformatics' classic ones. ;) 

The `snps_to_check_27_09_2019` folder contains the list of candidates found by DeepCOMBI and which were run through the pipeline.

Use a reference panel to calculate LD and validate discoveries through LD, for example the 1000 Genomes files, v5a. In this repo it is provided the list of CEU individuals `CEU_indv.txt` that PLINK will need to filter only the European individuals.

Download last release of the GWAS Catalog, both the associations file `gwas_catalog_v1.0.2-associations_e100_r2020-06-30.tsv` and the ancestry file `gwas_catalog-ancestry_r2020-06-30.tsv`. Ancestry file is needed to get an estimate of the majoritarian population in a GWAS study. Given that we are validating European data, we'll only keep European populations' GWAS. By July 2020 the files here were the last two ones available. Note that last versions of the Catalog are in hg38. Our WTCCC dataset is in hg37, hence why we need `liftOver`. Chain file `hg38ToHg19.over.chain.gz` is also provided here.

### Testing the discoveries of DeepCOMBI against GWAS Catalog. 
The basic premise behind this test is: If by 2007, when the WTCCC was published, we had applied the DeepCOMBI to the dataset, how many hit would have been truly discovered (and validated) in the future? So, we aim to do this. Hits found by DeepCOMBI to be tested can be found in `snps_to_check_27_09_2019`.

### DeepCOMBI results comparison against Lippert et al., (2013)'s LMM method
In order to compare the performance against some other algorithms, we evaluated the performance against the LMM method by [Lippert et al., 2013](https://www.nature.com/articles/srep01099). In supplementary table 2, where they report a list of 573 SNPs that they found with their LMM method. This file is provided here for convenience: `lmm_univar.txt`. Some of the SNPs in there were not found by LMM strictly, so we discarded those SNPs. The resulting file is `lmm_univar.filter.txt`, also provided here. Note that ~85% of the SNPs they found map in the HLA regions. This file was further crossed with the list of the SNPs that went through COMBI, as we further filtered the `lmm_univar.filter.txt` file to contain only for each disease those SNPs that run through our pipeline, to make a fair comparison and avoid ascertainment bias. See the [original COMBI article](https://www.nature.com/articles/srep36671) for more info on this. 
Then, we LD-pruned those SNPs from Lippert using `PLINK`'s command `--indep-pairwise 2 1 0.8` to get the independent signals that [Lippert et al., 2013](https://www.nature.com/articles/srep01099) obtained, because they reported hundreds of SNPs in LD in those supplementary tables, which wouldn't make a fair comparison. The final lists per diseases of SNPs used to compare are in the `lippert_020` folder.


### Validation step 

Download the script you need `combi_check_SNPs_1.sh` and:

Make a working folder and put there: 
1- The script ./combi_check_SNPs.sh 
2- Both GWAS Catalog files. Check and open the script to set the name of the correct Catalog files you are using. 
3- Another folder containing the COMBI output SNPs as a list (See `snps_to_check_27_09_2019`) for an example.

The script takes as a first argument the trait as reported in the EFO terms field of the catalog. For type 1 diabetes, you would write: "type I diabetes mellitus"

As a second argument add the name as you have it written in the outputs of the COMBI. In this case it would be "t1d". **Check the correct paths for the files, binaries, etc...**

Then, run it like this:

`./combi_check_SNPs.sh "type I diabetes mellitus" "t1d"`

And the script runs and performs the search... You can check the results in the folder created wit the suffix: `*_vs_GWASCat_July2020`. In there, there will be one folder per disease, and once there, there's an `./ld/` folder, where you can see the snps found in the files with the \*candidates\* prefix, either direct candidates (the proposed SNPs are in the catalog like that) of through ld, where a SNP in LD >0.2 is found. Check the found ones in the `final_all_${theDiseaseTested}` files, which congregates all of the found original SNPs.
