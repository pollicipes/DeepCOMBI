#!/bin/bash

#set -Eeuxo pipefail;

# gwas catalog associations from https://www.ebi.ac.uk/gwas/docs/file-downloads
associ="gwas_catalog_v1.0.2-associations_e100_r2020-06-30.tsv"
# gwas ancestries
ances="gwas_catalog-ancestry_r2020-06-30.tsv"
# folder for the project
folder="/home/jrodriguez/Projects/combi_nn"
# results
dis_folder="snps_to_check_27_09_2019_vs_GWASCat_July2020"
# trait as reported in catalog associations EFO classifications
trait="bipolar disorder"
#trait=$1
# trait name as reported from Bettina
disease="BipolarDisorder";
#disease=$2

# remove those SNPs that could not be found (empty files, X chromosome or weird chromosome names):
cd ${folder}/${dis_folder}/"${disease}"_snps/
find . -size  0 -print0 | xargs -0 rm --
echo "Calculating LD..."
# make them a single vcf.
vcf-concat rs*.vcf > all"${disease}"SNPs.vcf
# sort
vcf-sort all"${disease}"SNPs.vcf> all"${disease}"SNPs_sort.vcf

# make some dirs and clean up single vcfs.
mkdir -p ${folder}/${dis_folder}/"${disease}"_snps/vcfs_snps;
mkdir -p ${folder}/${dis_folder}/"${disease}"_snps/ld;
mv rs* vcfs_snps;

# select the European individuals (n=85) genotypes to calculate LD and output it in bed format for plink
vcftools --vcf all"${disease}"SNPs_sort.vcf --keep  ${folder}/CEU_indv.txt --plink --out "${disease}"_CEU;

# keep only the snps that have an rs number (do not consider structural variants)
grep -v '^#' all"${disease}"SNPs_sort.vcf | awk -F'\t' -v OFS=$'\t' '{if($3 ~ /^rs/){print $3}}' > snps_to_keep.rs

# filter bed file to keep only the selected snps above
plink --file "${disease}"_CEU --extract snps_to_keep.rs --make-bed --out "${disease}"_CEU_rs

# for each snp calculate the LD R² that it has with other snps from the trait present at the Catalog 1Mb around @ R² > 0.2. 
# will report all the LD relationships over the threshold between each of our combi_nn snps and all the GWAS Catalog snps for that trait

while read snp; do
    echo ${snp}
    plink --bfile "${disease}"_CEU_rs --ld-snp ${snp} --ld-window-kb 500 --ld-window-r2 0.2 --r2 --ld-window 99999 --silent --out ld/${snp}
done < <(cat b)

# How many SNPs were actually tested. Out of the bettina snps find those that failed and subtract them from the initial number.
fail=$(grep -i 'Error' ld/*log | wc -l)
initial=$(wc -l b | cut -f1 -d' ')
total=$(( initial - fail ))
echo "$total SNPS WERE ACTUALLY TESTED FOR ${disease}";

# Once the above is completed, find all those files that have ld with SNPs, other than with self.
# If the snp is multiallelic or is fixed will not be analyzed.
# This returns a plink ld format output file with the snps and the LD value
# Adds also a header
cat <(echo -e "CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2") <(grep -v 'CHR' ld/*ld | awk -v OFS=$'\t' '{if($3 != $6){print $2,$3,$4,$5,$6,$7,$8}}') > ld/ld_candidates_${disease}

# we also need to evaluate those combi_nn selected snps that are directly present in the catalog; without any other LD proxy
# compare those that are already in both data sets
comm -12 a b > ld/direct_candidates_${disease}

# For each of the types of snps, find their respective entries in the GWAS Catalog associations
# note that for those ld proxies the snp that is featured and searched in the Catalog is the one in LD (column 6 of the {snp}.ld file), and not the SNP found by combi_nn itself.
cat ${folder}/header_GWASCat <(grep -wf ld/direct_candidates_${disease} ${folder}/${disease}.assoc) > ld/direct_candidates_${disease}.assoc
cat ${folder}/header_GWASCat <(grep -wf <(cut -f6 ld/ld_candidates_${disease} | sort | uniq) ${folder}/${disease}.assoc) > ld/ld_candidates_${disease}.assoc

echo "All done!"

# END
# # CLEAN...
# mkdir -p ${folder}"/mid_files";
# cd ${folder}
# mv /list_studies_* info_studies* studies_with_Europeans* *.assoc *_snps2019_hg38.bed mid_files/;
