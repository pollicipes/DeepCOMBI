#!/bin/bash

# This controls the bash pipelines.
#set -Eeuxo pipefail;

### PART I ###

# gwas catalog associations from https://www.ebi.ac.uk/gwas/docs/file-downloads
associ="gwas_catalog_v1.0.2-associations_e100_r2020-06-30.tsv"
# gwas ancestries
ances="gwas_catalog-ancestry_r2020-06-30.tsv"
# folder for the project
folder="/home/jrodriguez/Projects/DeepCOMBI/"

# Data set to test:
group="lmm_lippert_020" #"lmm_lippert_STRONG_PRUNE" #"snps_to_check_27_09_2019"
r2thresh='0.2'
# results
#dis_folder="/snps_to_check_27_09_2019_vs_GWASCat_July2020/"
dis_folder="/"${group}"_vs_GWASCat_July2020/"

# trait as reported in catalog associations EFO classifications
#trait="bipolar disorder"
trait=$1
# trait name as reported from Bettina
#disease="BipolarDisorder";
disease=$2

echo "Starting..."

# Make the folder for benchmark
wd=${folder}"/"${dis_folder}"/"${disease}"/";
mkdir -p ${wd}

# Format headers and set up
# Save GWAS catalog header for later show with the results

head -n1 ${folder}/${associ} > ${folder}/header_GWASCat
echo "Identifying the list of studies for..." ${disease}

# Identify the list of studies for this trait from the GWAS Catalog association file
awk -F'\t' -v OFS='\t' -v t="${trait}" '$35 == t'  ${folder}/${associ} | cut -f37 | sort | uniq > ${wd}/list_studies_${disease}

echo "Searching population/ancestry info..."    
# Search the list of studies exclusively for the trait and their pop info
LC_ALL=C grep -wf ${wd}/list_studies_${disease} ${ances} > ${wd}/info_studies_${disease}

# Discard those studies performed exclusively in East Asians
# Remove also the WTCCC study itself from the list of studies returning SNPs for the trait in europeans
awk  -F'\t' -v OFS='\t' '$9 != "East Asian"' ${wd}/info_studies_${disease} |  awk -F'\t' -v OFS='\t' '$2 != "17554300"' | cut -f1 | sort | uniq > ${wd}/studies_with_Europeans_${disease}

# Full info for studies selected
LC_ALL=C grep -wf ${wd}/studies_with_Europeans_${disease} ${folder}"/"${associ} > ${wd}/${disease}.assoc

# SNPs coming from studies in Europeans
LC_ALL=C grep -wf ${wd}/studies_with_Europeans_${disease} ${folder}/${associ} | cut -f2,12,13,22 | sort -nk1,1 | awk '{print "chr"$2"\t"$3"\t"$3+1"\t"$4"\t"$1}' | sort | uniq > ${wd}/${disease}_snps2019_hg38.bed

# Format SNPs for liftOver. !!!REQUIRES THE chr* CHROMOSOME FORMAT!!!
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4}' ${wd}/${disease}_snps2019_hg38.bed | sort | uniq | grep -v '_\|;\|X\|:\|x' | sortBed > ${wd}/${disease}_snps2019_hg38_format.bed
 
# Get the positions for the input SNPs from DeepCOMBI
while read R; do 
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -Dhg19 -N -e "select concat(chrom,' ',chromEnd,' ',chromEnd+1,' ',name ) from snp142 where name='${R}'"; 
    # use this format to work with the Lippert test.
done < <(awk '{print $2}' ${folder}/${group}/all.${disease}.txt) > ${wd}/snps_combi_${disease}.l.bed

#done < <(awk '{print $2}' ${folder}/${group}/${disease}"_upto_10-4.txt") > ${wd}/snps_combi_${disease}.l.bed


# Remove weird HLA chromosome mappings
grep -v '_' ${wd}/snps_combi_${disease}.l.bed > ${wd}/snps_combi_${disease}.bed

# Format DeepCOMBI SNPs into a bed format
sed -i -e 's/chr//g' -e 's/ /\t/g' ${wd}/snps_combi_${disease}.bed;

# LIFTOVER CATALOG SNPS (requires the binary for liftOver, dowloadable from UCSC)
/home/jrodriguez/cluster_projects/metaloci/liftOver ${wd}/${disease}_snps2019_hg38_format.bed ${folder}/hg38ToHg19.over.chain.gz ${wd}/${disease}_hg38_to_hg19.bed ${wd}/unmapped_${disease}_hg38_to_hg19

# Now remove the chr* format.
sed -i 's/chr//g' ${wd}/${disease}_hg38_to_hg19.bed

# Concatenate them altogether
cat ${wd}/${disease}_hg38_to_hg19.bed  ${wd}/snps_combi_${disease}.bed | sortBed > ${wd}/full_${disease}.bed

# Make a dir to store VCFs...
mkdir -p ${wd}"/snps";
cd ${wd}"/snps";

# Get the corresponding VCFs for the SNPs in the GWAS Catalog. Needed for later calculating LD with our DeepCOMBI SNPs.
while read -r ch st end snp; do 

    ### If you want to use the tabix remote version, use this line below. But this is slower and genrally gives connectivity problems... 
    #tabix -f http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${ch}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz ${ch}:${st}-${end} > ${wd}"/snps/"${snp}.vcf; 

    ### Use better instead the local files, which you should download beforehand. Note that the current 1000G version is v5b (updated with some minor changes, irrelevant for this), 
    ### ...but we have in local database the v5a version which is fine.
    tabix -f -h ~/scratch/chr_files_1KG_CEU_MAF001/ALL.chr${ch}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${ch}:${st}-${end} > ${wd}"/snps/"${snp}.vcf; 
done < ${wd}/full_${disease}.bed

echo "REACHES; DONE!"

cd ${folder};

# A dir to store LD data.
mkdir -p ${wd}"/ld";

### PART II ###
## Its adviceable to run first part I, and then continue with this part II
## Part I is basically for downloading and getting base info for the SNPs
## Part II calculates LD between our DeepCOMBI and the GWAS SNPs.

# Remove those SNPs that could not be found (empty files, X chromosome or weird chromosome names):
find ${wd}"/snps/" -size 0 -print0 | xargs -0 rm -- 

echo "Calculating LD..."
# Put all the SNPs into a single vcf.
vcf-concat ${wd}"/snps/"rs*.vcf > ${wd}"/snps/"all${disease}SNPs.vcf
# sort
vcf-sort ${wd}"/snps/"all${disease}SNPs.vcf > ${wd}"/snps/"all${disease}SNPs_sort.vcf
# Select the European individuals (n=85) genotypes to calculate LD and output it in bed format for plink
vcftools --vcf ${wd}"/snps/"all${disease}SNPs_sort.vcf --keep  ${folder}/CEU_indv.txt --plink --out ${wd}"/ld/"${disease}_CEU;

# Keep only the snps that have an rs number (do not consider structural variants)
grep -v '^#' ${wd}"/snps/"all${disease}SNPs_sort.vcf | awk -F'\t' -v OFS=$'\t' '{if($3 ~ /^rs/){print $3}}' > ${wd}"/snps/snps_to_keep.rs"

# Filter bed file to keep only the selected snps above
plink --file ${wd}"/ld/"${disease}_CEU --extract ${wd}"/snps/snps_to_keep.rs" --make-bed --out ${wd}"/ld/"${disease}_CEU_rs

# For each snp calculate the LD R² that it has with other snps from the trait present at the Catalog 1Mb around @ R² > 0.2. 
# will report all the LD relationships over the threshold between each of our combi_nn snps and all the GWAS Catalog snps for that trait
while read snp; do
    echo ${snp}
    plink --bfile ${wd}"/ld/"${disease}_CEU_rs --ld-snp ${snp} --ld-window-kb 500 --ld-window-r2 ${r2thresh} --r2 --ld-window 99999 --silent --out ${wd}"/ld/"${snp}
done < <(awk -F'\t' '{print $4}' ${wd}/snps_combi_${disease}.bed)


# How many SNPs were actually tested. 
# Out of the DeepCOMBI SNPs, get those that failed and subtract them from the initial number.
fail=$(grep -i 'Error' ${wd}"/ld/"*.log | wc -l)
initial=$(wc -l ${wd}/snps_combi_${disease}.bed | cut -f1 -d' ')
total=$(( initial - fail ))
echo "$total SNPS WERE ACTUALLY TESTED FOR ${disease}";

# Once the above is completed, find all those files that have ld with SNPs, other than with self.
# If the SNP is multiallelic or is fixed, it will not be analyzed.
# This returns a plink ld format output file with the snps and the LD value
# Adds also a header
cat <(echo -e "CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2") <(grep -v 'CHR' ${wd}"/ld/"*.ld | awk -v OFS=$'\t' '{if($3 != $6){print $2,$3,$4,$5,$6,$7,$8}}') > ${wd}"/ld/ld_candidates_"${disease}

# We also need to evaluate those DeepCOMBI selected SNPs that 
# are directly present in the catalog; without any other LD proxy
# Compare those that are already in both data sets
comm -12 <(cut -f4 ${wd}/${disease}_hg38_to_hg19.bed | sort) <(cut -f4 ${wd}/snps_combi_${disease}.bed | sort) > ${wd}"/ld/direct_candidates_"${disease}

# For each of the types of SNPs, find their respective entries in the GWAS Catalog associations...
# ...only to complete more info.
# Note that for those ld proxies the SNP that is featured and searched in the Catalog is the one in LD (column 6 of the {snp}.ld file), and not the SNP found by DeepCOMBI itself.
## Basically, we are trying to identify which are the SNPs that appear in the Catalog which are in LD with our candidate.
cat ${folder}/header_GWASCat <(grep -wf ${wd}"/ld/direct_candidates_"${disease} ${wd}/${disease}.assoc) > ${wd}"/ld/direct_candidates_"${disease}.assoc
cat ${folder}/header_GWASCat <(grep -wf <(cut -f6 ${wd}"ld/ld_candidates_"${disease} | sort | uniq) ${wd}/${disease}.assoc) > ${wd}"/ld/ld_candidates_"${disease}.assoc

# We need to add these two lines in order to get the really validated ones.
## First, we get those candidates that were truly found in the catalog, which are in LD with some of our candidates.
cut -f22 ${wd}"/ld/ld_candidates_"${disease}.assoc | sed 1d | uniq > ${wd}"/ld/"ld.found.assoc
### Then, we search for them in the LD full list to get the original candidates, which is then merged with the direct candidates. 
cat <(grep -wf ${wd}"/ld/"ld.found.assoc ${wd}"ld/ld_candidates_"${disease} | cut -f3 | sort | uniq) ${wd}"/ld/direct_candidates_"${disease} | sort | uniq > ${wd}"/ld/final_all_"${disease}


echo "All done!"

# END
