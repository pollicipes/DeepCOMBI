#!/bin/bash

# This controls the bash pipelines
set -Eeuxo pipefail;

# gwas catalog associations from https://www.ebi.ac.uk/gwas/docs/file-downloads
associ="gwas_catalog_v1.0.2-associations_e100_r2020-06-30.tsv"
# gwas ancestries
ances="gwas_catalog-ancestry_r2020-06-30.tsv"
# folder for the project
folder="/home/jrodriguez/Projects/combi_nn/"
# results
dis_folder="/snps_to_check_27_09_2019_vs_GWASCat_July2020/"
# trait as reported in catalog associations EFO classifications
#trait="bipolar disorder"
trait=$1
# trait name as reported from Bettina
#disease="BipolarDisorder";
disease=$2

echo "Starting..."

wd=${folder}"/"${dis_folder}"/"${disease}"/";
# Move to the folder of benchmark
mkdir -p ${wd}

# headers and set up
# save gwas catalog header for later show on the results

head -n1 ${folder}/${associ} > ${folder}/header_GWASCat
echo "Identifying the list of studies for..." ${disease}

# identify the list of studies for this trait from the GWAS Catalog association file
awk -F'\t' -v OFS='\t' -v t="${trait}" '$35 == t'  ${folder}/${associ} | cut -f37 | sort | uniq > ${wd}/list_studies_${disease}

echo "Searching population/ancestry info..."    
# Search the list of studies exclusively for the trait and their pop info
grep -wf ${wd}/list_studies_${disease} ${ances} > ${wd}/info_studies_${disease}

# Discard those studies performed exclusively in East Asians
# Remove also the WTCCC study itself from the list of studies returning SNPs for the trait in europeans
#### CONSIDER REMOVING AFRICANS???
awk  -F'\t' -v OFS='\t' '$9 != "East Asian"' ${wd}/info_studies_${disease} |  awk -F'\t' -v OFS='\t' '$2 != "17554300"' | cut -f1 | sort | uniq > ${wd}/studies_with_Europeans_${disease}

# full info for studies selected
grep -wf ${wd}/studies_with_Europeans_${disease} ${folder}"/"${associ} > ${wd}/${disease}.assoc

# SNPs coming from studies in Europeans
grep -wf ${wd}/studies_with_Europeans_${disease} ${folder}/${associ} | cut -f2,12,13,22 | sort -nk1,1 | awk '{print "chr"$2"\t"$3"\t"$3+1"\t"$4"\t"$1}' | sort | uniq > ${wd}/${disease}_snps2019_hg38.bed

# FORMAT THESE GWAS CAT SNPS TO DO LIFTOVER. REQUIRES THE chr FORMAT
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4}' ${wd}/${disease}_snps2019_hg38.bed | sort | uniq | grep -v '_\|;\|X\|:\|x' | sortBed > ${wd}/${disease}_snps2019_hg38_format.bed
 
#### GET THE POSITION FOR THE INPUT SNPs FROM BETTINA. 
while read R; do 
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -Dhg19 -N -e "select concat(chrom,' ',chromEnd,' ',chromEnd+1,' ',name ) from snp142 where name='${R}'"; 
done < <(awk '{print $2}' ${folder}/snps_to_check_27_09_2019/${disease}_upto_10-4.txt | sed 1d) > ${wd}/snps_combi_${disease}.bed

# FORMAT BETTINA SNPS AND PUT THEM IN A BED FORMAT.
sed -i -e 's/chr//g' -e 's/ /\t/g' ${wd}/snps_combi_${disease}.bed;

## LIFTOVER CATALOG SNPS 
/home/jrodriguez/cluster_projects/metaloci/liftOver ${wd}/${disease}_snps2019_hg38_format.bed ${folder}/hg38ToHg19.over.chain.gz ${wd}/${disease}_hg38_to_hg19.bed ${wd}/unmapped_${disease}_hg38_to_hg19

# WE SHOULD NOW REMOVE THE chr FORMAT FROM THE CATALOG SNPs.
sed -i 's/chr//g' ${wd}/${disease}_hg38_to_hg19.bed

cat ${wd}/${disease}_hg38_to_hg19.bed  ${wd}/snps_combi_${disease}.bed | sortBed > ${wd}/full_${disease}.bed

mkdir -p ${wd}"/snps";

###GET THE VCFS FOR THE SNPS IN THE CATALOG
cd ${wd}"/snps";
while read -r ch st end snp; do 
    tabix -f -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${ch}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${ch}:${st}-${end} > ${wd}"/snps/"${snp}.vcf; 
done < ${wd}/full_${disease}.bed
cd ${folder};

# CALCULATING HOW MANY WERE FETCHED IN VCFs FROM THE CATALOG:
#grep -v '^#' rs*.vcf | cut -f3 | grep -v 'esv' | sort | uniq  > invcfs.snps
#cut -f4 ../BipolarDisorder_hg38_to_hg19.bed |  sort | uniq > orig.snps
#comm -12 orig.snps invcfs.snps | wc -l

mkdir -p ${wd}"/ld";

### PART II

# remove those SNPs that could not be found (empty files, X chromosome or weird chromosome names):
find ${wd}"/snps/" -empty -type f -delete

echo "Calculating LD..."
# make them a single vcf.
vcf-concat ${wd}"/snps/"rs*.vcf > ${wd}"/snps/"all${disease}SNPs.vcf
# sort
vcf-sort ${wd}"/snps/"all${disease}SNPs.vcf > ${wd}"/snps/"all${disease}SNPs_sort.vcf

# select the European individuals (n=85) genotypes to calculate LD and output it in bed format for plink
vcftools --vcf ${wd}"/snps/"all${disease}SNPs_sort.vcf --keep  ${folder}/CEU_indv.txt --plink --out ${wd}"/ld/"${disease}_CEU;

# keep only the snps that have an rs number (do not consider structural variants)
grep -v '^#' ${wd}"/snps/"all${disease}SNPs_sort.vcf | awk -F'\t' -v OFS=$'\t' '{if($3 ~ /^rs/){print $3}}' > ${wd}"/snps/snps_to_keep.rs"

# filter bed file to keep only the selected snps above
plink --file ${wd}"/ld/"${disease}_CEU --extract ${wd}"/snps/snps_to_keep.rs" --make-bed --out ${wd}"/ld/"${disease}_CEU_rs

# for each snp calculate the LD R² that it has with other snps from the trait present at the Catalog 1Mb around @ R² > 0.2. 
# will report all the LD relationships over the threshold between each of our combi_nn snps and all the GWAS Catalog snps for that trait
# exit;

while read snp; do
    echo ${snp}
    plink --bfile ${wd}"/ld/"${disease}_CEU_rs --ld-snp ${snp} --ld-window-kb 500 --ld-window-r2 0.2 --r2 --ld-window 99999 --silent --out ${wd}"/ld/"${snp}
done < <(awk -F'\t' '{print $4}' ${wd}/snps_combi_${disease}.bed)


# How many SNPs were actually tested. Out of the bettina snps find those that failed and subtract them from the initial number.
fail=$(grep -i 'Error' ${wd}"/ld/"*.log | wc -l)
initial=$(wc -l ${wd}/snps_combi_${disease}.bed | cut -f1 -d' ')
total=$(( initial - fail ))
echo "$total SNPS WERE ACTUALLY TESTED FOR ${disease}";

# Once the above is completed, find all those files that have ld with SNPs, other than with self.
# If the snp is multiallelic or is fixed will not be analyzed.
# This returns a plink ld format output file with the snps and the LD value
# Adds also a header
cat <(echo -e "CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2") <(grep -v 'CHR' ${wd}"/ld/"*.ld | awk -v OFS=$'\t' '{if($3 != $6){print $2,$3,$4,$5,$6,$7,$8}}') > ${wd}"/ld/ld_candidates_"${disease}

# we also need to evaluate those combi_nn selected snps that are directly present in the catalog; without any other LD proxy
# compare those that are already in both data sets

comm -12 <(cut -f4 ${wd}/${disease}_hg38_to_hg19.bed | sort) <(cut -f4 ${wd}/snps_combi_${disease}.bed | sort) > ${wd}"/ld/direct_candidates_"${disease}

# For each of the types of snps, find their respective entries in the GWAS Catalog associations
# note that for those ld proxies the snp that is featured and searched in the Catalog is the one in LD (column 6 of the {snp}.ld file), and not the SNP found by combi_nn itself.
cat ${folder}/header_GWASCat <(grep -wf ${wd}"/ld/direct_candidates_"${disease} ${wd}/${disease}.assoc) > ${wd}"/ld/direct_candidates_"${disease}.assoc
cat ${folder}/header_GWASCat <(grep -wf <(cut -f6 ${wd}"ld/ld_candidates_"${disease} | sort | uniq) ${wd}/${disease}.assoc) > ${wd}"/ld/ld_candidates_"${disease}.assoc

echo "All done!"
# END
