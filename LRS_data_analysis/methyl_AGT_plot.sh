#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=24:0:0
#$ -j y
#$ -N SNP_DNAm_plots
#$ -o /data/scratch/bt24076/
#$ -m beas


IN_DIR=""
OUT_DIR=""
REF="chr1.fa"
CPG="cpg_site.bed"
rs_id="rs41317993"
SNP="chr1:207788387"
cpg_site="cg10501210"
GT_FILE="/${rs_id}_genotype.tsv"
HP_FILE="/${rs_id}_haplotype_assigned.tsv"
PLOT=""



## STEP 1 ##
#load conda environment
module load miniforge/
mamba activate modkit

#modkit for methylation from bam files
for bam in "$IN_DIR"/*.bam; do
  sample_filename=$(basename "$bam")
  sample_id=$(echo "$sample_filename" | cut -d'-' -f1)

  #pileup
  modkit pileup --ref "$REF" \
    --cpg \
    -t "${NSLOTS}" \
    --partition-tag HP \
    --include-bed "$CPG" \
    --ignore h \
    --combine-strands \
    "$bam" \
    "$OUT_DIR" \
    --log-filepath "$OUT_DIR/${sample_id}-pileup.log" 
done

#deactivate modkit and environment
Deactivate modkit env
mamba deactivate
module unload miniforge/24.7.1 
module purge 

##STEP 2 ##
#load bcftools for SNP genotype
module load bcftools

#Create header for combined gt tsv file
echo -e "Sample_ID\tGenotype\tREF\tALT" > "$GT_FILE"

# Get genotype of SNP
for bam in "$IN_DIR"/*.bam; do
    filename=$(basename "$bam")
    sample=$(echo "$filename" | cut -d'-' -f1)  

    #call genotype and put into vcf file
    bcftools mpileup -f "$REF" -r "$SNP" "$bam" | \
    bcftools call -m -Ov -V indels -o "$OUT_DIR/${sample}_GT.vcf"

    # Extract genotype
    GT=$(bcftools query -f '[%GT]\n' "$OUT_DIR/${sample}_GT.vcf")
    REF_ale=$(bcftools query -f '[%REF]\n' "$OUT_DIR/${sample}_GT.vcf")
    ALT_ale=$(bcftools query -f '[%ALT]\n' "$OUT_DIR/${sample}_GT.vcf")

    #If no snp at location then must not have been called in bam
    if [[ -z "$GT" ]]; then
        GT="NA"
        REF_ale="NA"
        ALT_ale="NA"
    fi
    
    #add sample ID and GT into the tsv
    echo -e "${sample}\t${GT}\t${REF_ale}\t${ALT_ale}" >> "$GT_FILE"
    #remove the vcf 
    rm "$OUT_DIR/${sample}_GT.vcf"

done

## STEP 3 ##
#Extract hp information from phased vcfs
#Create header for combined HP tsv file
echo -e "Sample_ID\tHaplotype\tAllele" > "$HP_FILE"

for vcf in "$IN_DIR"/*.vcf.gz; do
    filename=$(basename "$vcf")
    sample=$(echo "$filename" | cut -d'-' -f1)  

    # Extract GT, REF, ALT for a single SNP position
    GT=$(bcftools query -r "$SNP" -f '[%GT]\n' "$vcf")
    REF_ale=$(bcftools query -r "$SNP" -f '[%REF]\n' "$vcf")
    ALT_ale=$(bcftools query -r "$SNP" -f '[%ALT]\n' "$vcf")

    # If SNP not found in this VCF, mark both haplotypes as NA
    if [[ -z "$GT" ]]; then
        echo -e "${sample}\t1\tNA" >> "$HP_FILE"
        echo -e "${sample}\t2\tNA" >> "$HP_FILE"
        continue
    fi

    # Assign alleles based on phased genotype
    if [[ "$GT" == "0|1" ]]; then
        echo -e "${sample}\t1\t${REF_ale}" >> "$HP_FILE"
        echo -e "${sample}\t2\t${ALT_ale}" >> "$HP_FILE"
    elif [[ "$GT" == "1|0" ]]; then
        echo -e "${sample}\t1\t${ALT_ale}" >> "$HP_FILE"
        echo -e "${sample}\t2\t${REF_ale}" >> "$HP_FILE"
    elif [[ "$GT" == "1|1" ]]; then
        echo -e "${sample}\t1\t${ALT_ale}" >> "$HP_FILE"
        echo -e "${sample}\t2\t${ALT_ale}" >> "$HP_FILE"
    elif [[ "$GT" == "0|0" ]]; then
        echo -e "${sample}\t1\t${REF_ale}" >> "$HP_FILE"
        echo -e "${sample}\t2\t${REF_ale}" >> "$HP_FILE"
    else
        echo -e "${sample}\t1\tNA" >> "$HP_FILE"
        echo -e "${sample}\t2\tNA" >> "$HP_FILE"
    fi
done

# #unload bcftool 
module unload bcftools

## STEP 4 ##
#load R and install libraries
module load R
Rscript /data/home/bt24076/scripts/00_Rlibs.R 

#and combine bed files to tsv and calulate average methylation between both haplotypes for each sample
Rscript /data/home/bt24076/scripts/01_Combine_beds.R  $OUT_DIR $PLOT 

#then comnine genotype, haplotype and methylation data for plotting 
Rscript /data/home/bt24076/scripts/02_Combine_Data.R $PLOT $rs_id

#finally plot results and send to plot directory-(genotype)
Rscript /data/home/bt24076/scripts/03_LRS_plots.R $PLOT "$PLOT/combined_methylation_genotype.tsv" $rs_id $cpg_site

module unload R 