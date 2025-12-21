# The Genome Analysis Toolkit (GATK) 
---
The GATK (Genome Analysis Toolkit) is one of the most used programs for genotype calling in sequencing data in model and non model organisms. 
Designed to analyze human genetic data and all its pipelines are optimized for the current purpose.
---

**Author:** Bhuwan Singh Bist

**Date:** 12/15/2025

###### Step 1. Haplotype Caller
###### Aligned BAM file(reads mapped to reference genome) and call all possible variants(SNPs and indels) from the samples.
###### Generates Genomic VCF(GVCF) which contains both variants and non-variants.
###### Since, my deduplicated bam file consists of three samples. I want to know the name of the samples.
- Now, I removed the deduplicated and also downsampled to ~28X for all the three samples.

#### Location of the deduplicated-downsampled sample.
```bash
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam
```

### Run the batch script 
### Run all samples in parallel (fast and efficient)
#### I would have run GenomicsdbImport rather than this if I had more than 20 samples or higher number of chromosomes in my dataset to optimize the resources and time.
```bash
#!/bin/bash -l
#SBATCH --job-name=HaplotypeCaller_parallel
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=bigmem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/HaplotypeCaller_%A.out
#SBATCH --error=logs/HaplotypeCaller_%A.err

# -------------------------------
# Load modules
# -------------------------------
module load gatk-4.1.2.0
module load samtools-1.22.1

# -------------------------------
# Input/output paths
# -------------------------------
BAM=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam
REFERENCE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/bHypOws1_hifiasm.bp.p_ctg.fasta
GVCF_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs

mkdir -p $GVCF_DIR

# -------------------------------
# Extract sample names
# -------------------------------
samples=$(samtools view -H $BAM | grep '^@RG' | awk '{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}' | sed 's/SM://')

# -------------------------------
# Run HaplotypeCaller for each sample in parallel
# -------------------------------
for sample in $samples; do
  echo " Launching HaplotypeCaller for: $sample"

  gatk --java-options "-Xmx16g" HaplotypeCaller \
      -R $REFERENCE \
      -I $BAM \
      -O $GVCF_DIR/${sample}.g.vcf.gz \
      -ERC GVCF \
      --sample-name $sample \
      --native-pair-hmm-threads 4 \
      > $GVCF_DIR/${sample}.log 2>&1 &

done

# Wait for all background processes to complete
wait

echo "✅ All parallel HaplotypeCaller jobs finished successfully."
```

#### Step 2. Combine Gentopyes.Joint Genotyping
```bash
#!/bin/bash -l
#SBATCH --job-name=CombineGVCFs
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=logs/CombineGVCFs_%A.out
#SBATCH --error=logs/CombineGVCFs_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# -------------------------------
# Load modules
# -------------------------------
module load samtools-1.22.1
# Make sure GATK 4.3 is in your PATH
export PATH=/shared/jezkovt_bistbs_shared/Guam_Rail/.../gatk-4.3.0.0:$PATH

# -------------------------------
# Define paths
# -------------------------------
WORKDIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs
REF=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/bHypOws1_hifiasm.bp.p_ctg.fasta
OUTDIR=$WORKDIR/Combined_GVCF

mkdir -p $OUTDIR

# -------------------------------
# Combine per-sample GVCFs
# -------------------------------
gatk --java-options "-Xmx120G" CombineGVCFs \
   -R $REF \
   --variant $WORKDIR/FMNH390989.g.vcf.gz \
   --variant $WORKDIR/HOW_N23-0063.g.vcf.gz \
   --variant $WORKDIR/HOW_N23-0568.g.vcf.gz \
   -O $OUTDIR/all_samples_combined.g.vcf.gz

# -------------------------------
# Completion message
# -------------------------------
echo "✅ CombineGVCFs job finished successfully at $(date)"

```

#### Step 3. Run Joint Genotyping in CombinedGVCFs.
```bash
gatk --java-options "-Xmx120G" GenotypeGVCFs \
   -R $REF \
   --variant /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/all_samples_combined.g.vcf.gz \
   -O /scratch/bistbs/GATK_Variant_Calling/Genotyped_VCF/all_samples_genotyped.vcf.gz
```

#### Step 4. Variant Filtration
```bash
module load gatk-4.1.2.0

RAW_VCF=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/all_samples_genotyped.vcf.gz

FILTERED_VCF=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/all_samples_genotyped_gatkfiltered.vcf.gz

gatk VariantFiltration \
  -V $RAW_VCF \
  -filter "DP < 1800" --filter-name "DP1800" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O $FILTERED_VCF

```

#### Step 5.A VCF Filtering for Population Genomics. We only keep SNPs for Population Genomics Analysis.
```bash
# Set variables
RAW_VCF="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/all_samples_genotyped_gatkfiltered.vcf.gz"
OUTDIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped"
OUTPREFIX="${OUTDIR}/Guam_rail_without_indels"

# Filter VCF
vcftools \
  --gzvcf $RAW_VCF \
  --minQ 30 \
  --remove-indels \
  --recode \
  --recode-INFO-all \
  --out $OUTPREFIX

# Generate individual missingness file
vcftools \
  --vcf ${OUTPREFIX}.recode.vcf \
  --missing-indv

```


# This will produce 'out.miss' file with % missing genotypes per individual
# Individuals with very high missingness can be removed in later filtering steps

# 1) Make biallelic SNP-only (keep your existing file unchanged)
bcftools view -v snps -m2 -M2 /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_without_indels.recode.vcf \
  -Oz -o /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_biallelic_snps.vcf.gz
tabix -p vcf /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_biallelic_snps.vcf.gz

```
#### Step 5.B Keep indels if you want to do SNpeff and VEP(Variant Effect Predictor).
```bash
INPUT_VCF="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/Guam_rail_without_indels.recode.vcf"
OUTDIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped"
OUT_VCF="${OUTDIR}/Guam_rail_biallelic_snps.vcf.gz"

# Filter for biallelic SNPs and compress
bcftools view -v snps -m2 -M2 $INPUT_VCF -Oz -o $OUT_VCF

# Index the new VCF
tabix -p vcf $OUT_VCF

```


