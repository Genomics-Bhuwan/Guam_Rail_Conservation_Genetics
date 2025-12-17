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
# -------------------------------
# Load GATK module
# -------------------------------
module load gatk-4.1.2.0

# -------------------------------
# Input/output paths
# -------------------------------
RAW_VCF=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/all_samples_genotyped.vcf.gz
FILTERED_VCF=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/all_samples_genotyped_gatkfiltered.vcf.gz

# -------------------------------
# Run VariantFiltration
# -------------------------------
gatk VariantFiltration \
   -V $RAW_VCF \
   
   # -------------------------------
   # Filter sites with extremely high depth
   # DP < 1800: remove sites where coverage is unusually high
   # High depth can indicate collapsed repeats or duplicated regions
   # -------------------------------
   -filter "DP < 1800" --filter-name "DP1800" \
   
   # -------------------------------
   # Strand bias filters
   # FS > 60.0: Fisher Strand test; removes sites with strong strand bias
   # SOR > 3.0: Strand Odds Ratio; another measure of strand bias
   # Sites failing these filters are likely sequencing/mapping artifacts
   # -------------------------------
   -filter "FS > 60.0" --filter-name "FS60" \
   -filter "SOR > 3.0" --filter-name "SOR3" \
   
   # -------------------------------
   # Mapping quality filter
   # MQ < 40.0: remove sites where reads are poorly aligned
   # Poorly mapped reads can introduce false variants
   # -------------------------------
   -filter "MQ < 40.0" --filter-name "MQ40" \
   
   # -------------------------------
   # Mapping quality rank sum test
   # MQRankSum < -12.5: tests if alternate alleles have lower mapping quality than reference
   # Extreme negative values indicate alignment bias against alternate alleles
   # -------------------------------
   -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
   
   # -------------------------------
   # Read position rank sum test
   # ReadPosRankSum < -8.0: tests if alternate alleles occur at biased positions in reads
   # Negative values indicate alt alleles mostly at read ends, often errors
   # -------------------------------
   -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
   
   # -------------------------------
   # Output filtered VCF
   # -------------------------------
   -O $FILTERED_VCF
```

#### Step 5.A VCF Filtering for Population Genomics. We only keep SNPs for Population Genomics Analysis.
```bash
module load vcftools

RAW_VCF=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/all_samples_genotyped_gatkfiltered.vcf.gz
VCFTOOLS_OUT=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_without_indels

# ------------------------------
# Filter with VCFtools
# ------------------------------
vcftools --gzvcf $RAW_VCF \                # Input VCF (can be gzipped)
         --minQ 30 \                       # Keep only sites with minimum quality score of 30 (high-confidence genotypes)
         --remove-indels \                 # Remove INDELs, keeping only SNPs
         --recode --recode-INFO-all \      # Produce a new VCF while keeping all INFO fields
         --out $VCFTOOLS_OUT               # Output prefix (final file: Dama_gazelle_without_indels.recode.vcf)

# ------------------------------
# Check missingness per individual
# ------------------------------
vcftools --vcf ${VCFTOOLS_OUT}.recode.vcf \  # Use the SNP-only VCF as input
         --missing-indv                     # Reports fraction of missing genotypes per sample


# This will produce 'out.miss' file with % missing genotypes per individual
# Individuals with very high missingness can be removed in later filtering steps

# 1) Make biallelic SNP-only (keep your existing file unchanged)
bcftools view -v snps -m2 -M2 /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_without_indels.recode.vcf \
  -Oz -o /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_biallelic_snps.vcf.gz
tabix -p vcf /scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/Without_Indels/Dama_gazelle_biallelic_snps.vcf.gz

```
#### Step 5.B Keep indels if you want to do SNpeff and VEP(Variant Effect Predictor).
```bash
module load vcf-tools

RAW_VCF=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/all_samples_genotyped_gatkfiltered.vcf.gz
VCFTOOLS_OUT=/scratch/bistbs/GATK_Variant_Calling/Combined_GVCF/Genotyped_VCF/With_Indels/Dama_gazelle_with_indels

# ------------------------------
# Filter with VCFtools (keep SNPs + INDELs)
# ------------------------------
vcftools \
  --gzvcf $RAW_VCF \
  --minQ 30 \
  --recode --recode-INFO-all \
  --out $VCFTOOLS_OUT

# ------------------------------
# Check missingness per individual
# ------------------------------
vcftools \
  --vcf ${VCFTOOLS_OUT}.recode.vcf \
  --missing-indv

```


