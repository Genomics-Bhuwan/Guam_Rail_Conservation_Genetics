#### Structural variants analysis using smoove pipeline
- https://github.com/brentp/smoove

#### Step 1 a. Download the smoove pipeline using wget
```bash
# Download the smoove binary
wget https://github.com/brentp/smoove/releases/download/v0.2.8/smoove

# Make it executable
chmod +x smoove

# Run the help command to see if it works
./smoove -h

```
#### Step 1 b. Install the dependencies
  ```bash
  mamba install -c bioconda -c conda-forge htslib gsort lumpy-sv samtools svtyper mosdepth
  ```

  #### Step 2. Running for small sample cohort for structural variants
```bash
  ./smoove call -x \
    --name Guam_Gazelle_Cohort \
    --fasta /shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/bHypOws1_hifiasm.bp.p_ctg.fasta \
    -p 20 \
    --genotype \
    /shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/*.bam
```
#### Step 3. Run the duphold to add the depth annotations(DHFFC).
# Optional: Run duphold to add depth annotations (DHFFC)
# Note: You would usually add -d to the initial 'call' to do this automatically
# If you didn't, you can run it now:
- Install the duphold for running.
---
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/smoove

wget https://github.com/brentp/duphold/releases/download/v0.2.1/duphold
chmod +x duphold
```

#### Step 3 a. Annotation included
```bash
./smoove duphold -x \
    --name Guam_rail_Annotated \
    --vcf Guam_Gazelle_Cohort-smoove.genotyped.vcf.gz \
    --fasta /shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/bHypOws1_hifiasm.bp.p_ctg.fasta \
    /shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/*.bam
```

#### Step 3 b. Run annotate with a GFF file for the Dama gazelle.
```bash
./smoove annotate --gff Guam_rail.gff Guam_Gazelle_Cohort-smoove.genotyped.vcf.gz | bgzip -c > Guam_rail_Final.vcf.gz
```

#### Step 4. Variant filtration for Smoove for high-quality genotypes and real depth changes using LUMPY
- Filter Smoove for high-quality genotypes and real depth changes
- MSHQ > 3: High quality heterozygotes
- DHFFC: deletions with DHFFC < 0.7 and duplications with DHFFC > 1.25
- These are stringently requested by the smoove pipeline.
```bash
bcftools view -i 'MSHQ > 3 && (SVTYPE == "DEL" && DHFFC < 0.7 || SVTYPE == "DUP" && DHFFC > 1.3 || SVTYPE == "INV")' \
    /Guam_Gazelle_Cohort-smoove.genotyped.vcf.gz | \
bcftools view -i 'QUAL >= 100' -O z -o Guam_rail_Smoove_Final_filtered.vcf.gz

# Always index after filtering
tabix -p vcf Dama_Gazelle_Smoove_Clean.vcf.gz
```



#### Step 5. Running Manta and Delly for structural variant calling
- DELLY- https://github.com/dellytools/delly
- singularity pull delly.sif docker://dellytools/delly:latest    : Use this for installing the delly.

#### Step 5 a. Preparation of the files.
- Create the exclusion file cause you have unplaced scaffolds and variable sex chromosomes.
- Avoid wasting time on hundreds of h1tg scaffolds.

```bash
#!/bin/bash
#SBATCH --job-name=delly_guam
#SBATCH --output=logs/delly_%a.out
#SBATCH --error=logs/delly_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=60:00:00
#SBATCH --array=0-2
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# --- 1. Absolute Paths ---
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/bHypOws1_hifiasm.bp.p_ctg.fasta"
BAM_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants"
OUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants"
DELLY_IMG="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/delly.sif"

# --- 2. Sample List ---
# These must match the start of your .bam filenames exactly
SAMPLES=("FMNH390989" "HOW_N23-0063" "HOW_N23-0568")
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# --- 3. Execution ---
mkdir -p logs

echo "Starting Delly for Sample: $SAMPLE"
echo "Using Reference: $REF"

# Note: I removed the -x (exclude) flag so you don't have to worry about chromosome names.
# Delly will process every contig found in the FASTA.
singularity exec --bind /shared:/shared "$DELLY_IMG" delly call \
    -g "$REF" \
    -o "${OUT_DIR}/${SAMPLE}.bcf" \
    "${BAM_DIR}/${SAMPLE}_downsampled.bam"

echo "Finished at $(date)"
```
#### Step 5 c.: Merge Sites
- We are now merging the findings from all 5 samples into a single unified list of candidate SVs.
```bash
singularity exec --bind /shared:/shared $DELLY_IMG delly merge -o sites.bcf \
    SRR17129394.bcf SRR17134085.bcf SRR17134086.bcf SRR17134087.bcf SRR17134088.bcf
```
#### Step 5 d.: Genotyping (Per Sample)
- We now go back to the original BAM files to see  which sample has which variant from the unified sites.bcf list.
```bash
for s in "${SAMPLES[@]}"; do
    echo "Genotyping: $s"
    singularity exec --bind /shared:/shared $DELLY_IMG delly call \
        -g $REF \
        -v sites.bcf \
        -x exclude.txt \
        -o "${s}.geno.bcf" \
        "${BAM_DIR}/${s}_mapped_sorted_RG_rmdup.bam"
done
```

#### Step 6 d.: Final Merge and Chromosome Filtering
- Finally, we merge the individual genotyped fiels and use bcftools to stringently keep only chromosomes 1-17.
```bash
# Merge all 5 samples into one multi-sample BCF
# Merge genotyped files
singularity exec --bind /shared:/shared $DELLY_IMG bcftools merge -m id -O b -o dama_merged.bcf \
    SRR17129394.geno.bcf SRR17134085.geno.bcf SRR17134086.geno.bcf SRR17134087.geno.bcf SRR17134088.geno.bcf

# Index the result
singularity exec --bind /shared:/shared $DELLY_IMG bcftools index dama_merged.bcf

# Filter for Chrs 1-17 and output VCF
singularity exec --bind /shared:/shared $DELLY_IMG bcftools view dama_merged.bcf \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 \
    -o dama_gazelle_final_shortread_1-17.vcf

```
- Summary of the above files created
- exclude.txt lists everuthing that you want to ignore.
- sites.bcf: master list of SVs found across the cohort of 5.
- *.geno.bcf: Individual files containing the presence/absence of SVs.
  
#### Step 6 e.: Final filtration
- Use bcftools and only keep the varaints that pass the quality. Rest eliminate.
- After filtration, look at the genotype field in the vcf file and see
- 0/0: Homozygous Reference (No SV).
- 0/1: Heterozygous SV.
- 1/1: Homozygous Alternative SV.
- USE IGV OR WALLY for genome visualtion for the SV: https://github.com/tobiasrausch/wally



##################################################################################################################################
#####  Step 7. WHAMg for structural variants calling##############################################################################
##################################################################################################################################
- Link to whamg pipeline: https://github.com/zeeev/wham

#### Step 7 a. Installation of the whamg
```bash
git clone --recursive  https://github.com/zeeev/wham.git; cd wham; make
```

#### Step 7 b. Runnig whamg
```bash
# 1. Define paths
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/Dama_gazelle_hifiasm-ULONT_primary.fasta"
OUT_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants/WHAM"
EXCLUDE_FILE="$OUT_DIR/exclude.txt"

# 2. Format the exclusion list (converts newlines to commas)
# This command takes the file, replaces newlines with commas, and removes the trailing comma.
EXCLUDE_LIST=$(tr '\n' ',' < "$EXCLUDE_FILE" | sed 's/,$//')

# 3. Define BAM files (comma-separated)
BAMS="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/SRR17129394_mapped_sorted_RG_rmdup.bam,\
/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/SRR17134085_mapped_sorted_RG_rmdup.bam,\
/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/SRR17134086_mapped_sorted_RG_rmdup.bam,\
/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/SRR17134087_mapped_sorted_RG_rmdup.bam,\
/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/SRR17134088_mapped_sorted_RG_rmdup.bam"

# 4. Run whamg
# -e: our newly formatted comma-separated exclusion list
# -a: reference genome
# -f: the list of BAMs
# -x: using 8 CPUs
# -z: force sampling (useful for de novo assemblies with many contigs)
whamg -e "$EXCLUDE_LIST" -a "$REF" -f "$BAMS" -x 8 -z \
    > "$OUT_DIR/Dama_Gazelle_joint_calls.vcf" \
    2> "$OUT_DIR/Dama_Gazelle_run.err"
 ```
- It is recommended that I look at the most interesting SVs in the genome browser like IGV(Integrative Genomics Viewer).
- What to look for: Look at your BAM files and the whamg vcf file. I want to see the clear cliff in coverage or a cluster of colored reads(discordant pairs) at the excat spots whamg made the call.

#### 7 C. Filtration of the noise or the artifacts
```bash
# 1. Define your file paths
INPUT_VCF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants/WHAM/Dama_Gazelle_joint_calls.vcf"
FINAL_VCF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants/WHAM/Dama_Gazelle_clean_final.vcf"

# 2. Run the combined filter
# EXPLANATION OF LOGIC:
# ABS(SVLEN)>=50       -> Keep only variants 50bp or larger (Standard SV definition).
# ABS(SVLEN)<=2000000  -> Remove variants larger than 2Mb (Usually assembly artifacts).
# INFO/A>=5            -> Total supporting reads across all 5 samples must be at least 5.
# INFO/CW[4]<0.2       -> The 5th value in the CW field (BND) must be less than 20%. 
#                         This removes false positives caused by repetitive sequences.

bcftools filter -i 'ABS(SVLEN)>=50 && ABS(SVLEN)<=2000000 && INFO/A>=5 && INFO/CW[4]<0.2' \
    $INPUT_VCF \
    -O v -o $FINAL_VCF

# 3. Quick Stats Comparison
echo "------------------------------------------------"
echo "Filtering Complete!"
echo "Raw Variants:      $(grep -v "^#" $INPUT_VCF | wc -l)"
echo "Filtered Variants: $(grep -v "^#" $FINAL_VCF | wc -l)"
echo "------------------------------------------------"
```
   
