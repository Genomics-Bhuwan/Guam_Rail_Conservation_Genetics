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
# 1. Fix the image path variable
export DELLY_IMG="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/DELLY/delly.sif"

# 2. Run the merge using your actual filenames
singularity exec --bind /shared:/shared $DELLY_IMG delly merge -o sites.bcf \
FMNH390989.bcf \
HOW_N23-0063.bcf \
HOW_N23-0568.bcf
```
#### Step 5 d.: Genotyping (Per Sample)
- We now go back to the original BAM files to see  which sample has which variant from the unified sites.bcf list.
```bash
# 1. Define paths and image
export DELLY_IMG="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/DELLY/delly.sif"
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/bHypOws1_hifiasm.bp.p_ctg.fasta"
BAM_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants"

# 2. Define the sample IDs
SAMPLES=("FMNH390989" "HOW_N23-0063" "HOW_N23-0568")

# 3. Run the Genotyping loop
for s in "${SAMPLES[@]}"; do
    echo "=========================================="
    echo "Genotyping Sample: $s"
    echo "=========================================="
    
    singularity exec --bind /shared:/shared $DELLY_IMG delly call \
        -g "$REF" \
        -v sites.bcf \
        -o "${s}.geno.bcf" \
        "${BAM_DIR}/${s}_downsampled.bam"
done
```

#### Step 6 d.: Final Merge and Chromosome Filtering
- Finally, we merge the individual genotyped fiels and use bcftools to stringently keep only chromosomes 1-17.
```bash
# 1. Merge the genotyped Guam Rail files
# Running bcftools directly since the module is loaded
bcftools merge \
    -m id \
    -O b \
    -o guam_rail_merged.bcf \
    FMNH390989.geno.bcf HOW_N23-0063.geno.bcf HOW_N23-0568.geno.bcf

# 2. Index the merged result
bcftools index guam_rail_merged.bcf

# 3. Output the final VCF (All contigs included)
bcftools view \
    guam_rail_merged.bcf \
    -o guam_rail_final_all_contigs.vcf

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
```bash
# Filter for PASS variants and output to VCF
bcftools view -f PASS guam_rail_merged.bcf -o guam_rail_final_filtered.vcf
```
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
#!/bin/bash
#SBATCH --job-name=wham_guam
#SBATCH --output=wham_%j.log
#SBATCH --error=wham_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=70:00:00
#SBATCH --partition=batch
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# --- 1. Fix TMPDIR to avoid Slurm errors ---
MY_TEMP_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/tmp_wham"
mkdir -p "$MY_TEMP_DIR"
export TMPDIR="$MY_TEMP_DIR"

# --- 2. Paths ---
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/bHypOws1_hifiasm.bp.p_ctg.fasta"
OUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants"
WHAMG_EXE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants/WHAM/wham/bin/whamg"

# --- 3. Guam Rail BAMs (Joint Calling) ---
# Ensure there are NO spaces between the filenames, only commas.
BAMS="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/FMNH390989_downsampled.bam,\
/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/HOW_N23-0063_downsampled.bam,\
/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/HOW_N23-0568_downsampled.bam"

# --- 4. Execution ---
mkdir -p "$OUT_DIR"

echo "Starting WHAM-G Joint Calling for Guam Rail"
echo "Job ID: $SLURM_JOB_ID"
echo "BAM files being processed: $BAMS"

# Note: The -e flag has been removed since there is no exclude list.
$WHAMG_EXE \
    -a "$REF" \
    -f "$BAMS" \
    -x 20 \
    -z > "$OUT_DIR/Guam_Rail_joint_calls.vcf"

echo "Finished at $(date)"
 ```
- It is recommended that I look at the most interesting SVs in the genome browser like IGV(Integrative Genomics Viewer).
- What to look for: Look at your BAM files and the whamg vcf file. I want to see the clear cliff in coverage or a cluster of colored reads(discordant pairs) at the excat spots whamg made the call.

#### 7 C. Filtration of the noise or the artifacts
```bash
# 1. Define your file paths
# 1. Define paths for Guam Rail
INPUT_VCF="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/WHAM/Guam_Rail_joint_calls.vcf"
FINAL_VCF="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/WHAM/Guam_Rail_clean_final.vcf"

# 2. Run the combined filter
# Logic: 
# - Size between 50bp and 2Mb
# - At least 5 supporting reads total (INFO/A)
# - CW[4]<0.2 to filter out BND artifacts in repeats
bcftools filter -i 'ABS(SVLEN)>=50 && ABS(SVLEN)<=2000000 && INFO/A>=5 && INFO/CW[4]<0.2' \
    $INPUT_VCF \
    -O v -o $FINAL_VCF

# 3. Quick Stats Comparison
echo "------------------------------------------------"
echo "Guam Rail WHAM Filtering Complete!"
echo "Raw Variants:      $(grep -v "^#" $INPUT_VCF | wc -l)"
echo "Filtered Variants: $(grep -v "^#" $FINAL_VCF | wc -l)"
echo "------------------------------------------------"----------"
```


   #### 7. Use SURVIVOR FOR merging the three vcf files and generate the consensus vcf file.
#### Step 7.a Installation of the survivor.

   ```bash
cd /shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/SURVIVOR

# Download and unpack
wget https://github.com/fritzsedlazeck/SURVIVOR/archive/master.tar.gz -O SURVIVOR.tar.gz
tar xzvf SURVIVOR.tar.gz

# Compile
cd SURVIVOR-master/Debug/
make

# Move the executable up to your main SURVIVOR directory for easier access
cp SURVIVOR ../../
cd ../../

```

#### Step 7. b. Create the sample files text file that tells SURVIVOR which vcfs to compare.
- Create the list of the samples
```bash
SMOOVE="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/smoove/FINAL/Guam_Gazelle_Cohort-smoove.genotyped.vcf.gz"
DELLY="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/DELLY/Final_Filtration/delly_final.vcf.gz"
WHAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Structural_Variants/WHAM/Final_Filtered_File/wham_final.vcf.gz"

# Create the input list
printf "$DELLY\n$WHAM\n$SMOOVE\n" > guam_rail_vcf_list.txt
```
#### Step 7 c: Run the Consensus MergeNow, execute the merge using the parameters defined in your documentation. 
- We will use:1000: 1kb breakpoint distance.2: Minimum support (at least 2 callers must agree).1: Agree on SV Type.1: Agree on Strand.0: No size-specific distance.50: Minimum SV length of 50bp.
```bash
./SURVIVOR merge guam_rail_vcf_list.txt 1000 2 1 1 0 50 Guam_Rail_Consensus_Final.vcf
```
#### Step 7 d. Optional: Draw  Venn Diagram. Since I have three callers and I likely want to know the overlap.
- i.e., how many were found by all three vs. just two.
- You can extract the support vector information for R as mentioned. 
```bash
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' Guam_Rail_Consensus_Final.vcf | sed -e 's/\(.\)/\1 /g' > guam_rail_overlap_stats.txt
```

#### Step 7 e. Summary of the Final Guam_Rail_Consensus_Final.vcf.
- This is the gold standard SV set.
- Every variant in this file has been double checked  by at least two different algorithms.
- Us this in R to create a venn diagram for your publication to show the consistency between SMOOVE, DELLY AND WHAMg.
- To see, how many high-confidence SVs survived the merging.
  
```bash
grep -v "^#" Guam_Rail_Consensus_Final.vcf | wc -l
```
