#### Estimating the nucleotide diversity of the rallidae.
#### I am attesting the codes used for the analysis for few samples. It will be unbearable to paste code for all samples in terms of space.
- a. Black rail
```bash
#!/bin/bash
#SBATCH --job-name=BlackRail_Theta
#SBATCH --output=BlackRail_Theta_%j.out
#SBATCH --error=BlackRail_Theta_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# Load required modules
module load gcc-6.3.0
module load angsd

# Add realSFS/thetaStat to PATH
export PATH=$PATH:/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc

# --- 1. SET DIRECTORIES (Updated for Black Rail) ---
# Reference Genome for Black Rail
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Black_rail/trimmed_fastq/GCA_022605575.1_bLatJam1.0.p_genomic.fna"

# Specific Black Rail BAM path
BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Black_rail/alignment/SRR18439748_rmdup.bam"

# Output Directories
BASE_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Nucleotide_diversity/Black_rail"
FINAL_OUT="$BASE_DIR/Final_output"

# Create directories if they don't exist
mkdir -p $BASE_DIR
mkdir -p $FINAL_OUT

# Define Sample ID based on the BAM file
SAMPLE="SRR18439748"

# --- 2. RUN ANGSD (Generate SAF) ---
echo "Starting ANGSD SAF generation for Black Rail ($SAMPLE)..."

angsd \
    -i $BAM \
    -anc $REF \
    -ref $REF \
    -out $BASE_DIR/$SAMPLE \
    -nThreads 10 \
    -GL 2 \
    -doSaf 1 \
    -minMapQ 30 \
    -minQ 30 \
    -baq 1 \
    -C 50 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -minInd 1 \
    -setMinDepth 5 \
    -setMaxDepth 80 \
    -doCounts 1

# --- 3. CALCULATE SFS AND THETAS ---
SAF_FILE="$BASE_DIR/${SAMPLE}.saf.idx"

if [ -f "$SAF_FILE" ]; then
    echo "Processing SFS for $SAMPLE ..."

    # Step A: Calculate Site Frequency Spectrum (SFS)
    # Note: For a single diploid individual, the SFS fold is usually size 2
    realSFS $SAF_FILE -P 10 > $FINAL_OUT/${SAMPLE}.sfs

    # Step B: Estimate thetas in binary format
    realSFS saf2theta $SAF_FILE \
        -sfs $FINAL_OUT/${SAMPLE}.sfs \
        -outname $FINAL_OUT/${SAMPLE}

    # Step C: Print per-site diversity
    thetaStat print $FINAL_OUT/${SAMPLE}.thetas.idx > $FINAL_OUT/${SAMPLE}_persite.theta.txt
    
    # Step D: Calculate windowed thetas (50kb window, 10kb step)
    thetaStat do_stat $FINAL_OUT/${SAMPLE}.thetas.idx \
        -win 50000 -step 10000 \
        -outnames $FINAL_OUT/${SAMPLE}.windowed.theta

    echo "Finished processing Black Rail."
else
    echo "Error: $SAF_FILE not found. Check ANGSD logs for errors during SAF generation."
fi
```
- b. Guam rail
  ```bash
  #!/bin/bash
#SBATCH --job-name=GuamRail_Thetas
#SBATCH --output=GuamRail_Thetas_%j.out
#SBATCH --error=GuamRail_Thetas_%j.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# Load required modules
module load gcc-6.3.0
module load angsd

# Add realSFS/thetaStat to PATH
export PATH=$PATH:/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc

# --- 1. SET DIRECTORIES ---
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/bHypOws1_hifiasm.bp.p_ctg.fasta"

# Directory where your script and SAF files will live
BASE_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Nucleotide_diversity/Guam_rail"
FINAL_OUT="$BASE_DIR/Final_output"

mkdir -p $FINAL_OUT

# List of your three BAM files
SAMPLES=("FMNH390989_downsampled" "HOW_N23-0063_downsampled" "HOW_N23-0568_downsampled")
BAM_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled"

# --- 2. LOOP THROUGH SAMPLES ---
for SAMPLE in "${SAMPLES[@]}"; do
    echo "--------------------------------------------"
    echo "Processing Sample: $SAMPLE"
    echo "--------------------------------------------"

    # Step A: Run ANGSD to generate SAF
    angsd \
        -i $BAM_DIR/${SAMPLE}.bam \
        -anc $REF \
        -ref $REF \
        -out $BASE_DIR/${SAMPLE} \
        -nThreads 10 \
        -GL 2 \
        -doSaf 1 \
        -minMapQ 30 \
        -minQ 30 \
        -baq 1 \
        -C 50 \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -minInd 1 \
        -setMinDepth 5 \
        -setMaxDepth 80 \
        -doCounts 1

    # Define SAF index path
    SAF_FILE="$BASE_DIR/${SAMPLE}.saf.idx"

    if [ -f "$SAF_FILE" ]; then
        # Step B: Calculate SFS
        realSFS $SAF_FILE -P 10 > $FINAL_OUT/${SAMPLE}.sfs

        # Step C: Estimate thetas
        realSFS saf2theta $SAF_FILE \
            -sfs $FINAL_OUT/${SAMPLE}.sfs \
            -outname $FINAL_OUT/${SAMPLE}

        # Step D: Print per-site diversity (Warning: Very large file)
        thetaStat print $FINAL_OUT/${SAMPLE}.thetas.idx > $FINAL_OUT/${SAMPLE}_persite.theta.txt
        
        # Step E: Windowed diversity (50kb window, 10kb step)
        thetaStat do_stat $FINAL_OUT/${SAMPLE}.thetas.idx \
            -win 50000 -step 10000 \
            -outnames $FINAL_OUT/${SAMPLE}.windowed.theta

        echo "Finished $SAMPLE"
    else
        echo "Error: SAF file for $SAMPLE not created."
    fi
done

echo "All Guam Rail samples processed."
```
- c. Okinawa rail
```bash
#!/bin/bash
#SBATCH --job-name=OkinawaRail_Theta
#SBATCH --output=OkinawaRail_Theta_%j.out
#SBATCH --error=OkinawaRail_Theta_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# Load required modules
module load gcc-6.3.0
module load angsd

# Add realSFS/thetaStat to PATH
export PATH=$PATH:/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc

# --- 1. SET DIRECTORIES ---
# Reference for Okinawa Rail
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Okinawa_rail/2_Adapter_trimming/GCA_027925045.1_Gokinawae_1.0_genomic.fna"

# Specific Okinawa Rail BAM path
BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Okinawa_rail/05_Remove_duplicates/DRR424032_rmdup.bam"

# Output Directories
BASE_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Nucleotide_diversity/Okinawa_rail"
FINAL_OUT="$BASE_DIR/Final_output"

mkdir -p $FINAL_OUT

# --- 2. RUN ANGSD (Generate SAF) ---
echo "Starting ANGSD SAF generation for Okinawa Rail..."

angsd \
    -i $BAM \
    -anc $REF \
    -ref $REF \
    -out $BASE_DIR/DRR424032 \
    -nThreads 10 \
    -GL 2 \
    -doSaf 1 \
    -minMapQ 30 \
    -minQ 30 \
    -baq 1 \
    -C 50 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -minInd 1 \
    -setMinDepth 5 \
    -setMaxDepth 80 \
    -doCounts 1

# --- 3. CALCULATE SFS AND THETAS ---
SAF_FILE="$BASE_DIR/DRR424032.saf.idx"
SAMPLE="DRR424032"

if [ -f "$SAF_FILE" ]; then
    echo "Processing SFS for $SAMPLE ..."

    # Step A: Calculate Site Frequency Spectrum (SFS)
    realSFS $SAF_FILE -P 10 > $FINAL_OUT/${SAMPLE}.sfs

    # Step B: Estimate thetas in binary format
    realSFS saf2theta $SAF_FILE \
        -sfs $FINAL_OUT/${SAMPLE}.sfs \
        -outname $FINAL_OUT/${SAMPLE}

    # Step C: Print per-site diversity
    thetaStat print $FINAL_OUT/${SAMPLE}.thetas.idx > $FINAL_OUT/${SAMPLE}_persite.theta.txt
    
    # Step D: Calculate windowed thetas (50kb window, 10kb step)
    thetaStat do_stat $FINAL_OUT/${SAMPLE}.thetas.idx \
        -win 50000 -step 10000 \
        -outnames $FINAL_OUT/${SAMPLE}.windowed.theta

    echo "Finished processing Okinawa Rail."
else
    echo "Error: $SAF_FILE not found. Check ANGSD logs."
fi
```
