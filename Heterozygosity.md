#### Estimating the Heterozygosity
- I am running ANGSD and SFS for Guam rail in three samples as given in code below.
```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_Het
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=angsd_het_%A.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# ----------------------------
# Paths
# ----------------------------
# ----------------------------
# Paths
# ----------------------------
ANGSD_EXE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/angsd/angsd"
BAM_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled"
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity"
REFERENCE="${BAM_DIR}/bHypOws1_hifiasm.bp.p_ctg.fasta"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# List of BAMs
SAMPLES=("FMNH390989_downsampled.bam" "HOW_N23-0063_downsampled.bam" "HOW_N23-0568_downsampled.bam")

# Get scaffold names from the reference FASTA
SCAFFOLDS=($(grep ">" "$REFERENCE" | sed 's/>//'))

# ----------------------------
# Run ANGSD scaffold by scaffold
# ----------------------------
for SAMPLE_BAM in "${SAMPLES[@]}"; do
    SAMPLE=$(basename "$SAMPLE_BAM" .bam)
    BAM_PATH="${BAM_DIR}/${SAMPLE_BAM}"

    for SCAF in "${SCAFFOLDS[@]}"; do
        echo "[$(date)] Running ANGSD for $SAMPLE scaffold $SCAF..."
        
        $ANGSD_EXE -P 24 \
            -i "$BAM_PATH" \
            -anc "$REFERENCE" \
            -ref "$REFERENCE" \
            -dosaf 1 \
            -gl 1 \
            -C 50 \
            -minQ 20 \
            -minmapq 30 \
            -out "${OUTPUT_DIR}/${SAMPLE}.${SCAF}" \
            -r "$SCAF"
    done
done

echo "[$(date)] All ANGSD runs completed."

```
#### Running Folded Site Frequency Spectrum
     ```bash
      REAL_SFS_EXE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/angsd/misc/realSFS"

INPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity"

OUTPUT_DIR="${INPUT_DIR}/realSFS_new"
mkdir -p "$OUTPUT_DIR"

SAMPLES=(
  "FMNH390989_downsampled"
  "HOW_N23-0063_downsampled"
  "HOW_N23-0568_downsampled"
)

for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing sample: $SAMPLE"

  for SAFIDX in ${INPUT_DIR}/${SAMPLE}.ptg*.saf.idx; do
    BASENAME=$(basename "$SAFIDX" .saf.idx)

    # Each scaffold in its own folder inside the new output directory
    RUN_DIR="${OUTPUT_DIR}/${SAMPLE}/${BASENAME}"
    mkdir -p "$RUN_DIR"
    cd "$RUN_DIR" || exit 1

    echo "  Running realSFS on $BASENAME"

    $REAL_SFS_EXE \
      "$SAFIDX" \
      -P 8 \
      -maxIter 100 \
      > "${BASENAME}.sfs"

  done
done
```
        

- The next step is to add the sample name and the scaffold number for each line in our output.
- This will make our work easier when we want to plot our results.

```bash

REAL_SFS_EXE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/angsd/misc/realSFS"

INPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity"
OUTPUT_DIR="${INPUT_DIR}/realSFS_new"
mkdir -p "$OUTPUT_DIR"

SAMPLES=(
  "FMNH390989_downsampled"
  "HOW_N23-0063_downsampled"
  "HOW_N23-0568_downsampled"
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"

    # Create sample-specific output directory
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUT"

    # Loop over SAF index files
    for SAFIDX in ${INPUT_DIR}/${SAMPLE}.ptg*.saf.idx; do
        # Skip empty SAF files
        if [ ! -s "$SAFIDX" ]; then
            echo "  Skipping empty SAF file: $SAFIDX"
            continue
        fi

        BASENAME=$(basename "$SAFIDX" .saf.idx)
        RUN_DIR="${SAMPLE_OUT}/${BASENAME}"
        mkdir -p "$RUN_DIR"
        cd "$RUN_DIR" || exit 1

        echo "  Running realSFS on $BASENAME"

        # Produce .est.ml file
        $REAL_SFS_EXE \
            "$SAFIDX" \
            -P 8 \
            -maxIter 100 \
            > "${BASENAME}.est.ml" 2> "${BASENAME}.log"

        echo "  Finished: ${BASENAME}.est.ml"
    done
done

echo "All realSFS runs completed."

```

#### Concatenate all realSFS output files for a given sample.
# Base directories
```bash
#!/bin/bash

# Base output directory containing the realSFS results
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/realSFS_new"

# Final concatenated output file
OUTPUT_FILE="${OUTPUT_DIR}/all_samples_est_ml_concatenated.txt"

# Remove output file if it already exists
[ -f "$OUTPUT_FILE" ] && rm "$OUTPUT_FILE"

# List of samples
SAMPLES=(
  "FMNH390989_downsampled"
  "HOW_N23-0063_downsampled"
  "HOW_N23-0568_downsampled"
)

# Loop over samples and their scaffold directories
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"

    for SCAF_DIR in "${SAMPLE_DIR}"/*; do
        # Ensure it is a directory
        if [ -d "$SCAF_DIR" ]; then
            # Expected est.ml file
            input_file="${SCAF_DIR}/$(basename "$SCAF_DIR").est.ml"

            if [ -f "$input_file" ]; then
                cat "$input_file" >> "$OUTPUT_FILE"
                echo "[$(date)] Added $input_file to concatenated file"
            else
                echo "[$(date)] WARNING: File $input_file not found, skipping."
            fi
        fi
    done
done

echo "[$(date)] Concatenation step completed!"
echo "Concatenated file: $OUTPUT_FILE"

```
