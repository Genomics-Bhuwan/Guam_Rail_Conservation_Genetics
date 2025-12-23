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

# Directory containing realSFS output
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/realSFS_new"

# List of samples
SAMPLES=(
  "FMNH390989_downsampled"
  "HOW_N23-0063_downsampled"
  "HOW_N23-0568_downsampled"
)

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
    echo "Processing sample: $SAMPLE"

    # Loop over scaffold directories for this sample
    for SCAF_DIR in "${SAMPLE_DIR}"/*; do
        if [ -d "$SCAF_DIR" ]; then
            input_file="${SCAF_DIR}/$(basename "$SCAF_DIR").est.ml"
            output_file="${SCAF_DIR}/$(basename "$SCAF_DIR").est.ml.annotated"

            if [ -f "$input_file" ]; then
                # Count number of lines
                num_lines=$(wc -l < "$input_file")

                # Annotate each line with line count, sample, and scaffold
                awk -v lines="$num_lines" -v sample="$SAMPLE" -v scaffold="$(basename "$SCAF_DIR")" \
                    '{print lines, sample, scaffold, $0}' "$input_file" > "$output_file"

                # Replace original file with annotated version
                mv "$output_file" "$input_file"

                echo "[$(date)] Annotated $input_file"
            else
                echo "[$(date)] WARNING: File $input_file not found, skipping."
            fi
        fi
    done
done

echo "[$(date)] All annotation completed!"


```

#### Concatenate all realSFS output files for a given sample.
# Base directories
```bash
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/realSFS_new"
OUTPUT_FILE="${OUTPUT_DIR}/all_samples_est_ml_concatenated.txt"

# Remove output file if it already exists
[ -f "$OUTPUT_FILE" ] && rm "$OUTPUT_FILE"

# List of samples
SAMPLES=(
  "FMNH390989_downsampled"
  "HOW_N23-0063_downsampled"
  "HOW_N23-0568_downsampled"
)

# Loop over samples and scaffold directories
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
    for SCAF_DIR in "${SAMPLE_DIR}"/*; do
        if [ -d "$SCAF_DIR" ]; then
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
