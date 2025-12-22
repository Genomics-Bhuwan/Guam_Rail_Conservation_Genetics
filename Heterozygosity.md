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
       REAL_SFS_EXE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity/angsd/misc/realSFS"
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/heterozygosity"

SAMPLES=("FMNH390989_downsampled" "HOW_N23-0063_downsampled" "HOW_N23-0568_downsampled")

for SAMPLE in "${SAMPLES[@]}"; do
    for i in {1..46}; do
        SCAF="SUPER_${i}"
        echo "Running realSFS for $SAMPLE scaffold $SCAF..."
        
```

- The next step is to add the sample name and the scaffold number for each line in our output.
- This will make our work easier when we want to plot our results.

```bash
#!/bin/bash

# ----------------------------
# Annotate realSFS output files with number of lines, sample name, and scaffold
# ----------------------------

output_dir="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity"

# List of BAM/sample names (without .bam)
samples=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

# Loop over each sample
for SAMPLE in "${samples[@]}"; do
    echo "[$(date)] Annotating $SAMPLE files..."

    # Loop over scaffolds 1 to 17
    for i in {1..17}; do
        input_file="${output_dir}/${SAMPLE}.${i}.est.ml"
        output_file="${output_dir}/${SAMPLE}.${i}.est.ml.annotated"

        if [ -f "$input_file" ]; then
            # Count number of lines
            num_lines=$(wc -l < "$input_file")

            # Annotate each line with line count, sample, and scaffold
            awk -v lines="$num_lines" -v sample="$SAMPLE" -v scaffold="$i" '{print lines, sample, scaffold, $0}' "$input_file" > "$output_file"

            # Replace original file with annotated version
            mv "$output_file" "$input_file"

            echo "[$(date)] Annotated $input_file"
        else
            echo "[$(date)] WARNING: File $input_file not found, skipping."
        fi
    done
done

echo "[$(date)] All annotation completed!"
Concatenate files

#!/bin/bash

# ----------------------------
# Concatenate all realSFS output files for a given sample
# ----------------------------

input_directory="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity"
output_directory="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity/concatenated"
mkdir -p "$output_directory"

# Output file
output_file="${output_directory}/all_samples_est_ml_concatenated.txt"

# Remove output file if it already exists
[ -f "$output_file" ] && rm "$output_file"

# List of samples
samples=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

# Loop over each sample and scaffold 1-17
for SAMPLE in "${samples[@]}"; do
    for i in {1..17}; do
        input_file="${input_directory}/${SAMPLE}.${i}.est.ml"
        if [ -f "$input_file" ]; then
            cat "$input_file" >> "$output_file"
        else
            echo "WARNING: $input_file not found, skipping."
        fi
    done
done

echo "All files concatenated into $output_file"
