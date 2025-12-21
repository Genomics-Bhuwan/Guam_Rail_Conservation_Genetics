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

#Set variables
# ANGSD executable (full path)
ANGSD_EXE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/Guam_rail_Population_Genomics/angsd
REAL_SFS_EXE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/Guam_rail_Population_Genomics/angsd/misc/realSFS"  # adjust if different

# Directories and reference files
bam_dir="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled"
output_dir="${bam_dir}/heterozygosity"
ancestral_fasta="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/bHypOws1_hifiasm.bp.p_ctg.fasta"
reference_fasta="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/bHypOws1_hifiasm.bp.p_ctg.fasta"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# List of BAM files
samples=("SRR17129394.bam" "SRR17134085.bam" "SRR17134086.bam" "SRR17134087.bam" "SRR17134088.bam")

# Loop over each sample
for bam_file in "${samples[@]}"; do
    SAMPLE=$(basename "$bam_file" .bam)
    
    # Loop over autosomes 1-17 (matching FASTA headers)
    for i in {1..17}; do
        CHR="$i"
        echo "[$(date)] Processing $SAMPLE chromosome $CHR..."
        
        # Run ANGSD to calculate SAF
        $ANGSD_EXE -P 24 \
            -i "${bam_dir}/${bam_file}" \
            -anc "$ancestral_fasta" \
            -ref "$reference_fasta" \
            -dosaf 1 \
            -gl 1 \
            -C 50 \
            -minQ 20 \
            -minmapq 30 \
            -out "${output_dir}/${SAMPLE}.${CHR}" \
            -r "$CHR"
        
        # Estimate folded SFS
        $REAL_SFS_EXE -fold 1 "${output_dir}/${SAMPLE}.${CHR}.saf.idx" > "${output_dir}/${SAMPLE}.${CHR}.est.ml"
    done
done

echo "[$(date)] All ANGSD runs completed."
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
