#### Calcualte the coverage of the samples.

```bash

module load samtools
module load parallel

# List of RGs
RG_LIST="FMNH390989 HOW_N23-0063 HOW_N23-0568"

# Parallel command to calculate average depth and genome coverage per RG
parallel -j 3 '
echo -n "{} "
samtools view -b -r {} rm_duplicates_BAM/Guam_Rail_merged.bam | \
samtools depth -aa - | \
awk '\''{sum+=$3; if($3>0) covered++; cnt++} END {if(cnt>0) printf "AvgDepth=%.2f, Coverage%%=%.2f\n", sum/cnt, covered/cnt*100; else print "0"}'\'' 
' ::: $RG_LIST
```
#### Result
```bash
HOW_N23-0568 28.0777
HOW_N23-0063 33.068
FMNH390989 58.4062
```
- *Step 1*: Decide target depth
- Pick the lowest depth among these three:
- I am picking the sample with the lowest depth and making other samples depth equivalent to that samples depth.
- HOW_N23-0568 → 28.08×


*Step 2*: Compute downsampling fractions
Sample	Depth	Fraction (F)
FMNH390989	58.41	28.08 / 58.41 ≈ 0.48
HOW_N23-0063	33.07	28.08 / 33.07 ≈ 0.85
HOW_N23-0568	28.08	28.08 / 28.08 = 1

*Step 3*: Downsample per read group using samtools view -s
mkdir -p rm_duplicates_BAM/downsampled

#### FMNH390989
```bash
#!/bin/bash

# Paths
BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/Guam_Rail_merged.bam"
OUTDIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/downsampled"
MERGED_BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam"

# Create output directory
mkdir -p $OUTDIR

# Downsample each sample
samtools view -b -s 0.48 -r FMNH390989 $BAM > $OUTDIR/FMNH390989_downsampled.bam
samtools view -b -s 0.85 -r HOW_N23-0063 $BAM > $OUTDIR/HOW_N23-0063_downsampled.bam
samtools view -b -r HOW_N23-0568 $BAM > $OUTDIR/HOW_N23-0568_downsampled.bam

# Merge downsampled BAMs
samtools merge $MERGED_BAM $OUTDIR/*.bam

# Index merged BAM
samtools index $MERGED_BAM

echo "Downsampling, merging, and indexing complete!"

samtools index rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam
```

✅ Result: BAM with roughly equal coverage (~28×) for all three samples.
