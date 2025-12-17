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
