# Dama Gazelle Whole-Genome Resequencing (WGS) Variant Calling Pipeline Tutorial

**Author:** Bhuwan Singh Bist

**Affiliation:** Jezkova lab

**Date:** 2025-10-03

This tutorial demonstrates a complete workflow for WGS variant calling in the Dama gazelle. The pipeline includes:

1. Quality control (FastQC)
2. Adapter trimming (Trim Galore + Cutadapt)
3. Mapping reads to a haplotype-resolved Ruminant Telomere-to-Telomere(T2T)reference genome assembly of Dama gazelle(Nanger dama)
4. BAM processing (Picard & Samtools)
5. Variant calling and filtering (GATK & VCFtools)

> **Note:** All code blocks are for **demonstration purposes only** and are not executed in this document.

---

# Dama Gazelle WGS Variant Calling Pipeline

This tutorial demonstrates the workflow for whole-genome sequencing (WGS) variant calling in the Dama gazelle. Each step includes the commands in a copyable code block.

---

## 1. Quality control using FastQC

```bash
#!/bin/bash -l
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Fastqc

module load fastqc

THREADS=24
INDIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq"
OUTDIR="${INDIR}/1_Fastqc"

cd $INDIR
fastqc *_1.fastq *_2.fastq -o $OUTDIR -t $THREADS
```

## 2. Adapter trimming with Trim Galore

```bash
module load trimgalore-0.6.7
# Ensure cutadapt v5.1 is available

THREADS=24
INDIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq"
OUTDIR="${INDIR}/2_Adapter_trimming"

cd $INDIR

for R1 in *_1.fastq; do
    SAMPLE=${R1%_1.fastq}
    R2="${SAMPLE}_2.fastq"
    echo "Processing sample: $SAMPLE"
    trim_galore --paired --cores $THREADS --output_dir $OUTDIR $R1 $R2
done
```

## 3a. Index the reference

```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/

module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17

REFERENCE=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/Dama_gazelle_hifiasm-ULONT_primary.fasta

samtools faidx $REFERENCE
bwa index $REFERENCE
```

## 3b. Align reads and process BAMs

```bash
# Change directory to where the T2T reference is stored.
cd /scratch/bistbs/   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load required modules                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# T2T consertium generated reference assembly#
# Indexing with samtools and BWA             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

REFERENCE=/scratch/bistbs/Dama_gazelle_hifiasm-ULONT_primary.fasta

# Index with samtools (for downstream use)
samtools faidx $REFERENCE

# Index with BWA (for mapping)
bwa index $REFERENCE
```
## 3C. Align reads and process BAMs using BWA-MEM

```bash
# Change directory to where the T2T reference is stored.
cd /scratch/bistbs/  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load required modules                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17 

# Directories
INPUT_DIR=/scratch/bistbs/
OUTPUT_DIR=/scratch/bistbs/4_Aligning_BWA/Alignment_output
REFERENCE=/scratch/bistbs/Dama_gazelle_hifiasm-ULONT_primary.fasta

# Make output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Number of threads (matches ntasks-per-node)
THREADS=24

# Loop through all paired-end files
for R1 in $INPUT_DIR/*_1_val_1.fq; do
    BASE=$(basename $R1 _1_val_1.fq)
    R2="${INPUT_DIR}/${BASE}_2_val_2.fq"
    
    echo "Processing sample: $BASE"
    echo "Read1: $R1"
    echo "Read2: $R2"
    
# Align reads with BWA-MEM and output BAM (mapped only) and Process and filter unmapped reads with samtools
#bF 4-binary format and remove reads that are 4. 4 means reads are unmapped.
    bwa mem -t $THREADS $REFERENCE $R1 $R2 \
        | samtools view -bF 4 - > $OUTPUT_DIR/${BASE}_mapped.bam
done
```    
## 3d. Sort the mapped SAM files
```bash
#!/bin/bash -l
#SBATCH --job-name=Sort_SAM
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/logs/SortSAM_%A_%a.out
#SBATCH --error=/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/logs/SortSAM_%A_%a.err
#SBATCH --array=1-5

# Load Java
module load java-20

# Paths
PICARD_JAR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/picard.jar"
INPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1"
OUTPUT_DIR="${INPUT_DIR}/Sorted_BAMs"
SCRATCH="/scratch/bistbs_new/tmp_sortbam"
LOG_DIR="${INPUT_DIR}/logs"

# Create directories (must exist before job submission)
mkdir -p "$OUTPUT_DIR" "$SCRATCH" "$LOG_DIR"

# List all BAM files
bam_files=($(ls -1 "${INPUT_DIR}"/*mapped.bam))
NUM_BAMS=${#bam_files[@]}

# Check array index
if [ $SLURM_ARRAY_TASK_ID -gt $NUM_BAMS ] || [ $SLURM_ARRAY_TASK_ID -lt 1 ]; then
    echo "❌ SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of range (1-$NUM_BAMS)"
    exit 1
fi

# Pick the BAM for this task
this_bam="${bam_files[$((SLURM_ARRAY_TASK_ID-1))]}"
base=$(basename "$this_bam" .bam)

# Unique TMP per task
TASK_TMP="${SCRATCH}/job_${SLURM_ARRAY_JOB_ID}_task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TASK_TMP"

echo "📂 Processing BAM: $this_bam"
echo "🏷 Output BAM: ${OUTPUT_DIR}/${base}_sorted.bam"
echo "🖥 Host: $HOSTNAME"
echo "🗂 TMP_DIR: $TASK_TMP"

# Run Picard SortSam
java -Xmx120g -jar "$PICARD_JAR" SortSam \
    I="$this_bam" \
    O="${OUTPUT_DIR}/${base}_sorted.bam" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="$TASK_TMP"

# Check success
if [ $? -eq 0 ]; then
    echo "✅ Finished sorting ${base}.bam"
else
    echo "❌ Error sorting ${base}.bam"
    exit 1
fi

```

## 4. Add or Replace Groups

```bash
#!/bin/bash -l
#SBATCH --job-name=AddReadGroups
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/AddRG_%A_%a.out
#SBATCH --error=logs/AddRG_%A_%a.err
#SBATCH --array=1-5   # <-- One job per BAM file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Load modules and variables        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load java-20

PICARD_JAR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/picard.jar"
INPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/Sorted_BAMs"
OUTPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/ReadGroups"
SCRATCH="/scratch/bistbs_new/tmp_addRG"

mkdir -p "$OUTPUT_DIR" "$SCRATCH" logs

# Get all BAM files
bam_list=(${INPUT_DIR}/*_mapped_sorted.bam)
this_bam=${bam_list[$((SLURM_ARRAY_TASK_ID-1))]}

# Safety check
if [ ! -f "$this_bam" ]; then
    echo "❌ No BAM file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

base=$(basename "$this_bam" .bam)

# Extract sample info from filename
# Example: SRR17129394_mapped_sorted.bam → RGID=17129394, RGSM=SRR17129394
RGSM=$(echo "$base" | cut -d'_' -f1)
RGID=$(echo "$base" | cut -d'_' -f1)
RGPU="${RGSM}_unit1"

echo "📂 Processing: $this_bam"
echo "🧬 RGID=$RGID, RGPU=$RGPU, RGSM=$RGSM"

# Run Picard AddOrReplaceReadGroups
java -Xmx12g -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$this_bam" \
    O="${OUTPUT_DIR}/${base}_RG.bam" \
    RGID="$RGID" \
    RGLB="lib1" \
    RGPL="illumina" \
    RGPU="$RGPU" \
    RGSM="$RGSM" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="$SCRATCH"

echo "✅ Finished adding read groups to ${base}.bam"
```
---
