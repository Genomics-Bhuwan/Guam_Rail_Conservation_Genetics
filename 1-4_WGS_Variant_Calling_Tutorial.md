# Guam rail Whole-Genome Resequencing (WGS) Variant Calling Pipeline Tutorial

**Author:** Bhuwan Singh Bist

This tutorial demonstrates a complete workflow for WGS variant calling in the Guam rail. The pipeline includes:
1. Quality control (FastQC)
2. Adapter trimming (Trim Galore + Cutadapt)
3. Mapping reads to a haplotype-resolved Ruminant Telomere-to-Telomere(T2T)reference genome assembly of Dama gazelle(Nanger dama)
4. BAM processing (Picard & Samtools)
5. Variant calling and filtering (GATK & VCFtools)

---

# Guam rail e WGS Variant Calling Pipeline

This tutorial demonstrates the workflow for whole-genome sequencing (WGS) variant calling in the Dama gazelle. Each step includes the commands in a copyable code block.

---

## 1. Quality control using FastQC

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_guamrail
#SBATCH --output=fastqc_guamrail.out
#SBATCH --error=fastqc_guamrail.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=40G
#SBATCH --time=12:00:00

module load fastqc

THREADS=24

INDIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis"
OUTDIR="${INDIR}/fastqc_results"

mkdir -p $OUTDIR
cd $INDIR

# MATCH YOUR FILE FORMAT: *.fastq-*.gz
fastqc *.fastq-*.gz -o $OUTDIR -t $THREADS

```

## 2. Adapter trimming with Trim Galore

```bash
#!/bin/bash
#SBATCH --job-name=trimgalore_guamrail
#SBATCH --output=trimgalore_guamrail_%A_%a.out
#SBATCH --error=trimgalore_guamrail_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=200:00:00
#SBATCH --array=1-4   # Number of samples in the list

module load trimgalore-0.6.7

THREADS=8
INDIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis"
OUTDIR="${INDIR}/trimmed_fastq"
mkdir -p $OUTDIR

cd $INDIR

# List of all sample prefixes
SAMPLES=("FMNH390989" "HOW_N23-0063" "HOW_N23-0568" "KSW5478")

# Pick the sample for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Match R1 and R2 files with any trailing number
R1=$(ls ${SAMPLE}_1.fastq-*.gz)
R2=$(ls ${SAMPLE}_2.fastq-*.gz)

echo "Processing sample: $SAMPLE"
echo "R1 file: $R1"
echo "R2 file: $R2"

# Run Trim Galore in paired-end mode
trim_galore --paired --cores $THREADS --output_dir $OUTDIR $R1 $R2

```

## 3a. Index the reference

```bash
- Change to your indexing directory
- Load required modules
module load samtools-1.22.1
module load bwa-0.7.17

cd /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Indexing_reference_genome_assembly/

- Set path to your reference genome assembly
REFERENCE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Indexing_reference_genome_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta

- Index reference genome with samtools (creates .fai file)
samtools faidx $REFERENCE

- Index reference genome with bwa (creates .bwt, .pac, .ann, .amb, .sa files)
bwa index $REFERENCE

echo "Indexing of $REFERENCE completed!"

```


## 3b. Align reads and process BAMs using BWA-MEM

```bash
# Change directory to where the reference is stored.
#!/bin/bash
#SBATCH --job-name=align_BWA_array
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Alignment_%A_%a.log
#SBATCH --array=0-3                  # Array job for 4 samples (0,1,2,3)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --partition=batch

# load the modules.
module load samtools-1.22.1
module load bwa-0.7.17

# Directories for the files.
INPUT_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/trimmed_fastq
OUTPUT_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem
REFERENCE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Indexing_reference_genome_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta

mkdir -p $OUTPUT_DIR

# Define the sample list

SAMPLES=(
"FMNH390989_1.fastq-003.gz_val_1.fq.gz FMNH390989_2.fastq-002.gz_val_2.fq.gz FMNH390989"
"HOW_N23-0063_1.fastq-001.gz_val_1.fq.gz HOW_N23-0063_2.fastq-002.gz_val_2.fq.gz HOW_N23-0063"
"HOW_N23-0568_1.fastq-004.gz_val_1.fq.gz HOW_N23-0568_2.fastq-005.gz_val_2.fq.gz HOW_N23-0568"
"KSW5478_1.fastq-004.gz_val_1.fq.gz KSW5478_2.fastq-001.gz_val_2.fq.gz KSW5478"
)

# Select the sample for this array job.
SAMPLE_PAIR=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R1=$(echo $SAMPLE_PAIR | cut -d' ' -f1)
R2=$(echo $SAMPLE_PAIR | cut -d' ' -f2)
BASE=$(echo $SAMPLE_PAIR | cut -d' ' -f3)

echo "Starting alignment for sample: $BASE"
echo "Read1: $INPUT_DIR/$R1"
echo "Read2: $INPUT_DIR/$R2"

# Align the reads with BWA-MEM
# Sort and Index the BAM.

bwa mem -t $SLURM_CPUS_PER_TASK $REFERENCE $INPUT_DIR/$R1 $INPUT_DIR/$R2 \
    | samtools view -bF 4 - \
    | samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTPUT_DIR/${BASE}_mapped_sorted.bam

## Index the sorted BAM
samtools index $OUTPUT_DIR/${BASE}_mapped_sorted.bam

echo "Finished alignment and indexing for $BASE"

```
``

##### 4. Add or Replace Groups

```bash
#!/bin/bash -l
#SBATCH --job-name=AddReadGroups
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/AddRG_%A_%a.out
#SBATCH --error=logs/AddRG_%A_%a.err
#SBATCH --array=1-4

# Load modules
module load java-20
module load samtools-1.22.1

# Path to the directory
BASE_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem"

PICARD_JAR="${BASE_DIR}/Add_RG/picard.jar"
INPUT_DIR="${BASE_DIR}"
OUTPUT_DIR="${BASE_DIR}/Add_RG"
SCRATCH="/shared/jezkovt_bistbs_shared/tmp_addRG"

mkdir -p "$OUTPUT_DIR" "$SCRATCH" logs

# Get the BAM this task.
bam_list=(${INPUT_DIR}/*_mapped_sorted.bam)
this_bam=${bam_list[$((SLURM_ARRAY_TASK_ID-1))]}

if [ ! -f "$this_bam" ]; then
    echo " No BAM file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

base=$(basename "$this_bam" .bam)

# Read the group fields.

RGSM=$(echo "$base" | cut -d'_' -f1)
RGID="$RGSM"
RGPU="${RGSM}_unit1"

echo " Processing: $this_bam"
echo " RGID=$RGID | RGSM=$RGSM | RGPU=$RGPU"
echo " CPUs allocated: $SLURM_CPUS_PER_TASK"

# Add teh Read Groups
java -Xmx90g \
  -XX:ParallelGCThreads=14 \
  -Dsamjdk.use_async_io_read_samtools=true \
  -Dsamjdk.use_async_io_write_samtools=true \
  -Dsamjdk.use_async_io_write_tribble=true \
  -Dsamjdk.compression_level=2 \
  -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$this_bam" \
    O="${OUTPUT_DIR}/${base}_RG.bam" \
    RGID="$RGID" \
    RGLB="lib1" \
    RGPL="ILLUMINA" \
    RGPU="$RGPU" \
    RGSM="$RGSM" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="$SCRATCH"

# Index the output BAM
samtools index "${OUTPUT_DIR}/${base}_RG.bam"

echo " Finished adding read groups to ${base}.bam"

```
#### 5. Combine all the BAMs together
```bash
java -Xmx90G -jar /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/picard.jar MergeSamFiles \
    I=FMNH390989_RG.bam \
    I=HOW_N23-0063_RG.bam \
    I=HOW_N23-0568_RG.bam \
    I=KSW5478_RG.bam \
    O=Guam_Rail_merged.bam \
    USE_THREADING=true \
    CREATE_INDEX=true
```

#### 6. Remove the Duplicates from the combined BAMs.

```bash
#!/bin/bash -l
#SBATCH --job-name=MarkDuplicates
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/MarkDup_%A.out
#SBATCH --error=logs/MarkDup_%A.err

# Load Java
module load java-20

# Picard JAR
PICARD_JAR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/picard.jar"

# Input BAM
INPUT_BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/Guam_Rail_merged.bam"

# Output directory
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM"

# Temporary directory (INSIDE shared, will be removed)
TMP_DIR="${OUTPUT_DIR}/tmp_MarkDuplicates"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" logs

# Output files
MARKED_BAM="${OUTPUT_DIR}/Guam_Rail_merged_rmdup.bam"
METRICS_FILE="${OUTPUT_DIR}/Guam_Rail_merged_rmdup.metrics"

echo "Running Picard MarkDuplicates"
echo "Input BAM : $INPUT_BAM"
echo "Output BAM: $MARKED_BAM"

java -Xmx100g -jar "$PICARD_JAR" MarkDuplicates \
    I="$INPUT_BAM" \
    O="$MARKED_BAM" \
    METRICS_FILE="$METRICS_FILE" \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="$TMP_DIR" \
    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
    MAX_RECORDS_IN_RAM=50000 \
    CREATE_INDEX=true

echo "MarkDuplicates finished"

# REMOVE temporary directory
rm -rf "$TMP_DIR"
echo "Temporary directory removed"
```

---
