#### Variant calling using DeepVariant
```bash
#!/bin/bash -l
#SBATCH --time=70:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=DeepVariant

APPTAINER=/usr/bin/apptainer
DEEPVARIANT_IMAGE=docker://google/deepvariant:1.9.0

# ----------------User Variables------------------ #
INPUT_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/
OUTPUT_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/DeepVariant
REF_GENOME=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/bHypOws1_hifiasm.bp.p_ctg.fasta
BAM_FILE=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam

# ----------------Create output directories------------------ #
mkdir -p ${OUTPUT_DIR}/logs

# ----------------Run DeepVariant---------------- #
echo "Starting DeepVariant at $(date)"

$APPTAINER run \
  -B ${INPUT_DIR}:/input \
  -B ${OUTPUT_DIR}:/output \
  ${DEEPVARIANT_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/${REF_GENOME} \
    --reads=/input/${BAM_FILE} \
    --output_vcf=/output/${BAM_FILE%.bam}.dv.vcf \
    --output_gvcf=/output/${BAM_FILE%.bam}.dv.g.vcf \
    --num_shards=$(nproc) \
    --vcf_stats_report=true \
    --disable_small_model=true \
    --logging_dir=/output/logs

echo "DeepVariant finished at $(date)"
```
