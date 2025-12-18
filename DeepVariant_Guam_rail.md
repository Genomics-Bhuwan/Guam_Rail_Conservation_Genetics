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
# Keep these as absolute paths for the host system
INPUT_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG
OUTPUT_DIR=${INPUT_DIR}/DeepVariant
REF_GENOME=bHypOws1_hifiasm.bp.p_ctg.fasta
# Since the BAM is in a subdirectory of INPUT_DIR, we use the relative path from INPUT_DIR
BAM_FILE=rm_duplicates_BAM/Guam_Rail_merged_downsampled.bam

# ----------------Run DeepVariant---------------- #
$APPTAINER run \
  -B ${INPUT_DIR}:/input \
  -B ${OUTPUT_DIR}:/output \
  ${DEEPVARIANT_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/${REF_GENOME} \
    --reads=/input/${BAM_FILE} \
    --output_vcf=/output/Guam_Rail_merged_downsampled.dv.vcf \
    --output_gvcf=/output/Guam_Rail_merged_downsampled.dv.g.vcf \
    --num_shards=$(nproc) \
    --vcf_stats_report=true \
    --logging_dir=/output/logs
```
