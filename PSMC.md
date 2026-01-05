#### How to do PSMC analyses?
- The PSMC is a pairwise sequential markovian coalescent model used for infering the effective population size or genetic diversity of the species based on the demographic history of the species.
- It is based on coalescent theory helping understand how genetic diversity is shaped by the history of a population in a pairwise fashion.
- mu <- 1.45e-8
- g  <- 3.4

```bash
mkdir psmc
```

##### Consensus building
- The first step is to generate the consensus sequence from the bam files.
- Although the command for this step is relatively old, it still functions effectively.
- However, it is important to ensure that the command is still operational and functioning correctly before proceeding.
```bash
 #!/bin/bash -l
#SBATCH --job-name=PSMC_Guamrail
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --array=0-2
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_%a.log
#SBATCH --error=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# --------------------------
# Load modules
# --------------------------
module load samtools-1.22.1
module load bcftools-1.15

# --------------------------
# Reference genome
# --------------------------
REF=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/bHypOws1_hifiasm.bp.p_ctg.fasta

# --------------------------
# Output directory (PSMC FASTQs)
# --------------------------
OUTDIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail
mkdir -p $OUTDIR
cd $OUTDIR

# --------------------------
# BAM files
# --------------------------
BAMS=(
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/FMNH390989_downsampled.bam
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/HOW_N23-0063_downsampled.bam
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/HOW_N23-0568_downsampled.bam
)

BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename $BAM .bam)

echo "Processing sample: $SAMPLE"
echo "Log file: PSMC_${SLURM_ARRAY_TASK_ID}.log"
echo "Error file: PSMC_${SLURM_ARRAY_TASK_ID}.err"

# --------------------------
# Generate consensus FASTQ
# --------------------------
bcftools mpileup -Ou -f $REF $BAM | \
bcftools call -c | \
vcfutils.pl vcf2fq -d 10 -D 100 > ${SAMPLE}.fq

echo "Consensus FASTQ generated:"
echo "${OUTDIR}/${SAMPLE}.fq"
```



1. `bcftools mpileup -C50 -uf <reference_genome> <bam_file>`: This command generates a textual pileup format of the input BAM file (`<bam_file>`) using the given reference genome (`<reference_genome>`). The `C50` option applies a coefficient to adjust the base alignment quality, and the `u` flag outputs the results in the uncompressed BCF format, which is required for piping to `bcftools`. The `f` flag specifies the reference genome file.
2. `bcftools call -c`: This command performs variant calling on the input data received from the `bcftools mpileup` command (indicated by `` as input). The `c` option uses the consensus caller, which is suitable for calling a diploid consensus sequence.
3. `vcfutils.pl vcf2fq -d 10 -D 100`: This command is part of the `bcftools` package and converts the output from `bcftools call` (in VCF format) to a FastQ format. The `d 10` and `D 100` options set the minimum and maximum depth thresholds for filtering variants, respectively.
4. `gzip > <output.fq.gz>`: This part of the command compresses the final output using `gzip` and saves it as a `.fq.gz` file (`<output.fq.gz>`).


- With this consensus sequence, we created the input file to run PSMC.

  ##### The slurm script will convert the fastq file into psmcfa
```bash

#!/bin/bash -l
#SBATCH --job-name=PSMC_preprocess
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --array=0-2
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/PSMC_pre_%A_%a.log
#SBATCH --error=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/PSMC_pre_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# --------------------------
# PSMC binaries (YOUR paths)
# --------------------------
PSMC_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/psmc/psmc
FQ2PSMCFA_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/psmc/utils/fq2psmcfa
SPLITFA_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/psmc/utils/splitfa
PSMC_PLOT_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/psmc/utils/psmc_plot.pl

# --------------------------
# Input FASTQ directory
# --------------------------
FQDIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail

# --------------------------
# FASTQ files (3 Guam rail samples)
# --------------------------
FQS=(
FMNH390989_downsampled.fq
HOW_N23-0063_downsampled.fq
HOW_N23-0568_downsampled.fq
)

# --------------------------
# Select FASTQ based on array ID
# --------------------------
FQ=${FQS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename $FQ .fq)

OUTDIR=${FQDIR}/PSMC_results/${SAMPLE}
mkdir -p $OUTDIR
cd $OUTDIR

echo "========================================"
echo "Processing sample: $SAMPLE"
echo "Input FASTQ: ${FQDIR}/${FQ}"
echo "Output directory: $OUTDIR"
echo "========================================"

# --------------------------
# Step 1: FASTQ → PSMCFA
# --------------------------
$FQ2PSMCFA_BIN -q20 ${FQDIR}/${FQ} > ${SAMPLE}.psmcfa

# --------------------------
# Step 2: Split PSMCFA
# --------------------------
$SPLITFA_BIN ${SAMPLE}.psmcfa > ${SAMPLE}_split.psmcfa

echo "PSMC preprocessing complete for $SAMPLE"

```

##### Step 3: Main PSMC
- Psmc: command to run PSMC tool.
- N25: sets effective population size to 25. Ne is used for calculating the time to the most recent common ancestor of the population.
- t15: flag sets the scaled mutation rate per generation(t) to 15. It is the product of mutation rate per base pair per generation and the effective population size.
- r5: sets the scaled This flag sets the scaled recombination rate per generation (r) to 5.
- The scaled recombination rate is the product of the recombination rate per base pair per generation and the effective population size.
- p "4+25*2+4+6": This flag sets the time intervals (p) for the PSMC model. The specified pattern, "4+25*2+4+6", means that there are 4 intervals of equal size at the start, followed by 25 intervals with twice the size of the previous intervals, and then 4 more intervals of equal size, and finally 6 more intervals of increasing size. This allows the model to have higher time resolution near the present and lower resolution in the more distant past.
o <output.psmc>: This flag specifies the output file name for the PSMC results. Replace <output.psmc> with the desired output file name.
<input.psmcfa>: This is the input file in PSMCFA format, which contains the sequence data to be analyzed. Replace <input.psmcfa> with the name of the input file.

```bash
#!/bin/bash -l
#SBATCH --job-name=PSMC_run
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --array=0-2
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/PSMC_results/Results/PSMC_run_%A_%a.log
#SBATCH --error=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/PSMC_results/Results/PSMC_run_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# --------------------------
# PSMC binary
# --------------------------
PSMC_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail/psmc/psmc

# --------------------------
# Input/output directories
# --------------------------
FQDIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/PSMC_Guamrail
OUTDIR=${FQDIR}/PSMC_results/Results
mkdir -p $OUTDIR

# --------------------------
# Samples
# --------------------------
SAMPLES=(
FMNH390989_downsampled
HOW_N23-0063_downsampled
HOW_N23-0568_downsampled
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "========================================"
echo "Running PSMC for sample: $SAMPLE"
echo "Input PSMCFA: ${FQDIR}/PSMC_results/${SAMPLE}.psmcfa"
echo "Output PSMC: ${OUTDIR}/${SAMPLE}.psmc"
echo "========================================"

# --------------------------
# Run main PSMC
# --------------------------
$PSMC_BIN -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o ${OUTDIR}/${SAMPLE}.psmc \
    ${FQDIR}/PSMC_results/${SAMPLE}.psmcfa

echo "Main PSMC done for $SAMPLE."

```
##### Step 4: Bootstrap the samples(100)
- I ran individual batch script for each of the 5 samples using command below. 
```bash
#!/bin/bash -l
#SBATCH --job-name=88_run                  # Job name
#SBATCH --time=100:00:00                   # Walltime: 100 hours
#SBATCH --cpus-per-task=16                 # Number of CPU cores per task
#SBATCH --mem=128G                         # Memory per node
#SBATCH --partition=batch                  # Partition/queue
#SBATCH --output=psmc_run_88_%A_%a.log    # Output log file
#SBATCH --error=psmc_run_88_%A_%a.err     # Error log file
#SBATCH --mail-type=END,FAIL               # Email notifications for job end/failure
#SBATCH --mail-user=bistbs@miamioh.edu    # Your email

# Load necessary modules (if any)
# module load psmc

# Path to PSMC binary
PSMC_BIN=/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc

# Sample directory
SAMPLE_DIR=/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17134088
SAMPLE=$(basename "$SAMPLE_DIR")

# Number of bootstrap replicates
BOOTSTRAPS=100

# Number of parallel processes (should not exceed CPUs)
PARALLEL=16

# Change to the sample directory
cd "$SAMPLE_DIR" || { echo "Cannot cd into $SAMPLE_DIR"; exit 1; }

echo "Starting bootstrapping for $SAMPLE using $PARALLEL cores ..."

# Run bootstrap replicates in parallel
seq $BOOTSTRAPS | xargs -P $PARALLEL -I{} \
    $PSMC_BIN -N25 -t15 -r5 -b -p "4+25*2+4+6" \
    -o ${SAMPLE}_round-{}.psmc ${SAMPLE}_split.psmcfa

# Combine main PSMC with bootstrap replicates
cat ${SAMPLE}.psmc ${SAMPLE}_round-*.psmc > ${SAMPLE}.combined.psmc

echo "Bootstrapping completed for $SAMPLE"

```

#####  Step 5. Final Plotting
###### Step 5.A. Generate separate plots for each sample (one plot file per sample)
```bash
# Path to psmc_plot.pl
PSMC_PLOT_BIN=/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/utils/psmc_plot.pl

# Generation time & mutation rate
GEN=3.4
MU=2.96e-09

# List of your sample directories
SAMPLES=(
"/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17129394"
"/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17134085"
"/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17134086"
"/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17134087"
"/scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/SRR17134088"
)

# Loop over all samples
for SAMPLE_DIR in "${SAMPLES[@]}"
do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "Plotting PSMC for $SAMPLE ..."

    cd "$SAMPLE_DIR"

    # Required input file: SAMPLE.combined.psmc
    $PSMC_PLOT_BIN -g $GEN -u $MU -X 1000000 \
        $SAMPLE $SAMPLE.combined.psmc

    echo "Finished plotting for $SAMPLE"
done

echo "All PSMC plots completed."

```

#### Plotting using R
```bash
################################################################################
# Guam Rail PSMC plotting
# MAIN plot + per-individual bootstrap plots
################################################################################

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------
# Working directory
# ------------------------------------------------------------------
setwd("F:/Collaborative_Projects/Guam_Rail/PSMC")

# ------------------------------------------------------------------
# SOURCE PSMC FUNCTIONS
# ------------------------------------------------------------------
source("plotPsmc.r")

# ------------------------------------------------------------------
# Parameters (Guam Rail)
# ------------------------------------------------------------------
mu <- 1.45e-8
g  <- 3.4

# ------------------------------------------------------------------
# Sample IDs (must match filenames exactly)
# ------------------------------------------------------------------
ordered_ids <- c(
  "FMNH390989_downsampled",
  "HOW_N23-0063_downsampled",
  "HOW_N23-0568_downsampled"
)

# Clean labels for the legend (removes "_downsampled")
legend_labels <- gsub("_downsampled", "", ordered_ids)

# Map colors to the CLEANED labels
ind_colors_named <- setNames(
  c("#0072B2", "#D55E00", "#009E73"),
  legend_labels
)

# ------------------------------------------------------------------
# Identify main PSMC files
# ------------------------------------------------------------------
all_files <- list.files(pattern = "\\.psmc$")
psmc_main_files <- all_files[!grepl("\\.combined\\.psmc$", all_files)]

psmc_main_files <- sapply(ordered_ids, function(id) {
  f <- grep(paste0("^", id, "\\.psmc$"), psmc_main_files, value = TRUE)
  if (length(f) == 0) stop("Missing main PSMC file for ", id)
  f[1]
}, USE.NAMES = FALSE)

# ------------------------------------------------------------------
# Read MAIN PSMC files
# ------------------------------------------------------------------
main_list <- vector("list", length(psmc_main_files))

for (i in seq_along(psmc_main_files)) {
  f <- psmc_main_files[i]
  message("Reading main PSMC: ", f)
  
  res <- psmc.result(file = f, mu = mu, g = g, i.iteration = 25)
  df  <- bind_rows(res, .id = "iter")
  
  df$SampleID <- ordered_ids[i]
  df$Label    <- legend_labels[i] # Uses the cleaned label here
  
  main_list[[i]] <- df
}

main_df <- bind_rows(main_list)

main_df_main <- main_df %>%
  filter(iter == "1") %>%
  group_by(SampleID) %>%
  slice(9:n()) %>%
  ungroup()

# ------------------------------------------------------------------
# MAIN PLOT (no bootstraps)
# ------------------------------------------------------------------
p_main <- ggplot(main_df_main, aes(x = YearsAgo, y = Ne, color = Label)) +
  geom_step(linewidth = 1.6, direction = "hv") +
  scale_color_manual(values = ind_colors_named, name = "Sample") +
  scale_x_log10(
    limits = c(1e4, 5e6),                  # X-axis 10k–5Mya
    breaks = c(1e4, 5e4, 1e5, 5e5, 1e6, 5e6),
    labels = c("10 Kya","50 Kya","100 Kya","500 Kya","1 Mya","5 Mya")
  ) +
  scale_y_log10(
    limits = c(7e3, 3.5e5),                # Y-axis 7k–350k
    breaks = c(7e3,1e4,2e4,5e4,1e5,2e5,3e5,3.5e5),
    labels = c("7k","10k","20k","50k","100k","200k","300k","350k")
  ) +
  annotation_logticks(sides = "bl") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12)
  ) +
  labs(
    x = "Time",
    y = "Effective population size (Ne)"
  ) +
  annotate(
    "text",
    x = 5e6,
    y = 8e3,                                # slightly above 7k for label
    label = "(μ = 1.45e-8, g = 3.4)",
    hjust = 1,
    vjust = 0,
    size = 5,
    fontface = "bold"
  )

print(p_main)

ggsave("GuamRail_PSMC_Main_NoBootstraps_5Mya_7k-350k.pdf",
       p_main, width = 14, height = 8)
ggsave("GuamRail_PSMC_Main_NoBootstraps_5Mya_7k-350k.jpeg",
       p_main, width = 14, height = 8, dpi = 300)

# ------------------------------------------------------------------
# BOOTSTRAP PLOTS (if combined files exist)
# ------------------------------------------------------------------
combined_dir <- "All_combined_PSMC"

for (i in seq_along(ordered_ids)) {
  
  sample_id <- ordered_ids[i]
  clean_id  <- legend_labels[i] # Cleaned version for title
  
  combined_file <- file.path(combined_dir,
                             paste0(sample_id, ".combined.psmc"))
  
  if (!file.exists(combined_file)) {
    message("No bootstrap file for ", sample_id, " — skipping.")
    next
  }
  
  message("Reading bootstrap file: ", combined_file)
  
  res_all <- psmc.result(file = combined_file,
                         mu = mu, g = g, i.iteration = 25)
  
  df_all <- bind_rows(res_all, .id = "iter") %>%
    group_by(iter) %>%
    slice(9:n()) %>%
    ungroup()
  
  df_main_ind <- df_all %>% filter(iter == "1")
  df_boot_ind <- df_all %>% filter(iter != "1")
  
  p_boot <- ggplot() +
    geom_step(
      data = df_boot_ind,
      aes(x = YearsAgo, y = Ne, group = iter),
      color = "grey60", alpha = 0.35,
      linewidth = 0.6, direction = "hv"
    ) +
    geom_step(
      data = df_main_ind,
      aes(x = YearsAgo, y = Ne),
      color = ind_colors_named[clean_id], # Reference the color by cleaned name
      linewidth = 1.6, direction = "hv"
    ) +
    scale_x_log10(
      limits = c(1e4, 5e6),                  # X-axis 10k–5Mya
      breaks = c(1e4, 5e4, 1e5, 5e5, 1e6, 5e6),
      labels = c("10 Kya","50 Kya","100 Kya","500 Kya","1 Mya","5 Mya")
    ) +
    scale_y_log10(
      limits = c(7e3, 3.5e5),                # Y-axis 7k–350k
      breaks = c(7e3,1e4,2e4,5e4,1e5,2e5,3e5,3.5e5),
      labels = c("7k","10k","20k","50k","100k","200k","300k","350k")
    ) +
    annotation_logticks(sides = "bl") +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA, linewidth = 1),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text  = element_text(size = 12)
    ) +
    labs(
      title = paste("PSMC bootstraps:", clean_id),
      x = "Years ago",
      y = "Effective population size (Ne)"
    )
  
  print(p_boot)
  
  ggsave(paste0(sample_id, "_PSMC_with_bootstraps_5Mya_7k-350k.pdf"),
         p_boot, width = 14, height = 8)
  ggsave(paste0(sample_id, "_PSMC_with_bootstraps_5Mya_7k-350k.jpeg"),
         p_boot, width = 14, height = 8, dpi = 300)
}

message("DONE: Guam Rail PSMC plots (5 Mya, 7k–350k) generated.")
```

