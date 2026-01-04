#### How to do PSMC analyses?
- The PSMC is a pairwise sequential markovian coalescent model used for infering the effective population size or genetic diversity of the species based on the demographic history of the species.
- It is based on coalescent theory helping understand how genetic diversity is shaped by the history of a population in a pairwise fashion.

---
Assumptions
---
- I already have one BAM per individual in */scratch/bistbs/Population_Genomic_Analysis/PSMC/*, named exactly as the sample IDs (e.g. SRR17129394.bam, etc.).
- Reference FASTA is available and indexed at REF_FA=/scratch/bistbs/Population_Genomic_Analysis/PSMC/Dama_gazelle_hifiasm-ULONT_primary.fasta 
- bcftools, samtools, psmc and vcfutils.pl are available via module load.
- The script will produce gziped fq, psmcfa, psmc output and bootstrap results in the results/ directory.

##### Find our the generation time and mutation rate of your species.
Generation time= 5.85 years (https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.10876)
We use generation time GEN=5.85 years and mutation rate MU=2.96E-09 per site per year. (science.sciencemag.org/content/364/6446/eaav6202/suppl/DC1).

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
#SBATCH --array=0-2                    # One task per sample (3 samples)
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/PSMC_Guamrail/PSMC_%a.log
#SBATCH --error=/shared/jezkovt_bistbs_shared/Guam_Rail/PSMC_Guamrail/PSMC_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# Load modules
module load samtools-1.22.1
module load bcftools-1.15

# --------------------------
# Paths
# --------------------------
WORKDIR=/shared/jezkovt_bistbs_shared/Guam_Rail/PSMC_Guamrail
mkdir -p $WORKDIR
cd $WORKDIR

REF=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/bHypOws1_hifiasm.bp.p_ctg.fasta

# List of BAM files
BAMS=(
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/FMNH390989_downsampled.bam
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/HOW_N23-0063_downsampled.bam
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/HOW_N23-0568_downsampled.bam
)

# Select BAM based on SLURM_ARRAY_TASK_ID
BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename $BAM .bam)

echo "Processing sample $SAMPLE ..."

# --------------------------
# Step 1: Generate consensus FASTQ
# --------------------------
bcftools mpileup -Ou -f $REF $BAM | \
bcftools call -c | \
vcfutils.pl vcf2fq -d 10 -D 100 > ${SAMPLE}.fq

echo "Consensus FASTQ generated: ${WORKDIR}/${SAMPLE}.fq"

```

1. `bcftools mpileup -C50 -uf <reference_genome> <bam_file>`: This command generates a textual pileup format of the input BAM file (`<bam_file>`) using the given reference genome (`<reference_genome>`). The `C50` option applies a coefficient to adjust the base alignment quality, and the `u` flag outputs the results in the uncompressed BCF format, which is required for piping to `bcftools`. The `f` flag specifies the reference genome file.
2. `bcftools call -c`: This command performs variant calling on the input data received from the `bcftools mpileup` command (indicated by `` as input). The `c` option uses the consensus caller, which is suitable for calling a diploid consensus sequence.
3. `vcfutils.pl vcf2fq -d 10 -D 100`: This command is part of the `bcftools` package and converts the output from `bcftools call` (in VCF format) to a FastQ format. The `d 10` and `D 100` options set the minimum and maximum depth thresholds for filtering variants, respectively.
4. `gzip > <output.fq.gz>`: This part of the command compresses the final output using `gzip` and saves it as a `.fq.gz` file (`<output.fq.gz>`).


- With this consensus sequence, we created the input file to run PSMC.

  ##### The slurm script will convert the fastq file into psmcfa
```bash
#!/bin/bash -l
#SBATCH --job-name=PSMC_run
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=psmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --array=0-4  # 5 samples

# Binaries
PSMC_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/psmc/psmc
FQ2PSMCFA_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/psmc/utils/fq2psmcfa
SPLITFA_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/psmc/utils/splitfa
PSMC_PLOT_BIN=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/rm_duplicates_BAM/downsampled/psmc/utils/psmc_plot.pl

# Mutation rate and generation time
MU=2.96e-09
GEN=5.85

# FASTQ files
FQS=(SRR17129394.fq.gz SRR17134085.fq.gz SRR17134086.fq.gz SRR17134087.fq.gz SRR17134088.fq.gz)

for FQ in "${FQS[@]}"
do
    SAMPLE=$(basename $FQ .fq.gz)
    OUTDIR=PSMC_results/${SAMPLE}
    
    mkdir -p $OUTDIR

    echo "========================================"
    echo "Processing sample: $SAMPLE"
    echo "Saving results in: $OUTDIR"
    echo "========================================"

    # Step 1: Convert FASTQ → PSMCFA
    $FQ2PSMCFA_BIN -q20 $FQ > $OUTDIR/${SAMPLE}.psmcfa

    # Step 2: Split
    $SPLITFA_BIN $OUTDIR/${SAMPLE}.psmcfa > $OUTDIR/${SAMPLE}_split.psmcfa

done
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
$PSMC_BIN -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o $OUTDIR/${SAMPLE}.psmc \
    $OUTDIR/${SAMPLE}.psmcfa

echo "Main PSMC done."
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
GEN=5.85
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

##### 5.B. Multiple sample plot or the PSCM curves

```bash
# Step 1 — Create a combined plotting directory
mkdir -p /scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/combined_plot

# Step 2 — Copy all combined PSMC files into this folder
cd /scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/combined_plot

cp ../SRR17129394/SRR17129394.combined.psmc .
cp ../SRR17134085/SRR17134085.combined.psmc .
cp ../SRR17134086/SRR17134086.combined.psmc .
cp ../SRR17134087/SRR17134087.combined.psmc .
cp ../SRR17134088/SRR17134088.combined.psmc .

# Step 3 — Create names file (recommended)
This file assigns custom colors & sample labels.

Create the file names.txt:

SRR17129394   Addra_1
SRR17134085   Addra_2
SRR17134086   Addra_3
SRR17134087   Mhorr_1
SRR17134088   Mhorr_2

(You can customize labels however you want.)

⚡ Step 4 — Run ONE command to generate the combined multi-sample plot

Use psmc_plot.pl with multiple inputs:

PSMC_PLOT_BIN=/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/utils/psmc_plot.pl

GEN=5.85
MU=2.96e-09

cd /scratch/bistbs/Population_Genomic_Analysis/PSMC/PSMC_results/combined_plot

$PSMC_PLOT_BIN \
   -g $GEN -u $MU -X 1000000 \
   -M names.txt \
   combined_plot \
   SRR17129394.combined.psmc \
   SRR17134085.combined.psmc \
   SRR17134086.combined.psmc \
   SRR17134087.combined.psmc \
   SRR17134088.combined.psmc

🎉 Output

Inside combined_plot/ you will get:

combined_plot.eps
combined_plot.pdf
combined_plot.svg
combined_plot.png

This one figure contains all 5 PSMC curves, each in a different color, with labels from your names.txt.

Want species-colored curves?

I can also generate:

Addra = blue

Mhorr = red

Custom line widths

Custom legend names

Export in publication-ready format

Just tell me!
````

