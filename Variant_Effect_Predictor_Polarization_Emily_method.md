################################# RAILS POLARIZING ############################

#### Virginia rail https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR18506199 and Eurasian coot https://trace.ncbi.nlm.nih.gov/Traces/?run=ERR12765175 as an outgroup.

#### Finding deleterious mutation in the Guam rail.
- Usually kept three species but here I am keeping two outgroups.
- Virginia rail:https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR18506199  and Eurasian coot (https://trace.ncbi.nlm.nih.gov/Traces/?run=ERR12765175)
#### Step 1. Download the short-read sequences of Virginia rail  and Eurasian coot.
- Adapter trimming, Map these with BWA-MEM using the reference genome assembly of Guam rail.
- Sort, rmduplicates and calculate the depth.
- Once the depth is calculated downsample them to the depth of the sample having the lowest depth.
- Virginia rail had the depth of 18X and I downsampled the depth of Eurasian rail from 42X to 18X.
- Merged them.

```bash
#!/bin/bash

#SBATCH --job-name=merge_outgroups
#SBATCH --output=merge_%j.log
#SBATCH --error=merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=95G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

# Load the correct Java version available on your cluster
module purge
module load java-20

# Define Variables
PICARD="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/picard.jar"
VIRGINIA="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/Virginia_rail/Alignment_Results/SRR18506199_final.bam"
COOT="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/Eurasian_coot_outgroup/trimmed_fastq/Alignment_Results/Processed_BAM/ERR12765175_coot_downsampled.bam"
OUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/Outgroups_Merged_Eurasian_coot_Virginia_18X"

# Create output directory
mkdir -p $OUT_DIR

# Run Picard MergeSamFiles
# Using java-20 will resolve the UnsupportedClassVersionError
java -Xmx90G -jar $PICARD MergeSamFiles \
    I=$VIRGINIA \
    I=$COOT \
    O=$OUT_DIR/Outgroups_Merged_Eurasian_coot_Virginia_18X.bam \
    USE_THREADING=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

echo "Merging complete at $(date)"
```

#### Step 2. Consensus generation for Virginia rail and Eurasian coot outgroups.

```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_Rail
#SBATCH --time=54:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/angsd_%j.out
#SBATCH --error=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/angsd_%j.err

# 1. Load the ANGSD module 
module load angsd

# 2. Define Input and Output
# Your merged outgroup BAM file
INPUT_BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/Outgroups_Merged_Eurasian_coot_Virginia_18X/Outgroups_Merged_Eurasian_coot_Virginia_18X.bam"
# Where the ancestral fasta will be saved
OUT_PREFIX="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus"

echo "Starting ANGSD Consensus Generation for Dama Gazelle Outgroups..."

# 3. Run ANGSD
# -i: Input BAM
# -P: Number of threads (CPUs to use)
# -doFasta 2: Generates a consensus sequence (The 'Ancestral' map)
# -doCounts 1: Counts the bases at each site (Needed for doFasta)
angsd \
    -i "$INPUT_BAM" \
    -P 12 \
    -doFasta 2 \
    -doCounts 1 \
    -out "$OUT_PREFIX"

echo "Process Complete. Your ancestral file is: ${OUT_PREFIX}.fa.gz"
```
#### Step 3. Conversion of vcf file to .bed format
- We need to take every mutation SNP that exists in your Guam rail biallelic SNPs and put their address into a simple list.
```bash
# Define your paths

IN_VCF="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Guam_rail_biallelic_SNPs.vcf.gz"
OUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Plink_For_Next_Step"
PLINK_EXE="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/plink"
# Run the conversion
$PLINK_EXE --vcf "$IN_VCF" \
--make-bed \
--allow-extra-chr \
--out "$OUT_DIR/Guam_Final_Binary"
```
- After running the plink, run below command
---
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $4}' Guam_Final_Binary.bim > Guam_SNPs_Coordinates.bed
---

#### Step 4. Extract ancestral alleles using bedtools
- Since, I already have .bed file and ANGSD consensus fasta.
```bash
module load bedtools-2.28
# 2. Define your paths
FASTA="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/rail_outgroup_consensus.fa"
BED="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Plink_For_Next_Step/Guam_SNPs_Coordinates.bed"
OUT_FILE="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Plink_For_Next_Step/ancestral_alleles.out"

# 3. Extract the ancestral base
# Note: Ensure the FASTA is unzipped first
bedtools getfasta -fi $FASTA -bed $BED -fo $OUT_FILE
```
#### Step 5. Reformat the Ancestral Alleles.
```bash
# 1. Extract just the DNA letters from the bedtools output
grep -v ">" ancestral_alleles.out > alleles_only.txt

# 2. Get the SNP IDs (Chromosome:Position) from your FIXED bed file
# Note: We use the 3rd column (actual position) to match PLINK's @:# format
awk '{print $1":"$3}' Guam_SNPs_Coordinates.bed > positions_only.txt

# 3. Paste them together to make the reference file
paste positions_only.txt alleles_only.txt > ancestral_alleles.txt
```

#### Step 6. Prepare the VCF IDs.
- PLINK needs the IDs in your vcf to match the IDs in your ancestral_alleles.txt
- A file named temp_Guam.vcf where every SNP is named Chromosome: Position.
```bash
# Navigate to your Plink folder
cd /shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Plink_For_Next_Step/

# Run the conversion
./plink --vcf /shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Guam_rail_biallelic_SNPs.vcf \
  --set-missing-var-ids @:# \
  --allow-extra-chr \
  --recode vcf \
  --const-fid 0 \
  --real-ref-alleles \
  --out temp_Guam
```

#### Step 7. Polarize (Big Switch)
- Now, tell PLINK to force the ancestral letter from your list to be the Reference(0) allele.
```bash
 1. Create a list of positions to keep
awk '{print $1}' ancestral_alleles.txt > ancestral_positions.txt

# 2. Run the polarization
./plink --vcf temp_Guam.vcf \
  --extract ancestral_positions.txt \
  --a2-allele ancestral_alleles.txt 2 1 \
  --recode vcf \
  --allow-extra-chr \
  --const-fid 0 \
  --out temp_Guam_polarized
```    

#### Step 8. Final Cleanup(Removing Mismatches)
- Remove the SNPs where the ancestral letter doesn't match the Guam_rail letters (e.g., Ancestral has "A", but Dama only has "C/G")
```bash
# 1. Identify the mismatches from the log file
grep 'Warning' temp_Guam_polarized.log | awk '{print $7}' | sed 's/.//;s/.$//' > mismatches.txt

# 2. Create the final, clean, polarized VCF
./plink --vcf temp_Guam_polarized.vcf \
  --exclude mismatches.txt \
  --allow-extra-chr \
  --const-fid 0 \
  --recode vcf \
  --out Guam_rail_POLARIZED_Final

# 3. Clean up the temporary files to save space
rm temp_Guam.vcf temp_Guam_polarized.vcf alleles_only.txt positions_only.txt
```

#### Step 9. Variant Annotation using VEP or SnpEff.
- The final vcf has the ALT allele (the '1') is officially the Derived Mutation.
- This is exactly what VEP needs to tell you if the mutation is harmful.
```bash
singularity exec -B /shared Final_VEP/ensembl-vep/ensembl-vep_latest.sif vep \
--input_file Guam_rail_POLARIZED_Final.vcf \
--output_file Final_VEP/Final_Annotation/Guam_rail_Annotated_Individuals.txt \
--fasta /shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/bHypOws1_hifiasm.bp.p_ctg.fasta \
--gff Final_VEP/Guam_rail.withExons.gff.gz \
--dir_cache Final_VEP/vep_cache \
--species Guam_rail \
--cache \
--offline \
--protein \
--biotype \
--pick \
--fork 8 \
--no_stats \
--buffer_size 10000 \
--force_overwrite \
--warning_file Final_VEP/vep_warnings.txt \
--individual all
```

#### Step 10. Filter sites by the Consequences
#### Step 10.a Filter site by the consequences
```bash
# Filter for LoF, missense, synonymous, and intergenic
# 1. Filter for Missense variants
filter_vep -i Final_VEP/Final_Annotation/Guam_rail_Annotated_Individuals.txt \
-o Final_VEP/Final_Annotation/missense_sites.txt --filter "Consequence is missense_variant"

# 2. Filter for Synonymous variants
filter_vep -i Final_VEP/Final_Annotation/Guam_rail_Annotated_Individuals.txt \
-o Final_VEP/Final_Annotation/synonymous_sites.txt --filter "Consequence is synonymous_variant"

# 3. Filter for Loss of Function (LoF) and Splicing variants
filter_vep -i Final_VEP/Final_Annotation/Guam_rail_Annotated_Individuals.txt \
-o Final_VEP/Final_Annotation/lof_sites.txt \
--filter "Consequence is transcript_ablation or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is inframe_insertion or Consequence is inframe_deletion or Consequence is splice_region_variant"

# 4. Filter for Intergenic variants
filter_vep -i Final_VEP/Final_Annotation/Guam_rail_Annotated_Individuals.txt \
-o Final_VEP/Final_Annotation/intergenic_sites.txt --filter "Consequence is intergenic_variant"
```
#### Step 10. b. Create ID Lists
- use awk to turn the VEP locations into a format vcftools understands(Chromosome and Position).
```bash
cat missense_sites.txt | awk '{ print $1 }' | awk '{sub(/\:/," ",$1)};1' > missense_IDs.txt
cat synonymous_sites.txt | awk '{ print $1 }' | awk '{sub(/\:/," ",$1)};1' > synonymous_IDs.txt
cat lof_sites.txt | awk '{ print $1 }' | awk '{sub(/\:/," ",$1)};1' > lof_IDs.txt
cat intergenic_sites.txt | awk '{ print $1 }' | awk '{sub(/\:/," ",$1)};1' > intergenic_IDs.txt
```
#### Step 10.c Extract the Genotypes from the polarized vcf
- Go to the polarized vcf and pull out the only SNPs taht fall into these categories.
```bash
for i in missense synonymous lof intergenic
do
    vcftools --vcf /shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/VEP_Polarization/rail_outgroup_consensus/Plink_For_Next_Step/Guam_rail_POLARIZED_Final.vcf \
    --recode \
    --recode-INFO-all \
    --positions ${i}_IDs.txt \
    --out Guam_rail_${i}_snps
done
```
#### Step 10.d. Convert to PLINK(The A-transpose format)
- This is the most important part for the R analysis.
- Use Plink --export A-transposae.
- This creates a .traw file where: Rows= SNPs, Columns - five individuals(SRR IDs); Values = 0, 1 or 2 (number of derived alleles).
```bash
for i in missense synonymous lof intergenic
do
    /home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/plink \
        --vcf Guam_rail_${i}_snps.recode.vcf \
        --export A-transpose \
        --allow-extra-chr \
        --out $OUTDIR/Guam_rail_${i}_genotypes
done
```

#### Step 11. Visualization for the plots
```bash
#### This is the final script for plotting the LOF and Missense mutation.

# 1. Load Libraries
library(tidyverse)
library(patchwork)

# 2. Set Working Directory
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/VEP/Polarization_VEP")

# 3. Define the files
files <- list(
  "Missense" = "Dama_Gazelle_missense_genotypes.traw",
  "LoF"      = "Dama_Gazelle_lof_genotypes.traw"
)

# 4. Process Data: Clean IDs and Assign Subspecies
get_genotype_metrics <- function(file_path, label) {
  dat <- read.table(file_path, header = TRUE, check.names = FALSE)
  geno <- dat[, 7:ncol(dat)]
  
  df <- data.frame(
    # Truncate IDs at the underscore
    SampleID  = names(geno) %>% str_remove("_.*"), 
    Total_Alt = colSums(geno, na.rm = TRUE),
    Hom_Alt   = colSums(geno == 2, na.rm = TRUE),
    Het       = colSums(geno == 1, na.rm = TRUE),
    Category  = label
  )
  
  # Assign Subspecies groups
  df <- df %>%
    mutate(Subspecies = case_when(
      str_detect(SampleID, "94$|85$|86$") ~ "Addra",
      str_detect(SampleID, "87$|88$")     ~ "Mhorr",
      TRUE                                ~ "Unknown"
    ))
  
  return(df)
}

all_metrics <- map2_df(files, names(files), get_genotype_metrics)

# 5. Plotting Function: Large Font for Molecular Ecology Resources
create_plot <- function(df, y_var, y_label, letter, cat_name) {
  
  plot_color <- ifelse(cat_name == "Missense", "#4C9A2A", "#F28E2B")
  
  p <- ggplot(df %>% filter(Category == cat_name), 
              aes(x = reorder(SampleID, !!sym(y_var)), y = !!sym(y_var))) +
    geom_col(fill = plot_color, color = "black", width = 0.7, alpha = 0.8) +
    geom_point(size = 3.5, color = "black") + 
    facet_grid(. ~ Subspecies, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    labs(title = paste0(letter, "   ", cat_name), y = y_label, x = "") +
    theme(
      # Axis Titles (Heterozygotes, etc.) - Large and Bold
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15), color = "black"),
      axis.title.x = element_blank(),
      
      # Axis Values (Sample names and Numbers) - High Legibility
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black", face = "bold"),
      
      # Subspecies Headers (Addra/Mhorr)
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
      
      # Panel Lettering (A, B, C...)
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      
      # Plot Layout Adjustments
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(linewidth = 1, color = "black"),
      panel.spacing = unit(1.5, "lines"),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p)
}

# 6. Generate the 6 Specific Panels
pA <- create_plot(all_metrics, "Het", "Heterozygotes", "A", "Missense")
pB <- create_plot(all_metrics, "Het", "Heterozygotes", "B", "LoF")
pC <- create_plot(all_metrics, "Hom_Alt", "Derived Homozygotes", "C", "Missense")
pD <- create_plot(all_metrics, "Hom_Alt", "Derived Homozygotes", "D", "LoF")
pE <- create_plot(all_metrics, "Total_Alt", "Total Derived Alleles", "E", "Missense")
pF <- create_plot(all_metrics, "Total_Alt", "Total Derived Alleles", "F", "LoF")

# 7. Join and Save the Final Combined Comparison
combined_comparison <- (pA + pB) / (pC + pD) / (pE + pF)

# Final high-resolution export
ggsave("Dama_Gazelle_MER_HighRes_Final.jpg", 
       plot = combined_comparison, 
       width = 12, 
       height = 16, 
       dpi = 600)

ggsave("Dama_Gazelle_MER_HighRes_Final.pdf", 
       plot = combined_comparison, 
       width = 12, 
       height = 16)

print("Figures with increased font sizes saved successfully.")
```

