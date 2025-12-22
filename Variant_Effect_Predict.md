#### Deleterious variants analyses using *VEP*
- I am using this tutorial from SMSC Workshop 2024 from my amazing instructors.

##### Introduction
- Ensembl's Variant Effect Predictor, __[VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)__, is a tool for annotating and analyzing genetic variants, with a focus on deleterious variants. 
- There are prebuild databases for thousands of species, and using the tool on these provides some extra features, such as prediction score for deleteriousness of non-synonymous sites, calculated from large sets of homologous proteins.
- However, it is possible to use the tool on novel genomes as well, using just the fasta sequence and a gene annotation file.
- This will provide one of the following _impact factors_ for all our variants:

##### Run VEP  
###### 1. Preparing the reference genome and annotation

- To run VEP on a non-model organism, you need to have a reference genome and an annotation file.
- The annotation file usually include gene models in GFF or GTF format.
- The annotation file further needs to be indexed with tabix.

##### Files and their paths is in:
#### a) Create a directory in your directory to work in
```bash


#### b) Indexing the annotation file
Vep requires that exons are annotated. Since our annotation file lacks exons, we will duplicate the CDS annotation and just replace the "CDS" with "exon" before we index the file.
```bash
# Tabix should be included in the vep module
- I cloned ensembl-vep from:
- Install it from Singularity. It is easy as compared to manual intallation or other cause of issues with the perl and other dependencies.

##### OR if you have singularity. You can directly run using Singularity.
```bash
singularity pull docker://ensemblorg/ensembl-vep

singularity exec ../ensembl-vep_latest.sif \
  vep -i homo_sapiens_GRCh38.vcf \
      --cache --offline \
      -o homo_sapiens_GRCh38.vep.vcf

```
---bash
module load ensembl-vep
module unload gcc/8.5.0                                                                                            
module load htslib
```
#### Sort file and add exons before compressing
#### Prepare Guam_rail GFF with exons, sort, compress and index
- Taylor already did this for us.
```bash
/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail.withExons.gff.gz
```

#### c) Prepare the vcf file
- If we want we could just use the vcf file as it is, but we can also do some more filtering.
- Here I decided to remove indels, and only keep bi-allelic sites (SNPs with two alleles, not more) with no missing data:
-  I had already done this but keeping below code jsut for reference.
-   I am using the vcf file already accounted for missingness and only keeping biallelic SNPs and also removing the INDELS.
-   It's upto you to decide if you want to include indels or not. I decided not.
- Use vcftools to filter if needed.
  


### 2. Run the variant predictor
- This code should be placed in a script and run as a job.
- It took around ~30 min for the full genome using 8 cores.
- For this test chromosome it only took a couple of minutes. 
```bash
singularity exec vep.sif vep \
    --dir_cache /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/dummy_vep_cache \
    --fasta /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/bHypOws1_hifiasm.bp.p_ctg.fasta \
    --input_file /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_biallelic_snps.vcf.gz \
    --custom file=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail.withExons.gff.gz,short_name=GUAM_RAIL_EXONS,format=gff \
    --output_file /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_vep_biallelic_snps.txt \
    --offline \
    --everything \
    --fork 8 \
    --no_stats \
    --quiet \
    --buffer_size 10000
```

##### Analyze the output ################################################################################

### 1. Familiarize yourself with the output
- If VEP worked, it will produce two or three output files: vep_rail.txt, vep_rail.txt_summary.html and possibly vep_rail.txt_warnings.txt. 
- The first file contains the raw output.
- The second file contains a neat summary suitable to open in a web browser (if the hydra cluster doesn't support opening this file, one option is to save the file locally on your laptop).
- (A copy of the summary file can be found [HERE](https://github.com/SmithsonianWorkshops/SMSC_2024_Conservation_Genomics/blob/main/Deleterious_variants/vep_GuamRails_summary.html)).
- The third file shows any warnings, for example if there are annotations in the gff file not recognized by VEP.

- All the information in the summary can also be found in the raw output. We can use a variety of unix tools to check it out:
##### Count the number of annotated variants as annotated by VEP 
##### (Grep -v removes the header lines: the ones that start with #)
- The total annotated variants is:12520208
```bash
 grep -v "^#" /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_vep_biallelic_snps.txt | wc -l 
```
##### How does this relate to the number variants in your input file?
- Count the number of unique annotated variants
- (cut -f1 extracts the first column only, uniq take only unique lines, and wc -l counts the lines)
```bash
grep -v "^#" /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_vep_biallelic_snps.txt | cut -f1 | uniq | wc -l

```
##### Many variants are annotated several times, for example if a gene has multiple transcripts, or if two genes are overlapping.
- Check what annotations are classified as having a low impact, and how many there are of each type
- (sort will sort the output in alphabetical order, and uniq -c will count all unique lines)

```bash
grep "IMPACT=LOW" /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_vep_biallelic_snps.txt  |cut -f7 |sort |uniq -c
```
- Now you can do the same with for the other impact types.

##### 2. Extract subsets of data
- As we are interested in the deleterious variants, we will mainly focus on the MODERATE and the HIGH impact category.
- However, to have something potentially neutral to compare with, we will keep the LOW category for a little bit longer.
  
```bash
module unload gcc-8.5.0
module load vcf-tools
```
```bash
VEP_FILE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_vep_biallelic_snps.txt"
VCF_FILE="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Guam_rail_Population_Genomics/VEP/Guam_rail_biallelic_snps.vcf.gz"
```
#### Extract unique positions per impact
```bash
grep "IMPACT=HIGH" "$VEP_FILE" | cut -f2 | sort -u > high_all.pos.txt
grep "IMPACT=MODERATE" "$VEP_FILE" | cut -f2 | sort -u > moderate_all.pos.txt
grep "IMPACT=LOW" "$VEP_FILE" | cut -f2 | sort -u > low_all.pos.txt
```

#### Define HIGH-only (all HIGH variants)
```bash
cp high_all.pos.txt high_only.pos.txt
```
#### Define MODERATE-only (remove HIGH overlap)
```bash
join -v1 moderate_all.pos.txt high_only.pos.txt > moderate_only.pos.txt
```
#### Define LOW-only (remove HIGH and MODERATE overlap)
```bash
join -v1 low_all.pos.txt high_only.pos.txt | join -v1 - moderate_only.pos.txt > low_only.pos.txt
```
#### Convert to tab-separated format for VCF extraction
```bash
sed -i 's/:/\t/g' high_only.pos.txt
sed -i 's/:/\t/g' moderate_only.pos.txt
sed -i 's/:/\t/g' low_only.pos.txt
```
#### Extract VCF subsets using vcftools
```bash
vcftools --gzvcf "$VCF_FILE" --positions high_only.pos.txt --recode --out Guam_rail_high_impact
vcftools --gzvcf "$VCF_FILE" --positions moderate_only.pos.txt --recode --out Guam_rail_moderate_impact
vcftools --gzvcf "$VCF_FILE" --positions low_only.pos.txt --recode --out Guam_rail_low_impact
```
#### Quick visual check for moderate impact 
```bash
grep -v "##" Guam_rail_moderate_impact.recode.vcf | \
awk '{out=$1"\t"$2; for(i=10; i<=11; i++){split($i,s,":"); out=out"\t"s[1]}; print out}' | less
```
#### Calculate the allele frequencies with VCFtools
- We want to know the site frequency spectrum(SFS) for each impact class(Low, Moderate, High).

```bash
- Create SFS output table with header
echo -e "Type\tNumber\tRef_freq\tAlt_freq" > SFS.txt

- Loop over impact classes
for type in "low" "moderate" "high"
do
    - Calculate allele frequencies per site
    vcftools --vcf Guam_rail_${type}_impact.recode.vcf --freq --out Guam_rail_${type}_impact

    - Extract allele frequencies and summarize counts
    tail -n+2 Guam_rail_${type}_impact.frq | \
    cut -f5,6 | sed 's/:/\t/g' | cut -f2,4 | sort | uniq -c | \
    awk -v t=$type -v OFS="\t" '{print t,$0}' >> SFS.txt
done
```
```R
require(tidyverse)

file<-"SFS.txt"
SFS<-file %>% read.table(header=TRUE) %>% as_tibble()

# Order the types
SFS$Type<-factor(SFS$Type, levels=c("high","moderate","low"))

ggplot(SFS, aes(x=Alt_freq, y=Number, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```
What do we see here? The sizes are so different between the types so they are hard to compare! We can try making the bars using relative sizes instead:

```R
# With relative sizes
SFS_rel <- SFS %>% group_by(Type) %>% mutate(Rel=Number/sum(Number))

ggplot(SFS_rel, aes(x=Alt_freq, y=Rel, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```
Now it looks better! We see that the 'HIGH' category is shifted to the left, and the 'LOW' category has relatively more sites with higher alternative frequency. Can you explain why?

################################################################################################################################################################
#### b) Masked and realized load
From the lectures, we recall that masked load comes from heterozygous deleterious mutations, and realized load comes from homozygous deleterious mutations. First we can just look at different genotype counts from one individual:

```bash
# Create output directory
OUTPUT_DIR="./Heterozygosity_load"
mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="$OUTPUT_DIR/Load_table.txt"

# Add header
echo -e "Type\tInd\tLoad\tNumber" > "$OUTPUT_FILE"

# Loop over types
for type in "moderate" "high" "low"
do
    VCF_FILE="./Guam_rail_${type}_impact.recode.vcf"

    # Get sample names from VCF header
    samples=($(grep "^#CHROM" "$VCF_FILE" | cut -f10-))

    # Loop over each individual
    for idx in $(seq 0 $((${#samples[@]} - 1)))
    do
        ind=${samples[$idx]}
        col=$((idx + 10))  # VCF columns start at 1

        # Count heterozygous (masked) and homozygous alt (realized)
        grep -v "^#" "$VCF_FILE" | \
        cut -f$col | cut -f1 -d":" | \
        awk -v OFS="\t" -v t=$type -v i=$ind \
        '{if($1=="0/1" || $1=="0|1"){masked++} else if($1=="1/1" || $1=="1|1"){realized++}} \
        END{print t,i,"masked",masked; print t,i,"realized",realized}' >> "$OUTPUT_FILE"
    done
done

echo "Heterozygosity load table saved in: $OUTPUT_FILE"

```
##### This script is used for plotting the load(Masked and Realized load in Dama gazelle.
```R
library(tidyverse)
library(grid)

# Load the Load table
LOAD <- read.table("Load_table.txt", header=TRUE, fill=TRUE) %>% as_tibble()

# Replace missing Number values with 0
LOAD <- LOAD %>%
  mutate(Number = as.numeric(Number),
         Number = ifelse(is.na(Number), 0, Number))

# Capitalize Load and Type labels
LOAD <- LOAD %>%
  mutate(
    Load = ifelse(Load == "masked", "Masked", "Realized"),
    Type = case_when(
      Type == "high" ~ "High",
      Type == "moderate" ~ "Moderate",
      Type == "low" ~ "Low",
      TRUE ~ Type
    )
  )

# Colors for Load types
load_colors <- c("Masked" = "#1b9e77", "Realized" = "#d95f02")

# Combined plot with bounding box
p <- ggplot(LOAD, aes(x=Ind, y=Number, fill=Load)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
  geom_text(aes(label=Number), 
            position=position_dodge(width=0.8), vjust=-0.4, size=3.5) +
  facet_wrap(~Type, nrow=1, scales='free_y') +  # All types in one row
  scale_fill_manual(values = load_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(size=12, angle=45, hjust=1, face="bold"),
    axis.text.y = element_text(size=12),
    legend.title = element_text(size=14, face="bold"),
    legend.text = element_text(size=12),
    legend.position = "top",
    strip.text = element_text(size=14, face="bold"),
    panel.grid.major = element_line(color="grey80", linetype="dashed"),
    panel.grid.minor = element_line(color="grey90", linetype="dashed"),
    panel.border = element_rect(color="black", fill=NA, size=0.8), # bounding box for each panel
    plot.background = element_rect(color="black", size=1, fill=NA), # bounding box for entire plot
    panel.spacing = unit(1, "lines")
  ) +
  labs(x="Individual Sample", y="Number of Sites", fill="Load Type")

# Save combined figure with bounding box
ggsave("GeneticLoad_Combined_AllTypes_Box.pdf", plot = p, width=300, height=180, units="mm", dpi=300)
ggsave("GeneticLoad_Combined_AllTypes_Box.jpeg", plot = p, width=300, height=180, units="mm", dpi=600)

```

####  We decided to find out the funcational classification of the hihg impact variants using VEP.
- Functional classification of HIGH-impact variants (VEP)
- Extract functional consequences of HIGH-impact variants
```bash
#### Count HIGH-impact functional classes for Guam rail

grep "IMPACT=HIGH" Guam_rail_vep_biallelic_snps.txt | \
cut -f7 | tr "," "\n" | \
sort | uniq -c | sort -nr > GuamRail_HighImpact_FunctionalClasses.txt

# View the results
cat GuamRail_HighImpact_FunctionalClasses.txt

```
- Above result showed:
- stop_gained	770
- splice_donor_variant	448
- splice_acceptor_variant	190
- start_lost	184
- stop_lost	107
- splice_region_variant	46


  #### Next job is to see the chromosomal distribution of High impact variants.
  - Count HIGH-impact variants per chromosome
  - when I did cat on HighImpact_PerChromosome.txt, it showed at least all 17 autosomes had high impact variants.
```bash
grep "IMPACT=HIGH" Dama_gazelle_vep_biallelic_snps.txt | \
cut -f2 | cut -d: -f1 | \
sort | uniq -c | sort -nr > HighImpact_PerChromosome.txt
```
#### Extract the REALIZED load(Homozygous alternate) for all five individuals.
- Realized load = homozygous alternative(1/1 or 1|1).
- Convert VCF to BED (0-bed) for ROH intersection.
```bash
# Extract VCF header to get sample IDs
grep "#CHROM" Dama_gazelle_high_impact.recode.vcf > vcf_header.txt

# Loop over individuals (columns 10-14) to get 0-based BED
for col in 10 11 12 13 14; do
  ind=$(cut -f$col vcf_header.txt)
  
  grep -v "^#" Dama_gazelle_high_impact.recode.vcf | \
  awk -v c=$col '$c ~ /^1[\/|]1/' | \
  awk '{print $1"\t"$2-1"\t"$2}' > ${ind}_realized_high.bed
done

```
#### Intersect the realized load with ROH(All individuals).
- Intersect realized deleterious variants with ROH regions. 

```bash
module load bedtools-2.28

for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
  bedtools intersect \
    -a ${ind}_realized_high.bed \
    -b ${ind}.roh.bed \
    -u > ${ind}_realized_high_inROH.bed
done
```
#### Let's convert the PLINK ROH(.hom) to BED.
#### Convert all individuals from .hom to BED format.
- I got teh file Dama_gazelle_ROH.hom from the folder of ROH analysis after running PLINK.
```bash
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088; do
  awk -v OFS="\t" -v ind=$ind '$2==ind {print $4, $7-1, $8}' Dama_gazelle_ROH.hom > ${ind}.roh.bed
done
```
#### Time to Intersect Realized High Impact Variants with ROH
```
module load bedtools
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088; do
  bedtools intersect \
    -a ${ind}_realized_high.bed \
    -b ${ind}.roh.bed \
    -u > ${ind}_realized_high_inROH.bed
done
```
#### Let us summarize the realized load and ROH LOAD.
- Let's create a summary Table.
```bash
Summarize Realized Load and ROH Load
# Create summary table
echo -e "Individual\tSpecies\tTotal_Realized\tIn_ROH" >> GeneticLoad_ROH_Table.txt

for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
  # Correct species assignment
  if [[ "$ind" == "SRR17129394" || "$ind" == "SRR17134085" || "$ind" == "SRR17134086" ]]; then
      species="Addra"
  else
      species="Mhorr"
  fi

  # Count realized load and in ROH
  total=$(wc -l < ${ind}_realized_high.bed)
  inroh=$(wc -l < ${ind}_realized_high_inROH.bed)

  # Append to table
  echo -e "$ind\t$species\t$total\t$inroh" >> GeneticLoad_ROH_Table.txt
done

```
#### Create a combined summary file.
```bash
# Header
echo -e "Individual\tSpecies\tTotal_Realized\tIn_ROH" > GeneticLoad_ROH_Table.txt

# Loop over individuals
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    # Determine species
    species="Addra"
    [[ $ind == SRR17134087 || $ind == SRR17134088 ]] && species="Mhorr"

    # Count total realized load variants
    total=$(wc -l < ${ind}_realized_high.bed)
    # Count realized load in ROH
    inroh=$(wc -l < ${ind}_realized_high_inROH.bed)

    # Append to table
    echo -e "$ind\t$species\t$total\t$inroh" >> GeneticLoad_ROH_Table.txt
done
```
###############################################################################################################################

### IN ROH vs In High Impact variants.
```bash
module load bedtools-2.28
# -----------------------------
# Create unique output directory
# -----------------------------
OUTDIR="GuamRail_HighImpact_Results"
mkdir -p $OUTDIR
echo "All outputs will be saved in $OUTDIR"

# -----------------------------
# Step 1: Functional Classes of HIGH-impact variants
# -----------------------------
grep "IMPACT=HIGH" Guam_rail_vep_biallelic_snps.txt | \
cut -f7 | tr "," "\n" | \
sort | uniq -c | sort -nr > $OUTDIR/HighImpact_FunctionalClasses.txt
echo "HighImpact Functional Classes:"
cat $OUTDIR/HighImpact_FunctionalClasses.txt

# -----------------------------
# Step 2: HIGH-impact variants per chromosome
# -----------------------------
grep "IMPACT=HIGH" Guam_rail_vep_biallelic_snps.txt | \
cut -f2 | cut -d: -f1 | \
sort | uniq -c | sort -nr > $OUTDIR/HighImpact_PerChromosome.txt
echo "HighImpact variants per chromosome:"
cat $OUTDIR/HighImpact_PerChromosome.txt

# -----------------------------
# Step 3: Extract VCF header
# -----------------------------
grep "#CHROM" Guam_rail_high_impact.recode.vcf > $OUTDIR/vcf_header.txt

# -----------------------------
# Step 4: Convert VCF to 0-based BED for homozygous alternate (1/1 or 1|1)
# -----------------------------
for col in 10 11 12; do
    ind=$(cut -f$col $OUTDIR/vcf_header.txt)
    
    grep -v "^#" Guam_rail_high_impact.recode.vcf | \
    awk -v c=$col '$c ~ /^1[\/|]1/' | \
    awk '{print $1"\t"$2-1"\t"$2}' > $OUTDIR/${ind}_realized_high.bed
done

# -----------------------------
# Step 5: Convert PLINK ROH .hom to BED
# -----------------------------
for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    awk -v OFS="\t" -v ind=$ind '$2==ind {print $4, $7-1, $8}' Guam_rail_ROH.hom > $OUTDIR/${ind}.roh.bed
done

# -----------------------------
# Step 6: Intersect realized HIGH-impact variants with ROH
# -----------------------------
for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    bedtools intersect \
        -a $OUTDIR/${ind}_realized_high.bed \
        -b $OUTDIR/${ind}.roh.bed \
        -u > $OUTDIR/${ind}_realized_high_inROH.bed
done

# -----------------------------
# Step 7: Summarize Realized Load and ROH Load
# -----------------------------
echo -e "Individual\tSpecies\tTotal_Realized\tIn_ROH" > $OUTDIR/GeneticLoad_ROH_Table.txt

for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    species="Guam_rail"

    total=$(wc -l < $OUTDIR/${ind}_realized_high.bed)
    inroh=$(wc -l < $OUTDIR/${ind}_realized_high_inROH.bed)

    echo -e "$ind\t$species\t$total\t$inroh" >> $OUTDIR/GeneticLoad_ROH_Table.txt
done

echo "Summary table saved to $OUTDIR/GeneticLoad_ROH_Table.txt"
cat $OUTDIR/GeneticLoad_ROH_Table.txt
````

#### For moderate impact variants
```bash
module load bedtools-2.28

# -----------------------------
# Create unique output directory
# -----------------------------
OUTDIR="GuamRail_ModerateImpact_Results"
mkdir -p $OUTDIR
echo "All outputs will be saved in $OUTDIR"

# -----------------------------
# Step 1: Functional Classes of MODERATE-impact variants
# -----------------------------
grep "IMPACT=MODERATE" Guam_rail_vep_biallelic_snps.txt | \
cut -f7 | tr "," "\n" | \
sort | uniq -c | sort -nr > $OUTDIR/ModerateImpact_FunctionalClasses.txt
echo "ModerateImpact Functional Classes:"
cat $OUTDIR/ModerateImpact_FunctionalClasses.txt

# -----------------------------
# Step 2: MODERATE-impact variants per chromosome
# -----------------------------
grep "IMPACT=MODERATE" Guam_rail_vep_biallelic_snps.txt | \
cut -f2 | cut -d: -f1 | \
sort | uniq -c | sort -nr > $OUTDIR/ModerateImpact_PerChromosome.txt
echo "ModerateImpact variants per chromosome:"
cat $OUTDIR/ModerateImpact_PerChromosome.txt

# -----------------------------
# Step 3: Extract VCF header
# -----------------------------
grep "#CHROM" Guam_rail_moderate_impact.recode.vcf > $OUTDIR/vcf_header.txt

# -----------------------------
# Step 4: Convert VCF to 0-based BED for homozygous alternate (1/1 or 1|1)
# -----------------------------
for col in 10 11 12; do
    ind=$(cut -f$col $OUTDIR/vcf_header.txt)
    
    grep -v "^#" Guam_rail_moderate_impact.recode.vcf | \
    awk -v c=$col '$c ~ /^1[\/|]1/' | \
    awk '{print $1"\t"$2-1"\t"$2}' > $OUTDIR/${ind}_realized_moderate.bed
done

# -----------------------------
# Step 5: Convert PLINK ROH .hom to BED
# -----------------------------
for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    awk -v OFS="\t" -v ind=$ind '$2==ind {print $4, $7-1, $8}' Guam_rail_ROH.hom > $OUTDIR/${ind}.roh.bed
done

# -----------------------------
# Step 6: Intersect realized MODERATE-impact variants with ROH
# -----------------------------
for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    bedtools intersect \
        -a $OUTDIR/${ind}_realized_moderate.bed \
        -b $OUTDIR/${ind}.roh.bed \
        -u > $OUTDIR/${ind}_realized_moderate_inROH.bed
done

# -----------------------------
# Step 7: Summarize Realized Load and ROH Load
# -----------------------------
echo -e "Individual\tSpecies\tTotal_Realized\tIn_ROH" > $OUTDIR/GeneticLoad_ROH_Moderate_Table.txt

for ind in FMNH390989 HOW_N23-0063 HOW_N23-0568; do
    species="Guam_rail"

    total=$(wc -l < $OUTDIR/${ind}_realized_moderate.bed)
    inroh=$(wc -l < $OUTDIR/${ind}_realized_moderate_inROH.bed)

    echo -e "$ind\t$species\t$total\t$inroh" >> $OUTDIR/GeneticLoad_ROH_Moderate_Table.txt
done

echo "Summary table saved to $OUTDIR/GeneticLoad_ROH_Moderate_Table.txt"
cat $OUTDIR/GeneticLoad_ROH_Moderate_Table.txt
```


### *About finding out which allele is deleterious
- In a large population, deleterious variants will most likely be selected against, and will never reach high frequencies. 
- Therefore it is often safe to assume that the _minor_ allele is the deleterious variant.
- But in a small population, we know that also deleterious variants can reach high frequencies just due to drift.
- And what if we only have a couple of samples in the first place! With only 4 alleles in total, it is hard to tell which is the minor!
- Another option is to _polarize_ the variants (i.e. to determine the ancestral allele), and assume the ancestral allele is the 'normal' variant.
- This approach has been used in many conservation studies. Can you think of caveats with this method?
- Do we have enough data in the Dama gazelle project to polarize our variants?

---------


























### ** Repeat Analysis part 3. in R with the vcfR package
This code is adapted from the github [https://github.com/linneas/wolf-deleterious](https://github.com/linneas/wolf-deleterious) (Smeds et al. 2023).
When using vcfR it is more convenient to have all variants in a single vcf file:
```bash
rm -f all_three.pos.txt
for type in "low" "moderate" "high"
do
  awk -v t=$type '{print $0"\t"t}' ${type}_impact.pos.txt |uniq >>all_three.pos.txt
done
vcftools --vcf Dama_gazelle_SNPs.recode.vcf --positions all_three.pos.txt \
--recode --out Dama_gazelle_all_three
```
Read the vcf into R and perform the same analysis as above:
```R

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)

# Reading the two files
vcf_file="Dama_gazelle_all_three.recode.vcf"
type_file="all_three.pos.txt"
vcf <- read.vcfR(vcf_file)
vep_types <- type_file %>% read.table(header=FALSE) %>% as_tibble() %>%
        rename(CHROM=V1, POS=V2, Type=V3)

# Convert vcf to tidy format
tidy_vcf <- vcfR2tidy(vcf,
              #info_fields=c("AA"),   # This could be used if we have ancestral information
              format_fields=c("GT"),
              dot_is_NA=TRUE)

# Extracting relevant data from vcf and merge with VEP impact classes
tidy_gt <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% inner_join(tidy_vcf$fix) %>%
  inner_join(vep_types) %>% mutate(gt=paste(substr(gt_GT,1,1),substr(gt_GT,3,3), sep=""))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Building the SFS
SFS <- tidy_gt %>% select(CHROM, POS, Indiv, gt, Type) %>%
  group_by(CHROM, POS, gt, Type) %>% summarize(count=n()) %>%
  mutate(alt=case_when((gt=='11') ~ (count*2.0),
                       (gt=='01') ~ (count*1.0),
                       TRUE ~ 0)) %>%
  group_by(CHROM, POS, Type) %>%  summarize(totalt=sum(alt)) %>%
  group_by(Type, totalt) %>% summarize(count=n()) %>% ungroup() %>%
  group_by(Type) %>% mutate(Rel=count/sum(count))

# Order the types before plotting
SFS$Type<-factor(SFS$Type, levels=c("high","moderate","low"))

# Plot SFS
ggplot(SFS, aes(x=totalt, y=Rel, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(x="Alt alleles", y="Fraction of sites")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate the load

# Summarize the number of genotypes (skip 'low' mutations)
LOAD <- tidy_gt %>% filter(Type!="low") %>% select(Indiv, Type, gt) %>%
            group_by(Indiv,Type,gt) %>% summarize(count=n()) %>% ungroup() %>%
            mutate(Load=case_when(gt=="01" ~ "masked",
                                  gt=="11" ~ "realized",
                                  TRUE ~ NA)) %>% drop_na()

# Plot load
ggplot(LOAD, aes(x=Load, y=count, fill=Indiv)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.8) +
  facet_wrap(Type~., nrow=2, scales='free_y')
```
