#### Deleterious variants analyses using *VEP*
- I am using this tutorial from SMSC Workshop 2024 from my amazing instructors.

##### Introduction
- Ensembl's Variant Effect Predictor, __[VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)__, is a tool for annotating and analyzing genetic variants, with a focus on deleterious variants. 
- There are prebuild databases for thousands of species, and using the tool on these provides some extra features, such as prediction score for deleteriousness of non-synonymous sites, calculated from large sets of homologous proteins.
- However, it is possible to use the tool on novel genomes as well, using just the fasta sequence and a gene annotation file.
- This will provide one of the following _impact factors_ for all our variants:

| Impact factor | Description |
| --- | ----|
| LOW | synonymous variants |
| MODERATE | non-synonymous variants
| HIGH | non-sense variants (affect start, stop or reading frame) |
| MODIFIER | all other variants (for example intronic, UTRs, intergenic) |

- We will consider the classes 'MODERATE' and 'HIGH' as potentially deleterious. (This will not always be true, but we will never know the exact effect of all mutations, not even in model organisms, and that's just something we have to live with!)

- In this tutorial we will first annotate the SNP variants in the Dama gazelle, and then extract the classes LOW, MODERATE and HIGH to analyze them further.

##### Run VEP  
###### 1. Preparing the reference genome and annotation

- To run VEP on a non-model organism, you need to have a reference genome and an annotation file.
- The annotation file usually include gene models in GFF or GTF format.
- The annotation file further needs to be indexed with tabix.

##### Files and their paths is in:
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/VEP 
Reference genome: Dama_gazelle_hifiasm-ULONT_primary.fasta
Addra annotation: Addra_complete.genomic.gff
Mhorr annotation: Mohrr_complete.genomic.gff
VCF (biallelic SNPs): Dama_gazelle_biallelic_snps_autosomes.vcf
Individuals:
  Addra: SRR17129394, SRR17134085, SRR17134086
  Mhorr: SRR17134087, SRR17134088
- Since, we are doing VEP and using referecne genome assembly of Addra. Therefore, we will only use Addra.gff. We could have also used Mhorr.gff but for easiness we will use Addra.gff.
- VEP requires exons but my annotation has only CDS.
- We will duplicate the CDS as exons.

```
#### a) Create a directory in your scratch to work in
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
#### Prepare Addra GFF with exons, sort, compress and index
```bash
grep -v "#" Addra_complete.genomic.gff | \
sort -k1,1 -k4,4n -k5,5n -t$'\t' | \
awk -v OFS="\t" '{
    if ($3=="CDS") {
        print $0
        $3="exon"
        $8="."
        print $0
    } else { print $0 }
}' | bgzip -c > Addra.withExons.gff.gz

tabix -p gff Addra.withExons.gff.gz

```

- Do it for Mhorr gazelle
```bash
grep -v "#" /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mohrr_complete.genomic.gff | \
sort -k1,1 -k4,4n -k5,5n -t$'\t' | \
awk -v OFS="\t" '{
    if ($3=="CDS") {
        print $0
        $3="exon"
        $8="."
        print $0
    } else { print $0 }
}' | bgzip -c > /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mhorr.withExons.gff.gz

tabix -p gff /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mhorr.withExons.gff.gz

```

#### c) Prepare the vcf file
- If we want we could just use the vcf file as it is, but we can also do some more filtering.
- Here I decided to remove indels, and only keep bi-allelic sites (SNPs with two alleles, not more) with no missing data:
-  I had already done this but keeping below code jsut for reference.
-   I am using the vcf file already accounted for missingness and only keeping biallelic SNPs and also removing the INDELS.
-   It's upto you to decide if you want to include indels or not. I decided not.
- Use vcftools to filter if needed.
  
```bash
/scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_biallelic_snps_autosomes.vcf
```

### 2. Run the variant predictor
- This code should be placed in a script and run as a job.
- It took around ~30 min for the full genome using 8 cores.
- For this test chromosome it only took a couple of minutes. 
```bash
singularity exec vep.sif vep \
    --dir_cache /scratch/bistbs/Population_Genomic_Analysis/VEP/dummy_vep_cache \
    --fasta /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_hifiasm-ULONT_primary.fasta \
    --input_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_biallelic_snps_autosomes.vcf \
    --custom file=/scratch/bistbs/Population_Genomic_Analysis/VEP/Addra.withExons.gff.gz,short_name=ADDRA_EXONS,format=gff \
    --output_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt \
    --offline \
    --everything \
    --fork 8 \
    --no_stats \
    --quiet \
    --buffer_size 10000
```

##### Analyze the output

### 1. Familiarize yourself with the output
- If VEP worked, it will produce two or three output files: vep_rail.txt, vep_rail.txt_summary.html and possibly vep_rail.txt_warnings.txt. 
- The first file contains the raw output.
- The second file contains a neat summary suitable to open in a web browser (if the hydra cluster doesn't support opening this file, one option is to save the file locally on your laptop).
- (A copy of the summary file can be found [HERE](https://github.com/SmithsonianWorkshops/SMSC_2024_Conservation_Genomics/blob/main/Deleterious_variants/vep_GuamRails_summary.html)).
- The third file shows any warnings, for example if there are annotations in the gff file not recognized by VEP.

- All the information in the summary can also be found in the raw output. We can use a variety of unix tools to check it out:
##### Count the number of annotated variants as annotated by VEP 
##### (Grep -v removes the header lines: the ones that start with #)
- The total annotated variants is:9435988
```bash
 grep -v "^#" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt | wc -l 
```
##### How does this relate to the number variants in your input file?
- Count the number of unique annotated variants
- (cut -f1 extracts the first column only, uniq take only unique lines, and wc -l counts the lines)
- 19,638 variants have multiple annotations and 9416350 had unqiue variants.
```bash
grep -v "^#" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt | cut -f1 | uniq | wc -l

```
##### Many variants are annotated several times, for example if a gene has multiple transcripts, or if two genes are overlapping.
- Check what annotations are classified as having a low impact, and how many there are of each type
- (sort will sort the output in alphabetical order, and uniq -c will count all unique lines)
- 3745 synonymous_variant or low impact variant.
```bash
grep "IMPACT=LOW" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt  |cut -f7 |sort |uniq -c
```
- Now you can do the same with for the other impact types.

##### 2. Extract subsets of data
- As we are interested in the deleterious variants, we will mainly focus on the MODERATE and the HIGH impact category.
- However, to have something potentially neutral to compare with, we will keep the LOW category for a little bit longer.

##### First we extract the most severe category - nonsense variants
# Extract variants with a high predicted impact.

```bash
grep "IMPACT=HIGH" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt |cut -f2 |uniq >high_impact.pos.txt
```
- Now do the same for the non-synonymous variants or moderate effect variant.
```bash
grep "IMPACT=MODERATE" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt |cut -f2 |uniq >moderate_impact.pos.txt
```
- Are there any sites annotated as both HIGH and MODERATE? We can check this by joining the two files:
```bash
join <(sort high_impact.pos.txt) <(sort moderate_impact.pos.txt)
```
- If there are, we should remove the shared variants from the less severe impact class.
- Re-run extracting moderate variants, removing positions overlapping with high impact
- join -v1 will return all lines in file 1 not overlapping with file 2.
```bash
grep "IMPACT=MODERATE" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt |cut -f2 |uniq |join -v1 <(sort -) <(sort high_impact.pos.txt) >moderate_impact.pos.txt
```
- We can do the same for the LOW impact variants, removing any overlaps with either of the two previous files.
```bash
grep "IMPACT=LOW" /scratch/bistbs/Population_Genomic_Analysis/VEP/Dama_gazelle_vep_biallelic_snps.txt |cut -f2 |uniq |join -v1 <(sort -) <(sort high_impact.pos.txt) |join -v1 - <(sort moderate_impact.pos.txt) >low_impact.pos.txt
```

#### If we want to extract the variants from the vcf file using vcftools, we need tab separated positions files. 
- We can use the Stream EDitor sed to replace the ":" to tabs (\t).
```bash
sed -i 's/:/\t/g' low_impact.pos.txt
sed -i 's/:/\t/g' moderate_impact.pos.txt
sed -i 's/:/\t/g' high_impact.pos.txt
```

Create new vcf files with the variants we are interested in.
```bash
module unload gcc/8.5.0  
module load bio/vcftools/0.1.16
vcftools --vcf Dama_gazelle_biallelic_snps_autosomes.vcf  --positions low_impact.pos.txt --recode --out Dama_gazelle_low_impact
vcftools --vcf Dama_gazelle_biallelic_snps_autosomes.vcf  --positions moderate_impact.pos.txt --recode --out Dama_gazelle_moderate_impact
vcftools --vcf Dama_gazelle_biallelic_snps_autosomes.vcf  --positions high_impact.pos.txt --recode --out Dama_gazelle_high_impact
```

### 3. Analyze the alleles
- Before we start investigating genetic load, there is one important thing we need to think about: _Which_ of the two alleles in a site is the deleterious one??
- Vep doesn't provide us with this information.
- We know that a mutation in a certain position causes a non-synonymous variation, and that this could be harmful.
- But for all we know, it could be our reference individual who is carrying the deleterious allele, and the other individuals carrying the 'normal' variant.
- There are a few different ways to figure this out (see further below*), but for now we will assume that the REFERENCE allele is 'normal', and that the ALTERNATIVE allele is deleterious.
- A very convenient tool to count allele types and plot results is the vcfR package in R.
- As this course is not focusing on R, I'll show how we can look at genotypes and count alleles using different unix tools and vcftools. For those interested in the vcfR code, it can be found at the bottom of this page**.

#### a) Allele frequency spectrum
- A good method to see if our potentially deleterious sites are under more selective constraints than for example synonymous mutations, is to compare their allele frequency spectra. 
- With only give individuals we only have four alleles to work with, but let's give it a try!

```bash
# First we'll just look at the genotypes (remove everything else)
grep -v "##" Dama_gazelle_moderate_impact.recode.vcf | awk '{out=$1"\t"$2; for(i=10; i<=11; i++){split($i,s,":"); out=out"\t"s[1]}; print out}' |less
```
- This is a good start to just get a feeling for the data.
- Now we will use vcftools to calculate allele frequency for each site. This will create .frq output files that we can summarize into a little table.
``` bash
# Create a table header
echo "Type Number Ref_freq Alt_freq" |sed 's/ /\t/g' >SFS.txt
# Loop over the types, create a frequency table and summarize
for type in "low" "moderate" "high"
do
  vcftools --vcf Dama_gazelle_${type}_impact.recode.vcf --freq --out Dama_gazelle_${type}_impact
  tail -n+2 Dama_gazelle_${type}_impact.frq  |cut -f5,6 |sed 's/:/\t/g' | cut -f2,4 |sort |uniq -c |awk -v t=$type -v OFS="\t" '{print t,$0}' >>SFS.txt
done
# A lot of different unix tools here, including an awk script..
# Make sure you know what each step does!
```
- The file SFS.txt contains the site frequency spectra for all the three types of mutations.
- Below is some R code you can use for plotting, but you can use any tool you like (even Excel)

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


#### b) Masked and realized load
From the lectures, we recall that masked load comes from heterozygous deleterious mutations, and realized load comes from homozygous deleterious mutations. First we can just look at different genotype counts from one individual:

```bash
# First individual in column 10 (cut -f10 means extract the 10th column)
grep -v "##" Dama_gazelle_moderate_impact.recode.vcf |cut -f10 |cut -f1 -d":" |sort |uniq -c
```
To look at the next individual we can replace `cut -f10` with `cut -f11`. If we haven't noticed it before, some sites are _phased_, with a '|' between the alleles instead of the standard '/'. When we count for example heterozygous sites, we take both "0/1" and "0|1". Do you notice some striking difference between the individuals?

Now we use some more awk and unix tools to save heterozygous and homozygous alternative sites from high and moderate separately
```bash
# Create a new file with just a header
# Create a new file with header
echo -e "Type\tInd\tLoad\tNumber" > Load_table.txt

# Loop over types
for type in "moderate" "high"
do
    # Get the sample columns (all columns after the 9th in VCF)
    samples=($(grep "#CHROM" Dama_gazelle_${type}_impact.recode.vcf | cut -f10-))
    
    # Loop over each individual
    for idx in $(seq 0 $((${#samples[@]} - 1)))
    do
        ind=${samples[$idx]}
        col=$((idx + 10))  # VCF columns start at 1
        
        # Count masked (heterozygous) and realized (homozygous alt) load
        grep -v "^#" Dama_gazelle_${type}_impact.recode.vcf | \
        cut -f$col | cut -f1 -d":" | \
        awk -v OFS="\t" -v t=$type -v i=$ind -v masked=0 -v realized=0 \
        '{if($1=="0/1" || $1=="0|1"){masked++}else if($1=="1/1" || $1=="1|1"){realized++}} \
        END{print t,i,"masked",masked; print t,i,"realized",realized}' >> Load_table.txt
    done
done

```
##### This script is used for plotting the load(Masked and Realized load in Dama gazelle.
```R
# Load libraries
library(tidyverse)
library(grid)

# Load the Load table
LOAD <- read.table("Load_table.txt", header=TRUE) %>% as_tibble()

# Capitalize Load and Type labels
LOAD <- LOAD %>%
  mutate(
    Load = ifelse(Load == "masked", "Masked", "Realized"),
    Type = ifelse(Type == "high", "High", 
                  ifelse(Type == "moderate", "Moderate", Type))
  )

# Colors for Load types
load_colors <- c("Masked" = "#1b9e77", "Realized" = "#d95f02")

# Create the plot
p <- ggplot(LOAD, aes(x=Ind, y=Number, fill=Load)) +
  geom_bar(stat="identity", position=position_dodge(width=0.75), width=0.7) +
  geom_text(aes(label=Number), 
            position=position_dodge(width=0.75), 
            vjust=-0.4, size=4) +
  facet_wrap(Type ~ ., nrow=2, scales='free_y', strip.position = "top") +
  scale_fill_manual(values = load_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(size=12, face="bold"),
    axis.text.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    legend.text = element_text(size=13),
    legend.position = c(0.95, 0.95),          # top-right inside plot
    legend.justification = c("right", "top"),
    legend.direction = "vertical",            # vertical inside
    legend.key.size = unit(0.8, "lines"),
    legend.background = element_rect(fill = alpha('white', 0.6), color = NA),
    strip.text = element_text(size=14, face="bold"),
    strip.background = element_rect(fill="grey95", color=NA),
    panel.grid.major = element_line(color="grey80", linetype="dashed"),
    panel.grid.minor = element_line(color="grey90", linetype="dashed"),
    axis.ticks = element_line(size=0.8),
    panel.border = element_rect(color="black", fill=NA, size=0.8),
    panel.spacing.y = unit(6, "mm")
  ) +
  labs(x="Individual Sample", y="Number of Sites", fill="Load Type")

# Set publication-ready dimensions
fig_width <- 220  # mm
fig_height <- 260 # mm

# Save as PDF
ggsave("GeneticLoad_Figure_LoadTypeLegend_InsideTopRight.pdf", plot = p,
       width = fig_width, height = fig_height, units = "mm",
       dpi = 300, device = cairo_pdf)

# Save as high-quality JPEG
ggsave("GeneticLoad_Figure_LoadTypeLegend_InsideTopRight.jpeg", plot = p,
       width = fig_width, height = fig_height, units = "mm",
       dpi = 600)

```
- Do you see a difference between the individuals? Or a difference between 'HIGH' impact mutations and 'MODERATE' impact mutations?

- Remenber, this is just a very rough estimation of genetic load! Apart from finding the correct deleterious allele (which we ignored above), it might be necessary to account for differences in sequencing (if some individuals have more missing data, for example). 
- One way to do this is to calculate load as the number of deleterious alleles per genotyped sites.  
- This is the end of this tutorial! For the interested there are some extra information and code below.

####  We decided to find out the funcational classification of the hihg impact variants using VEP.
- Functional classification of HIGH-impact variants (VEP)
- Extract functional consequences of HIGH-impact variants

```bash
grep "IMPACT=HIGH" Dama_gazelle_vep_biallelic_snps.txt | \
cut -f7 | tr "," "\n" | \
sort | uniq -c | sort -nr > HighImpact_FunctionalClasses.txt

cat HighImpact_FunctionalClasses.txt
```
- Above result showed:
-  59 stop_gained
- 6 splice_donor_variant
- 27 splice_acceptor_variant
- 12 start_lost
- 9 stop_lost
- 2 non_coding_transcript_variant

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
---
- Variants for Moderate
---
####  We decided to find out the funcational classification of the Moderate impact variants using VEP.
- Functional classification of Moderate-impact variants (VEP)
- Extract functional consequences of Moderate-impact variants

```bash
grep "IMPACT=MODERATE" Dama_gazelle_vep_biallelic_snps.txt | \
cut -f7 | tr "," "\n" | \
sort | uniq -c | sort -nr > ModerateImpact_FunctionalClasses.txt

cat ModerateImpact_FunctionalClasses.txt

- result showed 3996 missense_variants

```
#### Next job is to see the chromosomal distribution of Moderate impact variants.
  - Count moderate-impact variants per chromosome
  - when I did cat on ModerateImpact_PerChromosome.txt, it showed at least all 17 autosomes had moderate impact variants.
```bash
grep "IMPACT=MODERATE" Dama_gazelle_vep_biallelic_snps.txt | \
cut -f2 | cut -d: -f1 | \
sort | uniq -c | sort -nr > MODERATEImpact_PerChromosome.txt
```
#### Extract the REALIZED load(Homozygous alternate) for all five individuals.
- Realized load = homozygous alternative(1/1 or 1|1).
- Convert VCF to BED (0-bed) for ROH intersection.
```bash
# Extract VCF header to get sample IDs
grep "#CHROM" Dama_gazelle_moderate_impact.recode.vcf > vcf_header.txt

# Loop over individuals (columns 10-14) to get 0-based BED
for col in 10 11 12 13 14; do
  ind=$(cut -f$col vcf_header.txt)
  
  grep -v "^#" Dama_gazelle_moderate_impact.recode.vcf | \
  awk -v c=$col '$c ~ /^1[\/|]1/' | \
  awk '{print $1"\t"$2-1"\t"$2}' > ${ind}_realized_moderate.bed
done

```
#### Intersect the realized load with ROH(All individuals).
- Intersect realized deleterious variants with ROH regions. 

```bash
# List of individual IDs
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    # 1. Generate per-individual ROH BED file
    awk -v id=$ind '$2 == id {print $4"\t"$7-1"\t"$8}' Dama_gazelle_ROH.hom > ${ind}.roh.bed

    # 2. Intersect realized moderate variants with ROH regions
    bedtools intersect \
        -a ${ind}_realized_moderate.bed \
        -b ${ind}.roh.bed \
        -u > ${ind}_realized_moderate_inROH.bed
done

```
#### Let's convert the PLINK ROH(.hom) to BED.
#### Convert all individuals from .hom to BED format.
- I got the file Dama_gazelle_ROH.hom from the folder of ROH analysis after running PLINK.
```bash
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088; do
  awk -v OFS="\t" -v ind=$ind '$2==ind {print $4, $7-1, $8}' Dama_gazelle_ROH.hom > ${ind}.roh.bed
done
```
#### Time to Intersect Realized Moderate Impact Variants with ROH
```
module load bedtools
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088; do
  bedtools intersect \
    -a ${ind}_realized_moderate.bed \
    -b ${ind}.roh.bed \
    -u > ${ind}_realized_moderate_inROH.bed
done
```
#### Let us summarize the realized load and ROH LOAD.
#### Output file
---
```bash
output="GeneticLoad_Moderate_ROH_Table.txt"

# Write header
echo -e "Individual\tSpecies\tTotal_Realized\tIn_ROH" > $output

# Loop over individuals
for ind in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    # Assign species
    species="Addra"
    [[ $ind == SRR17134087 || $ind == SRR17134088 ]] && species="Mhorr"

    # Count total realized moderate variants
    total=$(wc -l < ${ind}_realized_moderate.bed)
    # Count variants in ROH
    inroh=$(wc -l < ${ind}_realized_moderate_inROH.bed)

    # Append to the summary table
    echo -e "$ind\t$species\t$total\t$inroh" >> $output
done

echo "Summary table saved as $output"
$species\t$total\t$inroh" >> GeneticLoad_Moderate_ROH_Table.txt
done
```





#### Now, we will do Gene Ontology analysis for the high impact variants
- Step 1: Create VCFs per subspecies
- Suppose your VCF is Dama_gazelle_biallelic_snps_autosomes.vcf and samples are:
- Addra: SRR17129394, SRR17134085, SRR17134086
- Mhorr: SRR17134087, SRR17134088
- Use bcftools view to create per-subspecies VCFs:

#### Addra
```bash
echo -e "SRR17129394\nSRR17134085\nSRR17134086" > addra_samples.txt
bcftools view -S addra_samples.txt Dama_gazelle_biallelic_snps_autosomes.vcf -Oz -o addra.vcf.gz
bcftools index addra.vcf.gz
#### Mhorr
echo -e "SRR17134087\nSRR17134088" > mhorr_samples.txt
bcftools view -S mhorr_samples.txt Dama_gazelle_biallelic_snps_autosomes.vcf -Oz -o mhorr.vcf.gz
bcftools index mhorr.vcf.gz
```
- Step 2: Run VEP separately for Addra and Mhorr
#### Addra
  ```bash
  singularity exec vep.sif vep \
    --dir_cache /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/dummy_vep_cache\
    --fasta /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Dama_gazelle_hifiasm-ULONT_primary.fasta \
    --input_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/addra.vcf.gz \
    --custom file=/scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Addra.withExons.gff.gz,short_name=ADDRA_EXONS,format=gff \
    --output_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Addra_vep_biallelic_snps.txt \
    --offline \
    --everything \
    --fork 8 \
    --no_stats \
    --quiet \
    --buffer_size 10000
```
#### Mhorr
```bash
singularity exec vep.sif vep \
    --dir_cache /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/dummy_vep_cache\
    --fasta /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Dama_gazelle_hifiasm-ULONT_primary.fasta \
    --input_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mhorr.vcf.gz \
    --custom file=/scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mhorr.withExons.gff.gz,short_name=MHORR_EXONS,format=gff \
    --output_file /scratch/bistbs/Population_Genomic_Analysis/VEP/Gene_Ontology/Mhorr_vep_biallelic_snps.txt \
    --offline \
    --everything \
    --fork 8 \
    --no_stats \
    --quiet \
    --buffer_size 10000
    ```
#### Step 3: Extract high-impact genes for each subspecies
```bash
# Addra
awk -F'\t' '($14 ~ /IMPACT=HIGH/ && $4 != "-") {print $4}' Addra_vep_biallelic_snps.txt | sort | uniq > Addra_HighImpact_Genes.txt

# Mhorr
awk -F'\t' '($14 ~ /IMPACT=HIGH/ && $4 != "-") {print $4}' Mhorr_vep_biallelic_snps.txt | sort | uniq > Mhorr_HighImpact_Genes.txt

This gives you two clean gene lists ready for GO analysis: Use the name or ID of those gene list in the ShinyGo app.
```
---------
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
