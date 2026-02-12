################################# RAILS POLARIZING ############################
#### Step of variant filtration
```bash
Variant Filtering and Depth Thresholds

To ensure high-quality variant calls and reduce the inclusion of duplicated or paralogous regions, we applied hard filtering criteria using GATK v4.1.2.0.

First, joint genotyping was performed to generate a multi-sample VCF file.

We then calculated the genome-wide mean site coverage using the INFO/DP field in the VCF file with bcftools. The mean depth across all variant sites was 81.94×.

Following established best practices, we retained only sites with an overall coverage between 10× and twice the genome-wide mean coverage.

Therefore, variants with site-level depth (INFO/DP) <10× or >164× were excluded.

```
# Using Virginia Rail (Near) and Common Moorhen (Distant)
# 1. Assign Genotypes (Pseudo-haploidized)
# This uses the python script from the wolf project to pick 1 random read per site
d=5 # Minimum depth
for filt in "mac1" "mac2"
do
  for set in "GuamRail.autosomes"
  do
    s="$set.$filt"
    for i in "VirginiaRail" "CommonMoorhen"
    do
      # Pick one random allele based on depth weight
      python3 $scrdir/pseudo_haploidize.py -v outgroup/$i.$s.extract.vcf -d $d -o outgroup/$i.$s.weightedAF.DP$d
      
      # Fill missing sites with 'N'
      intersectBed -wao -a bed/$s.bed -b outgroup/$i.$s.weightedAF.DP$d.bed | \
      awk -v OFS="\t" '{if($8==0){$7="N"}; print}' | cut -f 7 > outgroup/$i.$s.weightedAF.DP$d.addN.txt
    done
  done
done

# 2. Combine and Define Ancestral State
```bash
#!/bin/bash
#SBATCH --job-name=rallidae_pipeline
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=120G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@domain.com
#SBATCH --output=logs/rallidae_pipeline_%j.out
#SBATCH --error=logs/rallidae_pipeline_%j.err

# ==============================================================================
# RALLIDAE GENOMIC PIPELINE
# ==============================================================================

# Directories
SCRDIR="./rallidae/scripts"          # custom scripts
REFDIR="./rallidae/reference"        # reference genome
VCFDIR="./rallidae/vcf"              # input VCFs
BEDDIR="./rallidae/bed"              # BED files
OUTDIR="./rallidae/outgroup"         # outgroup files
VEPDIR="./rallidae/vep"              # VEP output
AA_PROP_DIR="./rallidae/aa_prop"     # amino acid impact
LOGDIR="./rallidae/logs"

mkdir -p $BEDDIR $OUTDIR $VEPDIR $AA_PROP_DIR $LOGDIR

# Samples
FOCAL_PREFIX="Rallidae_Focal"        # your main dataset
X_PREFIX="Rallidae_Females"          # for X chromosome analysis

# ==============================================================================
# STEP 1: DOWNLOAD OUTGROUP FASTQ FILES
# ==============================================================================

# Outgroups: must be from different branches
# Example: Virginia Rail + Coot
# SRA_accession_outgroups.txt should have: SpeciesName \t SRA_ID
mkdir -p fastq
for ind in $(cut -f2 help_files/SRA_accession_outgroups.txt)
do
    echo "Downloading $ind"
    sbatch -p core -t 10:00:00 $SCRDIR/run_downloadFastq.sh $ind help_files/SRA_accession_outgroups.txt
done

# ==============================================================================
# STEP 2: MAP OUTGROUP READS TO REFERENCE
# ==============================================================================

# Map reads to your reference genome using bwa mem
for ind in $(cut -f2 help_files/SRA_accession_outgroups.txt)
do
    sbatch -p core -t 12:00:00 $SCRDIR/run_bwa_mem_outgroup.sh $REFDIR/rallidae_ref.fasta fastq/$ind.fastq.gz $OUTDIR/$ind
done

# ==============================================================================
# STEP 3: CALL SNPs AND GENERATE WEIGHTED ALLELE FILE
# ==============================================================================

# For each outgroup:
#  - Generate VCF with GATK HaplotypeCaller
#  - Randomly pick one allele per site using DP (to reduce reference bias)
D_MIN=5   # minimum depth for random sampling
for ind in $(cut -f2 help_files/SRA_accession_outgroups.txt)
do
    sbatch -p core -t 6:00:00 $SCRDIR/run_outgroup_genotype_randomDP.sh \
          -v $OUTDIR/$ind.raw.vcf -d $D_MIN -o $OUTDIR/$ind.weightedAF.DP$D_MIN
done

# Add Ns for missing sites:
# Intersect with BED of focal SNPs
for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        for ind in VirginiaRail Coot
        do
            intersectBed -wao -a $BEDDIR/$s.bed -b $OUTDIR/$ind.$s.weightedAF.DP$D_MIN.bed \
            | awk -v OFS="\t" '{if($8==0){$7="N"}; print}' | cut -f7 \
            > $OUTDIR/$ind.$s.weightedAF.DP$D_MIN.addN.txt
        done
    done
done

# ==============================================================================
# STEP 4: ASSIGN ANCESTRAL ALLELES USING TWO OUTGROUPS
# ==============================================================================

# Combine two outgroups per site and assign ancestral
for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        paste $BEDDIR/$s.bed \
              $OUTDIR/VirginiaRail.$s.weightedAF.DP$D_MIN.addN.txt \
              $OUTDIR/Coot.$s.weightedAF.DP$D_MIN.addN.txt \
              > $OUTDIR/2outgroups.$s.DP$D_MIN.txt

        # Assign ancestral allele: Python script (Linnea's method)
        python3 $SCRDIR/assign_ancestral.py \
                -i $OUTDIR/2outgroups.$s.DP$D_MIN.txt \
                -c 4 \
                -o $OUTDIR/Ancestral.2outgroups.$s.DP$D_MIN \
                2> $LOGDIR/stderr.2out.$s.DP$D_MIN.txt

        # Extract sites with 100% agreement between outgroups
        awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6==2){print $1,$2,$3,$4}' \
            $OUTDIR/Ancestral.2outgroups.$s.DP$D_MIN.ancestral.txt \
            > $OUTDIR/Pol.2out.$s.bed

        # Also as CHR:POS BASE for downstream scripts
        awk '($4!="N" && $5=="1.00" && $6==2){print $1":"$3"\t"$4}' \
            $OUTDIR/Ancestral.2outgroups.$s.DP$D_MIN.ancestral.txt \
            > $OUTDIR/Pol.2out.$s.txt
    done
done

# ==============================================================================
# STEP 5: VEP ANNOTATION AND CATEGORIES
# ==============================================================================

for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        extset="$chrset.SNPs.HF.$filt"

        # Extract VEP categories
        # Modifier
        grep "IMPACT=MODIFIER" $VEPDIR/$extset.txt \
            | grep -v "intergenic_variant" \
            | grep -v "intron" \
            | grep -v "downstream" \
            | grep -v "upstream" \
            | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' \
            | uniq > $VEPDIR/$s.firstExtract.modifier.bed

        # Synonymous
        grep "synonymous_variant" $VEPDIR/$extset.txt \
            | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' \
            | uniq > $VEPDIR/$s.firstExtract.synonymous.bed

        # Missense tolerated/deleterious
        grep "missense_variant" $VEPDIR/$extset.txt | grep "tolerated" \
            | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' \
            | uniq > $VEPDIR/$s.firstExtract.tolerated.missense.bed

        grep "missense_variant" $VEPDIR/$extset.txt | grep "deleterious" \
            | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' \
            | uniq > $VEPDIR/$s.firstExtract.deleterious.missense.bed

        # Nonsense
        grep "IMPACT=HIGH" $VEPDIR/$extset.txt \
            | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' \
            | uniq > $VEPDIR/$s.firstExtract.nonsense.bed
    done
done

# ==============================================================================
# STEP 6: MERGE CATEGORIES AND FILTER CHROMOSOMES
# ==============================================================================

for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        pref=`echo $chrset | cut -f1 -d"."`

        # Combine to single VEP final bed
        cat $VEPDIR/$s.firstExtract.synonymous.bed \
            $VEPDIR/$s.firstExtract.missense.bed \
            $VEPDIR/$s.firstExtract.nonsense.bed | sort -k1,1 -k2,2n \
            > $BEDDIR/$s.vepfinal.bed
    done
done

# ==============================================================================
# STEP 7: EXTRACT AND POLARIZE VEP SITES
# ==============================================================================

ANC="Pol.2out"
for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        pref=`echo $chrset | cut -f1 -d"."`

        # Extract VEP sites from genome-wide VCF
        intersectBed -header \
            -a $VCFDIR/$pref.SNPs.HF.$filt.vcf.gz \
            -b $BEDDIR/$s.vepfinal.bed \
            > $VCFDIR/$s.vepfinal.vcf

        # Polarize with ancestral allele
        perl $SCRDIR/addAAInfoToVCF_fromList.pl \
            $VCFDIR/$s.vepfinal.vcf \
            $OUTDIR/$ANC.$s.txt \
            $VCFDIR/$ANC.$filt/$s.vepfinal.vcf
    done
done

# ==============================================================================
# STEP 8: CALCULATE DERIVED ALLELE COUNTS AND INDIVIDUAL MUTATION LOAD
# ==============================================================================

for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        python3 $SCRDIR/count_derived_alleles.py \
            -v $VCFDIR/$ANC.$filt/$s.vepfinal.vcf \
            -o $VCFDIR/$ANC.$filt/$s.derived_counts.txt \
            -c VEP,SIFT
    done
done

# ==============================================================================
# STEP 9: OPTIONAL: ROH vs DELETERIOUS BURDEN
# ==============================================================================

# If you have ROH files per individual:
# - Intersect ROH with deleterious sites (missense + nonsense)
# - Summarize mutation load per ROH region

for filt in "mac1" "mac2"
do
    for chrset in "$X_PREFIX.chrX" "$FOCAL_PREFIX.chr1-38"
    do
        s="$chrset.$filt"
        python3 $SCRDIR/roh_deleterious_burden.py \
            -roh help_files/ROH.individuals.txt \
            -deleterious $BEDDIR/$s.vepfinal.bed \
            -o $VCFDIR/$ANC.$filt/$s.roh_burden.txt
    done
done

echo "Pipeline completed!"
```
