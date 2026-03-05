#### Estimating the nucleotide diversity of the rallidae.

#### 1. Step 1: Generate the SAF file
- For a single diploid individual, the "unfolded" SFS will tell you the ratio of homozygous to heterozygous sites.
- You don't need -doMaf or -SNP_pval here because you aren't looking for SNPs across a population;
- you are looking for sites where your one bird is heterozygous.

#### a. Generate Site Allele Frequency (SAF) for a single individual
```bash
# Define your directories for cleaner code
REF="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Clapper_rail/2_Adapter_trimming/GCA_028554615.1_bRalCre1.1_genomic.fna"
BAM="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Clapper_rail/05_Remove_duplicates/SRR23269683_rmdup.bam"
OUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Nucleotide_diversity"

angsd \
    -i $BAM \
    -anc $REF \
    -ref $REF \
    -out $OUT_DIR/SRR23269683 \
    -GL 2 \
    -doSaf 1 \
    -minMapQ 30 \
    -minQ 30 \
    -baq 1 \
    -C 50 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -minInd 1 \
    -setMinDepth 5 \
    -setMaxDepth 80 \
    -doCounts 1
 
```
Note: For a single individual, use -anc (ancestral) pointing to your reference genome to "unfold" the SFS.
#### b.Estimate the SFS (The Heterozygosity)The SFS for a single diploid individual will have 3 bins (0, 1, and 2 derived alleles). 
#### The middle bin (1) represents the heterozygous sites.Bash# Calculate the 1D SFS

```bash

realSFS $OUT_DIR/SRR23269683.saf.idx > $OUT_DIR/SRR23269683.sfs
```

#### c. Calculate Diversity Statistics (Thetas)Once you have the SFS, you can calculate the per-site diversity ($\theta$). 
- Estimate thetas
```bash
# 1. Estimate thetas in binary format
realSFS saf2theta $OUT_DIR/SRR23269683.saf.idx \
    -sfs $OUT_DIR/SRR23269683.sfs \
    -outname $OUT_DIR/SRR23269683

# 2. Print the per-site diversity to a readable text file
thetaStat print $OUT_DIR/SRR23269683.thetas.idx > $OUT_DIR/SRR23269683_persite.theta.txt
```
