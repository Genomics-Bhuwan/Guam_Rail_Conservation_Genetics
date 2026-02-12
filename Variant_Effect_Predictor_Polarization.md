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
for filt in "mac1" "mac2"
do
  for set in "GuamRail.autosomes"
  do
    s="$set.$filt"
    # Combine the two outgroups into one table
    paste bed/$s.bed outgroup/VirginiaRail.$s.weightedAF.DP$d.addN.txt \
                   outgroup/CommonMoorhen.$s.weightedAF.DP$d.addN.txt > outgroup/2outgroups_rails.$s.DP$d.txt

    # Run the assignment script
    # It checks if Virginia Rail and Moorhen agree
    python3 $scrdir/assign_ancestral.py -i outgroup/2outgroups_rails.$s.DP$d.txt -c 4 -o outgroup/Ancestral.2outgroups.$s.DP$d
    
    # Keep only sites where both outgroups agree 100%
    awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6==2){print $1,$2,$3,$4}' \
    outgroup/Ancestral.2outgroups.$s.DP$d.ancestral.txt > outgroup/Pol.Rails.$s.bed
  done
done
