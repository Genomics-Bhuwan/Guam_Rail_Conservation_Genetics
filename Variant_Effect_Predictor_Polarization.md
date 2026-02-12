################################# RAILS POLARIZING ############################
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
