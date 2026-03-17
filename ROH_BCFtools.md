#### Estimating the ROH using bcftools : https://samtools.github.io/bcftools/howtos/roh-calling.html
#### Species 1: South Island takahe: Using two individuals for this. I can use only also. 

#### Step 1. ROH for South Island takahe
```bash
bcftools roh \
  -G30 \
  --AF-dflt 0.4 \
  /home/bistbs/Guam_rail_ROH_Final_analysis/South_Island_takhe_biallelic/10_56X_Biallelc_South_Island_takahe/South_Island_takahe_PASS_filtered_biallelic.vcf.gz \
  -o /home/bistbs/Guam_rail_ROH_Final_analysis/South_Island_takhe_biallelic/10_56X_Biallelc_South_Island_takahe/bcftools/takahe_roh_noRec.txt
```
#### Step 2. Count the size distribution of the South Island takahe
```bash
grep "^RG" /home/bistbs/Guam_rail_ROH_Final_analysis/South_Island_takhe_biallelic/10_56X_Biallelc_South_Island_takahe/bcftools/takahe_roh_output.txt | awk '
{
    mb = $6 / 1000000;
    samp = $2;
    if (mb >= 0.1 && mb < 1) bin[samp, "1"]++;
    else if (mb >= 1 && mb < 5) bin[samp, "2"]++;
    else if (mb >= 5 && mb < 10) bin[samp, "3"]++;
    else if (mb >= 10) bin[samp, "4"]++;
    
    samples[samp] = 1;
}
END {
    printf "%-15s | %-10s | %-10s | %-10s | %-10s\n", "Sample", "0.1-1Mb", "1-5Mb", "5-10Mb", ">10Mb";
    print "---------------------------------------------------------------------------";
    for (s in samples) {
        printf "%-15s | %-10d | %-10d | %-10d | %-10d\n", s, bin[s,"1"], bin[s,"2"], bin[s,"3"], bin[s,"4"];
    }
}'
```
#### Species 2: Eurasian coot
#### Step 1. ROH
```bash
# 1. Create the output directory (if not already there)
mkdir -p /home/bistbs/Guam_rail_ROH_Final_analysis/Eurasian_coot_biallelic/bcftools

# 2. Run bcftools roh
bcftools roh \
  -G30 \
  --AF-dflt 0.4 \
  /home/bistbs/Guam_rail_ROH_Final_analysis/Eurasian_coot_biallelic/Eurasian_coot_biallelic_final.vcf.gz \
  -o /home/bistbs/Guam_rail_ROH_Final_analysis/Eurasian_coot_biallelic/bcftools/eurasian_coot_roh.txt
```

#### Step 2. Size distribution 
```bash
grep "^RG" /home/bistbs/Guam_rail_ROH_Final_analysis/Eurasian_coot_biallelic/bcftools/eurasian_coot_roh.txt | awk '
{
    mb = $6 / 1000000;
    samp = $2;
    if (mb >= 0.1 && mb < 1) bin[samp, "0.1-1"]++;
    else if (mb >= 1 && mb < 5) bin[samp, "1-5"]++;
    else if (mb >= 5 && mb < 10) bin[samp, "5-10"]++;
    else if (mb >= 10) bin[samp, "10+"]++;
    
    samples[samp] = 1;
}
END {
    printf "%-15s | %-10s | %-10s | %-10s | %-10s\n", "Sample ID", "0.1-1 Mb", "1-5 Mb", "5-10 Mb", ">10 Mb";
    print "---------------------------------------------------------------------------";
    for (s in samples) {
        printf "%-15s | %-10d | %-10d | %-10d | %-10d\n", s, bin[s,"0.1-1"], bin[s,"1-5"], bin[s,"5-10"], bin[s,"10+"];
    }
}'
```
#### Step 3. Okinawa rail
#### Step 1. ROH
```bash
# 1. Create the output directory
mkdir -p /home/bistbs/Guam_rail_ROH_Final_analysis/Okinawa_rail_final_analysis/bcftools

# 2. Run bcftools roh
# Using the optimized ALT default 0.4 
bcftools roh \
  -G30 \
  --AF-dflt 0.4 \
  /home/bistbs/Guam_rail_ROH_Final_analysis/Okinawa_rail_final_analysis/Okinawa_rail_biallelic_final.vcf.gz \
  -o /home/bistbs/Guam_rail_ROH_Final_analysis/Okinawa_rail_final_analysis/bcftools/okinawa_rail_roh.txt
```
#### Step 2. Size distribution
```bash
grep "^RG" /home/bistbs/Guam_rail_ROH_Final_analysis/Okinawa_rail_final_analysis/bcftools/okinawa_rail_roh.txt | awk '
{
    mb = $6 / 1000000;
    samp = $2;
    if (mb >= 0.1 && mb < 1) bin[samp, "0.1-1"]++;
    else if (mb >= 1 && mb < 5) bin[samp, "1-5"]++;
    else if (mb >= 5 && mb < 10) bin[samp, "5-10"]++;
    else if (mb >= 10) bin[samp, "10+"]++;
    
    samples[samp] = 1;
}
END {
    printf "%-15s | %-10s | %-10s | %-10s | %-10s\n", "Sample ID", "0.1-1 Mb", "1-5 Mb", "5-10 Mb", ">10 Mb";
    print "---------------------------------------------------------------------------";
    for (s in samples) {
        printf "%-15s | %-10d | %-10d | %-10d | %-10d\n", s, bin[s,"0.1-1"], bin[s,"1-5"], bin[s,"5-10"], bin[s,"10+"];
    }
}'
```
#### Step 4. Black rail
#### Step 1. ROH
```bash
# 1. Create the output directory
mkdir -p /home/bistbs/Guam_rail_ROH_Final_analysis/Blackrail_biallelic/bcftools

# 2. Run bcftools roh
bcftools roh \
  -G30 \
  --AF-dflt 0.4 \
  /home/bistbs/Guam_rail_ROH_Final_analysis/Blackrail_biallelic/Black_rail_biallelic_final.vcf.gz \
  -o /home/bistbs/Guam_rail_ROH_Final_analysis/Blackrail_biallelic/bcftools/black_rail_roh.txt
```

#### Step 2. Size distribution
```bash
grep "^RG" /home/bistbs/Guam_rail_ROH_Final_analysis/Blackrail_biallelic/bcftools/black_rail_roh.txt | awk '
{
    mb = $6 / 1000000;
    samp = $2;
    if (mb >= 0.1 && mb < 1) bin[samp, "0.1-1"]++;
    else if (mb >= 1 && mb < 5) bin[samp, "1-5"]++;
    else if (mb >= 5 && mb < 10) bin[samp, "5-10"]++;
    else if (mb >= 10) bin[samp, "10+"]++;
    
    samples[samp] = 1;
}
END {
    printf "%-15s | %-10s | %-10s | %-10s | %-10s\n", "Sample ID", "0.1-1 Mb", "1-5 Mb", "5-10 Mb", ">10 Mb";
    print "---------------------------------------------------------------------------";
    for (s in samples) {
        printf "%-15s | %-10d | %-10d | %-10d | %-10d\n", s, bin[s,"0.1-1"], bin[s,"1-5"], bin[s,"5-10"], bin[s,"10+"];
    }
}'
```

#### Step 5. Guam rail
#### Step 1. ROH
```bash
# 1. Create the output directory
mkdir -p /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/bcftools_roh

# 2. Run bcftools roh
bcftools roh \
  -G30 \
  --AF-dflt 0.4 \
  /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/Guam_rail_biallelic_snps.vcf.gz \
  -o /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/bcftools_roh/guam_rail_roh.txt
```
#### Step 2. Size distribution
```bash
grep "^RG" /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/rm_duplicates_BAM/GVCFs/Combined_GVCF/Genotyped/bcftools_roh/guam_rail_roh.txt | awk '
{
    mb = $6 / 1000000;
    samp = $2;
    if (mb >= 0.1 && mb < 1) bin[samp, "0.1-1"]++;
    else if (mb >= 1 && mb < 5) bin[samp, "1-5"]++;
    else if (mb >= 5 && mb < 10) bin[samp, "5-10"]++;
    else if (mb >= 10) bin[samp, "10+"]++;
    
    samples[samp] = 1;
}
END {
    printf "%-15s | %-10s | %-10s | %-10s | %-10s\n", "Sample ID", "0.1-1 Mb", "1-5 Mb", "5-10 Mb", ">10 Mb";
    print "---------------------------------------------------------------------------";
    for (s in samples) {
        printf "%-15s | %-10d | %-10d | %-10d | %-10d\n", s, bin[s,"0.1-1"], bin[s,"1-5"], bin[s,"5-10"], bin[s,"10+"];
    }
}'
```
#### Step 6. 
