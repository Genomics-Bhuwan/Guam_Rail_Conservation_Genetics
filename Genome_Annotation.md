- The whole tutorial is about Genome annotation using multiple methods. In doing so, we could compare the predcited number of genes.

#### Method 1.  Genome Annotation using *SwissProt*
```bash
#!/bin/bash -l
# Submit with: sbatch slurm_Guam.sh
#SBATCH --time=100:00:00          # Maximum runtime (5 days)
#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=20  
#SBATCH --mem=90G                
#SBATCH --partition=batch         
#SBATCH --mail-type=BEGIN,END     
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=SwissProt
#SBATCH --output=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Annotation_SwissProt/Guam_blastp_%j.log
#SBATCH --error=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Annotation_SwissProt/Guam_blastp_%j.err

# Purge old modules
module purge

# Load required modules
module load java-20
module load blast-2.13.0+

# Verify blastp path
echo "Using blastp at: $(which blastp)"

# Set InterProScan directory
IPR_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Annotation_SwissProt/

# Input protein file
Guam_PROT=$IPR_DIR/Guam_Bhuwan_proteins.fasta

# Output directory
Guam_OUT=$IPR_DIR/Guam_output
mkdir -p $Guam_OUT

# Run SwissProt Database
blastp -query $Guam_PROT \
       -db $IPR_DIR/swissprot_db \
       -evalue 1e-5 \
       -out $Guam_OUT/Guam_vs_swissprot.tsv \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads 20
```

#### Method 2.  Genome Annotation using *TrEMBL*

```bash
#!/bin/bash -l
#SBATCH --job-name=Guam_BLAST
#SBATCH --output=slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=128G
#SBATCH --time=300:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu

# Load BLAST module
module load blast/2.13.0+

# Working directory
WORK_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Annotation_TrEMBL

# Input protein file
GUAM_PROT=$WORK_DIR/Guam_Bhuwan_proteins.fasta

# Output directory
GUAM_OUT=$WORK_DIR/Guam_output
mkdir -p $GUAM_OUT

# Debug info
echo "Modules loaded:"
module list
echo "BLASTP path:"
which blastp
blastp -version

# Run BLASTP for Guam
echo "Running BLASTP for Guam..."
blastp -query $GUAM_PROT \
        -db $WORK_DIR/trembl_db \
        -evalue 1e-5 \
        -out $GUAM_OUT/Guam_vs_trembl.tsv \
        -outfmt 6 \
        -max_target_seqs 5 \
        -num_threads 24

echo "All BLAST jobs completed."
```
#### Method 3.  Genome Annotation using *InterProScan*

```bash
#!/bin/bash -l
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=InterProScan

module load java-20

IPR_DIR=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Anntation_InterProScan/interproscan-5.76-107.0
ADDRA_PROT=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Anntation_InterProScan/Guam_Bhuwan_proteins.fasta
ADDRA_OUT=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/GeMoMa/Gene_Anntation_InterProScan/Addra_output
mkdir -p $Guam_OUT

# Run InterProScan
$IPR_DIR/interproscan.sh -i $Guam_PROT \
                         -f tsv \
                         -dp \
                         -goterms \
                         -pa \
                         -cpu 24 \
                         -o $Guam_OUT/Guam_interproscan.out

```
