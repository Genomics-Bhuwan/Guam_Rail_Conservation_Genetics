#### To make the SNAIL PLOT; we have singularity container downloaded.

#### Step 1. Download the singularity using below for blobtools for snailplot.
```bash
singularity pull docker://quay.io/biocontainers/blobtools:1.1.1--py_0
```
#### Step 2. Map the reads to the assemlby with bbmap
- You have the reference genome assemlby of the species.
- You need to have the Illumina short-reads of the same individual.
- Use it to create a BAM file with coverage info.
```bash
cd /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools
```
#### 3. Map reads to assembly
- Make sure bbmap.sh is executable and run it using relative path.
```bash 
chmod +x bbmap/bbmap.sh

./bbmap/bbmap.sh \
  ref=bHypOws1_hifiasm.bp.p_ctg.fasta \
  in1=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/FMNH390989_1.fastq-003.gz \
  in2=/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/FMNH390989_2.fastq-002.gz \
  out=mapped.bam \
  threads=20
```


#### Step 4. Installation of the Diamond and run the hit using DIAMOND.
```bash
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.17/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

./diamond makedb -p 16 \
  --in uniprot/reference_proteomes.fasta.gz \
  --taxonmap uniprot/reference_proteomes.taxid_map \
  --taxonnodes taxdump/nodes.dmp \
  -d uniprot/reference_proteomes.dmnd
```

#### Step 5. Run the Diamond
-  Annotate contigs with taxonomic information using a protein database
```bash
./diamond blastx \
  -p 24 \
  --masking 0 \
  --very-sensitive \
  -d /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/snailplot/uniprot/reference_proteomes.dmnd \
  -q /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/bHypOws1_hifiasm.bp.p_ctg.fasta \
  -o Guam_Rail.blastx.tab \
  -f 6 qseqid staxids bitscore evalue \
  --max-target-seqs 1 \
  --evalue 1e-25
```

#### Step 6: Run BlobTools
- Combine assembly, coverage, and taxonomic info into a BlobTools database
```bash
singularity exec blobtools_1.1.1--py_0.sif \
  blobtools create \
    -i bHypOws1_hifiasm.bp.p_ctg.fasta \
    -b mapped.bam \
    -t Guam_Rail.blastx.tab \
    -o blobtools_db
```
#### Step 7. View Summary Statistics

- View blobtools output
```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools view \
  -i blobtools_db.blobDB.json \
  -o blobtools_db
```

#### Step 8. Visualize the standard BlobTools plots
```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools plot \
  -i blobtools_db.blobDB.json \
  -o blobtools_db
```
- The x axis represents the GC content of contigs. Y-axis has coverge from mapped reads and color: taxonomic assignment; size of dots: contig length.
- Output: blobtools_output.blobDB.json.

#### Step 9. Visualize the results
```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools view \
  -i blobtools_output.blobDB.json

singularity exec blobtools_1.1.1--py_0.sif blobtools plot \
  -i blobtools_output.blobDB.json
```
- Shows the table of coverage, GC and taxonomy.
- Plots the coverage vs. GC plots colored by taxonomy. 

#### Step 10. SNAIL PLOT
- Create a snailplot for a nice summary(GC, coverage, length, taxonomy)
```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools snail \
  -i blobtools_db.blobDB.json \
  -o blobtools_db/snailplot.jpeg \
  --threads 24
```
--- WHAT SNAILPLOT shows?
- Outerring: Contig length sorted: largest to smallest.
- Middle ring: coverage per contig.
- Inner ring: GC content.
- Colors: taxonomy per contig.
- Advantages over standard plots:
- Great for large genome and easy to spot the contamination.
- Visually summarizes GC + Coverage + taxonomy in one figure.


#### Step 11. Optional cleanup analysis.
- Filter the contaminant contigs.

  ```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools filter \
  -i blobtools_db.blobDB.json \
  -o blobtools_db_filtered \
  --keep "superkingdom:Chordata"  # Keep only vertebrate contigs

```

#### Step 12. Re-plot after filtering to check the final assembly.
- Export a table of contigs with taxon, coverage, and length for further analysis.
```bash
singularity exec blobtools_1.1.1--py_0.sif blobtools view \
  -i blobtools_db.blobDB.json \
  -o blobtools_db/export
```


--- Thanks---
