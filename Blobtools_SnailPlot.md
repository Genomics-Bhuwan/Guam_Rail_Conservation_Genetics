#### There are many method to make the SNAILPLOT.

#### Change the directory to assemlby-stats, use the perl secript to convert fasta into JSON for both assemlby as given below.
#### Installation of the repository.
```bash
#### Install the blobtoolkit using apptainer or singularity.
apptainer pull docker://genomehubs/blobtoolkit:latest
```

#### These will run blobtool snail analysis for producing the SNAILPLOT.
```bash
# Move to the directory containing your files
cd /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/
```

#### 1. Create the Guam_Rail_BlobDir
```bash
  blobtools create \
    --fasta /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/bHypOws1_hifiasm.bp.p_ctg.fasta \
    /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/Guam_Rail_BlobDir
```
#### 2. Add BUSCO Data
- This adds the assembly completeness metrics.
- For this, go to the BUSCO folder. Go to the aves_db10 inside the output folderand you will see .tsv table along with other files. We need other files as well.
```bash
apptainer exec /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/blobtoolkit_latest.sif \
  blobtools add \
    --busco /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/BUSCO_for_snailplot/full_table.tsv \
    /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/Guam_Rail_BlobDir
```

#### 3. Generate the Snail Plot PNG
- This renders the final image into your Guam_Rail_Plots folder.
```bash
apptainer exec /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/blobtoolkit_latest.sif \
  blobtools view \
    --plot \
    --view snail \
    --driver chromium \
    --out /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/Guam_Rail_Plots \
    /shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Blobtools/Guam_Rail_BlobDir
```
