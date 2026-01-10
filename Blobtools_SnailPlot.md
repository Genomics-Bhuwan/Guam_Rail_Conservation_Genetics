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

# Create the database
apptainer exec -B /shared:/shared blobtoolkit_latest.sif blobtools create \
    --fasta bHypOws1_hifiasm.bp.p_ctg.fasta \
    Guam_Rail_BlobDir
```
