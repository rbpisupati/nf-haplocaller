# Nextflow pipeline for GATK best practices, SNP calling

## Installation

This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub rbpisupati/nf-haplocaller
The pipeline is written mainly to run SNPmatch on GMI HPC mendel which is PBS system. Please change config file accordingly to run it on your system.

```bash
git clone https://github.com/rbpisupati/nf-haplocaller.git
```

## Running the pipeline

```bash
nextflow run nf-haplocaller/main.nf --reads "*bam" --fasta "TAIR10_wholeGenome.fasta" --outdir output_folder
```

Add `--notrim false` to include trimming of the reads. Also you can change the trimming parameters

## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)

## Citation
Cite the paper below if you use this pipeline.
Pisupati, R. *et al.*. Verification of *Arabidopsis* stock collections using SNPmatch, a tool for genotyping high-plexed samples.  *Nature Scientific Data*  **4**, 170184 (2017).
[doi:10.1038/sdata.2017.184](https://www.nature.com/articles/sdata2017184)
