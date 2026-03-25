
# 030_snakemake_snpcalling



Snakemake workflow for SNP calling and annotation on the TLE66 normal and tumor FASTQ files.



## Main output

- `100.final/snps.annotated.tsv`



## Run command

export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2026/bin:$PATH

cd $VSC_SCRATCH/030_snakemake_snpcalling_work

snakemake --snakefile $VSC_DATA/r1031919_Giang_Lam_Nguyen/030_snakemake_snpcalling/Snakefile --configfile $VSC_DATA/r1031919_Giang_Lam_Nguyen/030_snakemake_snpcalling/config.yaml --cores 2 -p

