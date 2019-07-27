# About

Materials for ChIP-Seq processing Practice, see [2019 Bioinf Summer School]( https://bioinf.me/education/summer/2019/program)

# Installation:
## Update Conda
```bash
conda update -n base -c defaults conda
```

## Create Conda Environment
```bash
cd pipeline
conda env create --name snakemake --file environment.yaml
```

```bash
conda activate snakemake
```

# Run Pipeline
```bash
conda activate snakemake
snakemake --use-conda --pr --ri --cores 4
```