# About

Materials for ChIP-Seq processing Practice, see [2019 Bioinf Summer School]( https://bioinf.me/education/summer/2019/program)

This is simplified version of [chipseq-smk-pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) adopted for workshop demonstration.

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

## Development environment

For editing Snakemake files we recommend to use PyCharm IDE 

* Download latest [PyCharm Community Edition](https://www.jetbrains.com/pycharm/download)
* Enable `Snakemake` support using SnakeCharm plugin, see installation details at [SnakeCharm Plugin Home Page](https://github.com/JetBrains-Research/snakecharm).

# Run Pipeline
```bash
conda activate snakemake
snakemake --use-conda --pr --ri --cores 4
```

# Downstream & Playground

Install conda environment
```bash
conda env create --name bio --file ./downstream.yaml
```

see `commands.txt` file