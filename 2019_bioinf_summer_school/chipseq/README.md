# About

Materials for ChIP-Seq processing workshop module, see [2019 Bioinf Summer School]( https://bioinf.me/education/summer/2019/program), 31 July - 1 August 2019

By Roman Chernyatchik(roman.chernyatchik@jetbrains.com) and Oleg Shpynov(oleg.shpynov@jetbrains.com).


This is simplified version of [chipseq-smk-pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) adopted for workshop demonstration.

# Installation:

The only tool required to launch the pipeline is `conda`.
* If `conda` is not installed,
follow the instructions at
[Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* If is installed, it is recommended to uptated it to the latest release:
```bash
conda update -n base -c defaults conda
```
Create a Conda environment for `snakemake`:
```bash
cd pipeline
conda env create --name snakemake --file environment.yaml
```

```bash
conda activate snakemake
```

# Development environment

For Snakemake files editing we recommend to use PyCharm IDE 

* Download latest [PyCharm Community Edition](https://www.jetbrains.com/pycharm/download)
* Enable `Snakemake` support using SnakeCharm plugin, see installation details at [SnakeCharm Plugin Home Page](https://github.com/JetBrains-Research/snakecharm).

# Run Pipeline
```bash
conda activate snakemake
snakemake --use-conda -pr --ri --cores 4
```

# Downstream & Playground

Install conda environment
```bash
conda env create --name bio --file ./downstream.yaml
```

see `commands.txt` file