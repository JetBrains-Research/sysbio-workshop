from pipeline_util import *

####### Step: RAW Reads QC ##################
rule step1_raw_qc_results:
    input:
         multiqc_fastq='multiqc/fastqc/multiqc.html',

rule fastqc:
    input: config['fastq_dir'] + '/{sample}.fastq.gz'
    output:
          html='qc/fastqc/{sample}_fastqc.html',
          zip='qc/fastqc/{sample}_fastqc.zip'
    log: 'logs/fastqc/{sample}.log'

    wrapper: '0.36.0/bio/fastqc' # https://bitbucket.org/snakemake/snakemake-wrappers/src/0.31.1/bio/fastqc/

rule multiqc_fastq:
    input: expand('qc/fastqc/{sample}_fastqc.zip', sample=fastq_names_wo_ext(config))
    output: 'multiqc/fastqc/multiqc.html'
    log: 'multiqc/fastqc/multiqc.log'

    wrapper: '0.36.0/bio/multiqc'
