from pipeline_util import *

######## Step: Visualization: Reads coverage ##################
rule step3_reads_coverage_results:
    input:
         bws=expand('bw/{sample}.bw', sample=fastq_aligned_names(config)),

rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample, [^/]*}.bam.bai'

    wrapper: '0.31.1/bio/samtools/index'

rule bam2bw:
    input:
         bam='bams/{filename}.bam',
         bai='bams/{filename}.bam.bai'
    output: 'bw/{filename, [^/]*}.bw'

    conda: 'envs/deeptools.env.yaml'
    threads: 4
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output}'
