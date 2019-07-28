from pipeline_util import *

######## Step: Alignment QC ##################
rule step2_alignment_results:
    input:
         bams=expand("bams/{sample}.bam", sample=fastq_aligned_names(config)),
         multiqc_bam_raw='multiqc/bam_raw/multiqc.html',
         multiqc_bam='multiqc/bam_filtered/multiqc.html',

# Indexes:
rule download_chrom_sizes:
    output: '{}.chrom.sizes'.format(config['genome'])
    shell:
        'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes'

rule download_fa:
    output: directory('fa')
    shell:
        "rsync -avz --partial --exclude='*.txt' "
        'rsync://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/chromosomes/ {output}'

rule bowtie2_index:
    input: 'fa'
    output: directory('bowtie2-index')
    conda: 'envs/bio.env.yaml'

    params:
        files_list=lambda wildcards: ','.join(glob('fa/*.fa.gz')),
        target='bowtie2-index/{genome}'.format(genome=config['genome'])
    shell: 'mkdir -p {output} && bowtie2-build {params.files_list} {params.target}'

# Align
rule bowtie2_align_single:
    input:
        sample=[config['fastq_dir'] + '/{sample}.fastq.gz'],
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam.raw")
    log: "logs/bowtie2/{sample}.log"

    threads: 4
    params:
        index=lambda wildcards, input: os.path.join(
            str(input.bowtie2_index_path), config['genome']
            #XXX config['work_dir'], str(input.bowtie2_index_path), config['genome']
        ),
        extra=''
    wrapper: "0.36.0/bio/bowtie2/align"


# Aligned bams qc
rule bam_raw_stats:
    input: 'bams/{sample}.bam.raw'
    output: 'qc/bam_raw_samtools_stats/{sample}.txt'
    wrapper: '0.36.0/bio/samtools/stats'

rule bam_raw_multiqc:
    input:
        expand(
            'logs/bowtie2/{sample}.log',
            sample=fastq_aligned_names(config)
        ),
        expand(
            'qc/bam_raw_samtools_stats/{sample}.txt',
            sample=fastq_aligned_names(config)
        )
    output: 'multiqc/bam_raw/multiqc.html'
    log: 'multiqc/bam_raw/multiqc.log'
    wrapper: '0.36.0/bio/multiqc'

# Filter aligned reads with good mapping score:
rule filter_sort_bam_single:
    input: '{anywhere}/{sample}.bam.raw'
    output: '{anywhere}/{sample}.bam'

    conda: 'envs/bio.env.yaml'
    shell:
        'samtools view -bh -q30 {input} > {output}.filtered &&'
        ' samtools sort {output}.filtered -o {output} &&'
        ' rm {output}.filtered'

# Filtered reads qc
rule bam_filtered_stats:
    input: 'bams/{sample}.bam'
    output: 'qc/bam_filtered_samtools_stats/{sample}.txt'
    wrapper: '0.36.0/bio/samtools/stats'

rule bam_filtered_multiqc:
    input:
        expand(
            'qc/bam_filtered_samtools_stats/{sample}.txt',
            sample=fastq_aligned_names(config)
        )
    output: 'multiqc/bam_filtered/multiqc.html'
    log: 'multiqc/bam_filtered/multiqc.log'
    wrapper: '0.36.0/bio/multiqc'

