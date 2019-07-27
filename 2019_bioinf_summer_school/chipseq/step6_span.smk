from pipeline_util import *

######## Step: Peak Calling: SPAN ##################
rule step6_span:
    input:
        span_peaks=expand('span/{sample}_{span_bin}_{span_fdr}_{span_gap}.peak',
            sample=fastq_aligned_names(config),
            span_bin=config['span_bin'],
            span_fdr=config['span_fdr'],
            span_gap=config['span_gap']
        ),

rule step6_span_tuned:
    input:
        span_tuned_peaks=expand('span/{sample}_{span_bin}_tuned.peak',
            span_bin=config['span_bin'],
            sample=fastq_aligned_names(config)
        )

rule download_span:
    output: 'bin/span-0.11.0.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.11.0.4882.jar'


def span_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control'] = f'bams/{control_sample}.bam'

    return dict(
        signal=f'bams/{sample}.bam',
        **control_args,
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_span:
    input: unpack(span_input_fun)
    output: 'span/{sample}_{bin}_{fdr}_{gap}.peak'
    log: 'logs/span/{sample}_{bin}_{fdr}_{gap}.log'

    conda: 'envs/java8.env.yaml'
    threads: 4
    params:
        span_params=config['span_params'],
        control_arg=lambda wildcards, input: "" if 'control' not in input else f" -c {input.control}"
    shell:
        'java -Xmx8G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
        '--peaks {output} --model span/fit/{wildcards.sample}_{wildcards.bin}.span --workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {wildcards.fdr} --gap {wildcards.gap} {params.span_params} &> {log}'

rule call_peaks_span_tuned:
    input:
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output,
        bam='bams/{sample}.bam'
    output: 'span/{sample}_{bin}_tuned.peak'
    log: 'logs/span/{sample}_{bin}_tuned.log'

    conda: 'envs/java8.env.yaml'
    threads: 4
    params:
        span_markup=config['span_markup']
    shell:
        'java -Xmx8G -jar bin/span-0.11.0.jar analyze --model span/fit/{wildcards.sample}_{wildcards.bin}.span '
        '--workdir span --threads {threads}  --labels {params.span_markup} --peaks {output} &> {log}'