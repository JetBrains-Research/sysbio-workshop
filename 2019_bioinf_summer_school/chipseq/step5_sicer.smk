from pipeline_util import *

######## Step: Peak Calling: SICER ##################
def sicer_all_peaks_input():
    files = []

    window_size=config['sicer_window']
    gap=config['sicer_gap']

    # XXX: change significance only via config, sicer rule takes the value from
    # config, not via wildcards:
    significance=config['sicer_significance']

    sample_2_config_dict = sample_2_control(config)
    for sample in fastq_aligned_names(config):
        if sample_2_config_dict[sample] is None:
            # w/o control
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-E{significance}.scoreisland')
        else:
            # with control
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-islands-summary-FDR{significance}')

    return files

rule step5_sicer_results:
    input: *sicer_all_peaks_input()


rule bam_to_pileup:
    input: 'bams/{sample}.bam'
    output: temp('bams/pileup/{sample}.bed')

    conda: 'envs/bio.env.yaml'
    shell: 'bedtools bamtobed -i {input} > {output}'


rule pileup_bed_effective_genome_fraction:
    input:
        pileup_bed=rules.bam_to_pileup.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        temp(str(rules.bam_to_pileup.output) + ".egf")
    run:
        value = effective_genome_fraction(
            config['genome'], input.chrom_sizes, input.pileup_bed
        )
        shell("echo -n '{value}' > {output}")


def sicer_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control_pileup'] = f'bams/pileup/{control_sample}.bed'

    return dict(
        # pileup_bed='bams/pileup/{sample}.bed',
        signal_pileup = f'bams/pileup/{sample}.bed',
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
        effective_genome_fraction=rules.pileup_bed_effective_genome_fraction.output,
    )


rule call_peaks_sicer:
    input: unpack(sicer_input_fun)
    output: 'sicer/{sample}-W{width}-G{gap}-{any_suffix}'
    log: 'logs/sicer/{sample}-W{width}-G{gap}-{any_suffix}.log'

    conda: 'envs/py27.env.yaml'
    shadow: "shallow"
    params:
        significance=config['sicer_significance'],
        signal_pileup_bed_fname=lambda wildcards, input: os.path.basename(input.signal_pileup),
        control_arg=lambda wildcards, input: os.path.basename(input.control_pileup) if 'control_pileup' in input else "",
        pileups_dir=lambda wildcards, input: os.path.split(str(input.signal_pileup))[0],
        peaks_file=lambda wildcards, output: os.path.basename(output[0]),
        fragment=config['sicer_fragment'],
        genome=config['genome'],
        script=lambda wildcards, input: "SICER.sh" if 'control_pileup' in input else "SICER-rb.sh"
    shell:
        # SICER.sh ["InputDir"] ["bed file"] ["control file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] [“FDR"]
        #
        # SICER-rb.sh ["InputDir"] ["bed file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] ["E-value"]
        'mkdir -p tmp_sicer &&'
        ' echo "Significance threshold: " {params.significance} 2>&1 >> ../{log} &&'
        ' cd tmp_sicer && '
        '  {params.script} ../{params.pileups_dir} {params.signal_pileup_bed_fname} {params.control_arg}'
        '    $(pwd) {params.genome} 1 {wildcards.width}'
        '    {params.fragment} $(cat "../{input.effective_genome_fraction}")'
        '    {wildcards.gap} {params.significance} 2>&1 >> ../{log} &&'
        ' mv {params.peaks_file} ../{output} 2>&1 >> ../{log}'

# rule call_peaks_sicer0:
#     input: unpack(sicer_input_fun)
#     output: 'sicer/{sample}-W{width}-G{gap}-islands-summary-FDR{fdr}'
#     log: 'logs/sicer/{sample}-W{width}-G{gap}-islands-summary-FDR{fdr}.log'
#     # output: 'sicer/{sample}-W{width}-G{gap}-E{escore}.scoreisland'
#     # log: 'logs/sicer/{sample}-W{width}-G{gap}-E{escore}_sicer.log'
#
#     conda: 'envs/py27.env.yaml'
#     shadow: "shallow"
#     params:
#         signal_pileup_bed_fname=lambda wildcards, input: os.path.basename(input.signal_pileup),
#         control_arg=lambda wildcards, input: os.path.basename(input.control_pileup) if 'control_pileup' in input else "",
#         pileups_dir=lambda wildcards, input: os.path.split(str(input.signal_pileup))[0],
#         peaks_file=lambda wildcards, output: os.path.basename(output[0]),
#         fragment=config['sicer_fragment'],
#         genome=config['genome'],
#         script=lambda wildcards, input: "SICER.sh" if 'control_pileup' in input else "SICER-rb.sh"
#     shell:
#         # SICER.sh ["InputDir"] ["bed file"] ["control file"]
#         #       ["OutputDir"] ["Species"] ["redundancy threshold"]
#         #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
#         #       ["gap size (bp)"] [“FDR"]
#         #
#         # SICER-rb.sh ["InputDir"] ["bed file"]
#         #       ["OutputDir"] ["Species"] ["redundancy threshold"]
#         #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
#         #       ["gap size (bp)"] ["E-value"]
#         'mkdir -p tmp_sicer &&'
#         ' cd tmp_sicer && '
#         '  {params.script} ../{params.pileups_dir} {params.signal_pileup_bed_fname} {params.control_arg}'
#         '    $(pwd) {params.genome} 1 {wildcards.width}'
#         '    {params.fragment} $(cat "../{input.effective_genome_fraction}")'
#         '    {wildcards.gap} {wildcards.escore} &> ../{log} &&'
#         ' mv {params.peaks_file} ../{output} &>> ../{log}'