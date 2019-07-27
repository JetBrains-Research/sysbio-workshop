import os
import re
from glob import glob

from match_control import find_control_for


def fastq_files(config):
    fq_dir = config['fastq_dir']
    return list(glob(os.path.join(fq_dir, '*.f*q'))) + glob(os.path.join(fq_dir, '*.f*q.gz'))


def fastq_names_wo_ext(config):
    # file name w/o ext and parent folders: supports *.fastq and *.fastq.gz
    return [_split_to_fname_and_ext(f)[0] for f in fastq_files(config)]


def fastq_aligned_names(config):
    return _paired_fastq_samples_names(config) + _single_fastq_samples_names(config)


def sample_2_control(config):
    fq_files = fastq_files(config)

    result = {}
    for fq_path in fq_files:
        ext = _split_to_fname_and_ext(fq_path)[1]
        control_path = find_control_for(fq_path, ext)

        if control_path:
            control_sample = _sample_by_fastq_file(control_path)
        else:
            control_sample = None

        sample = _sample_by_fastq_file(fq_path)
        result[sample] = control_sample

    return result


def _sample_by_fastq_file(fq_path):
    name = _split_to_fname_and_ext(fq_path)[0]

    if name[-2:] in ['_1', '_2']:
        # paired file
        return name[-2:]
    else:
        # single end file
        return name


def _split_to_fname_and_ext(path):
    # assumes file has valid ext
    # understands *.{ext} and *.{ext}.gz

    fname = os.path.basename(path)

    name, dot_ext = os.path.splitext(fname)
    if dot_ext == ".gz":
        name2, dot_ext2 = os.path.splitext(name)

        # remove first dot from ext:
        return name2, dot_ext2[1:] + dot_ext
    else:
        # remove first dot from ext:
        return name, dot_ext[1:]


def _single_fastq_samples_names(config):
    paired_samples = _paired_fastq_samples_names(config)
    return [name for name in fastq_names_wo_ext(config)
            if name[-2:] not in ['_1', '_2'] or name[:-2] not in paired_samples]


def _paired_fastq_samples_names(config):
    fq_names = set(fastq_names_wo_ext(config))
    paired_samples = [name[:-2] for name in fq_names if name[-2:] == '_1']
    return [sample for sample in paired_samples if sample + '_2' in fq_names]


def effective_genome_fraction(genome, chrom_sizes_path, pileup_bed):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""

    # Get chr names covered with reads (e.g. if data is filtered by chromosome name
    # or some chrs excluded during alignment
    chromosomes = set()
    with open(str(pileup_bed)) as f:
        for line in f:
            chr = line.split()[0]
            chromosomes.add(chr)

    # Sized of chromosomes covered with reads
    chrom_sizes = {}
    with open(str(chrom_sizes_path)) as f:
        for line in f:
            chromosome, size = line.split()
            chrom_sizes[chromosome] = int(size)

    # Normalization if not all genome chromosomes are covered
    chromosomes_length = sum([chrom_sizes.get(c, 0) for c in chromosomes])
    genome_length = sum(chrom_sizes.values())

    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return (size / genome_length) * (1.0 * chromosomes_length / genome_length)


def macs_species(genome):
    """Convert genome to macs2 species encoding"""
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))
