"""=================================================================================================
map stringtie transcripts to genome GFF file

Michael Gribskov     13 February 2024
================================================================================================="""
import sys
from gff import Gff, Dotdict


def seq_begin_sorter(data):
    """---------------------------------------------------------------------------------------------
    generator for gff object data ordered by sequence_id and begin position

    :param data: gff object     gff data to sort
    :return: int for            index postion in gff.data
    ---------------------------------------------------------------------------------------------"""
    data = data.data
    for entry in sorted(data, key=lambda l: (l['sequence'], l['begin'])):
        yield Dotdict(entry)

    return


def advance_sequence(ggen, sgen):
    """---------------------------------------------------------------------------------------------
    advance the pointer that is behind until both are on the chromosome. Assumes that g and s are
    seq_begin_sorter generators

    :param gptr: seq_begin_sorter (genome)
    :param sptr: seq_begin_sorter (stringtie)
    :return: bool, True if sequence records match, false otherwise
    ---------------------------------------------------------------------------------------------"""
    while g['sequence'] < s['sequence']:
        g = next(ggen)
    while g['sequence'] > s['sequence']:
        s = next(sgen)

    return True


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read both genes and transcripts into the same gff object
    # genome gff
    genome_file = 'data/phyllachora.gff'
    genome = Gff(file=genome_file, mode='GFF')
    n_genes = genome.read_feature(['gene'])
    sys.stderr.write(f'{n_genes} genes read from {genome_file}\n')

    # stringtie gtf
    stringtie_file = 'data/merged.gtf'
    genome.open(stringtie_file)
    genome.setmode('GTF')
    n_transcripts = genome.read_feature(['transcript'])
    sys.stderr.write(f'{n_transcripts} transcripts read from {stringtie_file}\n')

    # add attribute to data to store the source of the data
    genome.attribute_add('source', 'GFF', 0, n_genes)
    genome.attribute_add('source', 'GTF', n_genes, n_genes + n_transcripts)

    g_order = seq_begin_sorter(genome)
    for entry in g_order:
        print(f"{entry.sequence}\t{entry.begin}\t{entry.end}\t{entry.source}")

    exit(0)
