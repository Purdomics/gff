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


def advance_sequence(gen):
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


def overlap(gen):
    """---------------------------------------------------------------------------------------------
    find overlapping regions in data. features that begin and end at the same base are considered to
    overlap

    :param data: list           list of Dotdict, list should be sorted
    :param start: int           row to start at
    :return: list               overlapping rows
    ---------------------------------------------------------------------------------------------"""
    new = next(gen)
    region = [new]
    stop = new.end
    sequence = new.sequence

    for new in gen:
        if new.sequence == sequence and new.begin < stop:
            # add to current overlap region
            region.append(new)
            stop = max(stop, new.end)
        else:
            # start a new region
            yield region
            region = [new]
            stop = new.end
            sequence = new.sequence

    return


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

    # add the strand to the sequence name so that positive and negative strands will
    # not be merged - kind of a kludge
    for row in genome.data:
        row['sequence'] = row['sequence'] + row['strand']

    g_order = seq_begin_sorter(genome)
    n_overlap = 0
    for group in overlap(g_order):
        print(f'overlap group {n_overlap}')
        n_overlap += 1

        for entry in group:
            if entry.source == 'GFF':
                id = entry.ID
            else:
                id = entry.transcript_id

            print(f"{entry.sequence}\t{id}\t{entry.begin}\t{entry.end}\t{entry.strand}\t{entry.source}")

    exit(0)
