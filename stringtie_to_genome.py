"""=================================================================================================
map stringtie transcripts to genome GFF file

Michael Gribskov     13 February 2024
================================================================================================="""
import sys
import datetime
import argparse
from gff import Gff, Dotdict


def process_command_line():
    """---------------------------------------------------------------------------------------------

    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Compare genome annotation (GFF) and Stringtie transcripts (GTF)',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-g', '--gff',
                    help='GFF reference genome annotation',
                    type=str,
                    default='genome.gff')
    cl.add_argument('-s', '--gtf',
                    help='GTF Stringtie transcipts',
                    type=str,
                    default='merged.gtf')
    cl.add_argument('-o', '--output',
                    help='output GFF file',
                    type=str,
                    default='combined.gtf')

    return cl.parse_args()


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


def overlap(gen):
    """---------------------------------------------------------------------------------------------
    find overlapping regions in data. features that begin and end at the same base are considered to
    overlap

    :param gen: generator   provides one line of data at a time, use seq_begin_sorter()
    :return: list           overlapping rows
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
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')

    opt = process_command_line()
    sys.stderr.write(f'stringtie_to_genome {runstart}\n')
    sys.stderr.write(f'\tGFF annotation file: {opt.gff}\n')
    sys.stderr.write(f'\tGTF strintie merged file: {opt.gtf}\n')
    sys.stderr.write(f'\tGFF output file: {opt.output}\n')

    # read both genes and transcripts into the same gff object
    # genome gff
    genome = Gff(file=opt.gff, mode='GFF')
    n_genes = genome.read_feature(['gene'])
    sys.stderr.write(f'{n_genes} genes read from {opt.gff}\n')

    # stringtie gtf
    genome.open(opt.gtf)
    genome.setmode('GTF')
    n_transcripts = genome.read_feature(['transcript'])
    sys.stderr.write(f'{n_transcripts} transcripts read from {opt.gtf}\n')

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
                seqid = entry.ID
            else:
                seqid = entry.transcript_id

            print(f"{entry.sequence}\t{seqid}\t{entry.begin}\t{entry.end}\t{entry.strand}\t{entry.source}")

    exit(0)
