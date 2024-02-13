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
    # genome gff
    genome_file = 'data/phyllachora.gff'
    genome = Gff(file=genome_file, mode='GFF')
    n_genes = genome.read_feature(['gene'])
    sys.stderr.write(f'{n_genes} genes read from {genome_file}\n')

    # stringtie gtf
    stringtie_file = 'data/merged.gtf'
    transcripts = Gff(file=stringtie_file, mode='GTF')
    n_transcripts = transcripts.read_feature(['transcript'])
    sys.stderr.write(f'{n_transcripts} transcripts read from {stringtie_file}\n')

    s_order = seq_begin_sorter(transcripts)
    # for entry in s_order:
    #     print(f"{entry.sequence}\t{entry.begin}\t{entry.end}")

    g_order = seq_begin_sorter(genome)
    # for entry in g_order:
    #     print(f"{entry.sequence}\t{entry.begin}\t{entry.end}")

    # for each stringtie transcript, advance the genome pointer until the end coordinate is past the stringtie-begin
    # coordinate

    s = next(s_order)
    g = next(g_order)
    while True:
        if g.sequence < s.sequence:
            print(f"    gs   {s.sequence}\t{s.begin}\t{s.end}\t\t{g.sequence}\t{g.begin}\t{g.end}")
            g = next(g_order)
            continue
        elif g.sequence > s.sequence:
            print(f"    ss   {s.sequence}\t{s.begin}\t{s.end}\t\t{g.sequence}\t{g.begin}\t{g.end}")
            s = next(s_order)
            continue

        # if you reach here its g and s are on the same sequence

        if g.begin > s.end:
            print(f"    sa   {s.sequence}\t{s.begin}\t{s.end}\t\t{g.sequence}\t{g.begin}\t{g.end}")
            s = next(s_order)
            continue
        elif g.end < s.begin:
            print(f"    ga   {s.sequence}\t{s.begin}\t{s.end}\t\t{g.sequence}\t{g.begin}\t{g.end}")
            g = next(g_order)
            continue
        else:
            print(f"ovverlap {s.sequence}\t{s.begin}\t{s.end}\t\t{g.sequence}\t{g.begin}\t{g.end}")



    exit(0)
