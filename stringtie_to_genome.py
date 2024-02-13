"""=================================================================================================
map stringtie transcripts to genome GFF file

Michael Gribskov     13 February 2024
================================================================================================="""
import sys
from gff import Gff


def seq_begin_sorter(data):
    """---------------------------------------------------------------------------------------------
    generator for gff object data ordered by sequence_id and begin position

    :param data: gff object     gff data to sort
    :return: int for            index postion in gff.data
    ---------------------------------------------------------------------------------------------"""
    data = data.data
    for entry in sorted(data, key=lambda l: (l['sequence'], l['begin'])):
        yield entry

    return


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
    for entry in s_order:
        print(f"{entry['sequence']}\t{entry['begin']}\t{entry['end']}")

    exit(0)
