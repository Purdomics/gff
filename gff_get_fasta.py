"""=================================================================================================
get sequence ranges in a gff from the corresponding Fasta file

Michael Gribskov     15 February 2024
================================================================================================="""
import sys
import datetime
import argparse
from sequence.fasta import Fasta
from gff import Gff


def process_command_line():
    """---------------------------------------------------------------------------------------------

    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Get Fasta sequences based on GFF annotation',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-g', '--gff',
                    help='GFF reference genome annotation',
                    type=str,
                    default='genome.gff')
    cl.add_argument('-f', '--fasta',
                    help='Fasta formatted reference sequence',
                    type=str,
                    default='genome.fa')
    cl.add_argument('-o', '--output',
                    help='output GFF file',
                    type=str,
                    default='extracted.fa')

    return cl.parse_args()


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')

    opt = process_command_line()
    sys.stderr.write(f'gff_get_fasta {runstart}\n')
    sys.stderr.write(f'\tGFF annotation file: {opt.gff}\n')
    sys.stderr.write(f'\tFasta sequence file: {opt.fasta}\n')
    sys.stderr.write(f'\tFasta output file: {opt.output}\n\n')

    gff = Gff(file=opt.gff, mode='GFF')
    gff.read_feature(['gene'])

    fasta = Fasta(filename=opt.fasta)
    out = Fasta(filename=opt.output, mode='w')
    n_ext = 0
    for id, doc, seq in fasta:
        print(f'{id}')
        for gene in gff.get_by_sequence(id):
            print(f'{gene.sequence}\t{gene.begin}\t{gene.end}\t{gene.attribute}')
            out.id = f'ext_{n_ext}'
            out.doc = gene.attribute
            out.seq = seq[gene.begin-1:gene.end]
            out.fh.write(f'{out.format()}\n')
            n_ext += 1

    exit(0)
