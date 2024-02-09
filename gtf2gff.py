"""=================================================================================================
Should be able to do this within gff.py, but I'm in a hurry
written to convert braker output



Michael Gribskov     22 February 2021
================================================================================================="""
import sys


class Gtf:
    """=============================================================================================

    ============================================================================================="""

    def __init__(self, filename=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.fh = None
        self.filename = filename
        self.line = ''
        self.parsed = {}
        self.exon = 0

        if filename:
            self.opensafe()

    def next(self, skipcomment=True):
        """-----------------------------------------------------------------------------------------
        return the next line from the file

        :param skipcomment: boolean, skip comment lines (begin with #)
        :return: string, next line in file (empty string at eof)
        -----------------------------------------------------------------------------------------"""
        self.line = self.fh.readline()
        if skipcomment:
            while self.line.startswith('#'):
                self.line = self.fh.readline()

        return self.line

    def opensafe(self):
        """-----------------------------------------------------------------------------------------
        Open a file for reading with error trapping.  the filename should be in the BlastResult
        object. Exit status 1 if file open fails

        Usage
            gtf.opensafe()

        :return: filehandle, filehandle of opened file
        -----------------------------------------------------------------------------------------"""
        try:
            self.fh = open(self.filename, 'r')
        except (IOError, OSError):
            sys.stderr.write('Gtf::opensafe - error opening file {}'.format(self.filename))
            exit(1)

        return self.fh

    def parse(self):
        """-----------------------------------------------------------------------------------------
        <seqname>   The sequence ID
        <source>    a unique label indicating where the annotations came from - a prediction program or a public database.
        <feature>   The following feature types are likely to be present
            CDS - the stop codon is not included in the CDS for the terminal exon
            start_codon - first base of start codon
            stop_codon - first base of stop codon
            exon
        <start> <end>   Integer start and end coordinates of the feature <start> <= <end>.
        <score> score for a predicted feature, . means NA
        <strand> + / -
        <frame> 0 indicates that the first whole codon of the reading frame is located at 5'-most base.
        <attributes> text string with gene id, transcipt id, etc.

        :return:
        -----------------------------------------------------------------------------------------"""
        field = self.line.rstrip().split('\t', maxsplit=8)
        self.parsed = {'seqname':   field[0],
                       'source':    field[1],
                       'feature':   field[2],
                       'start':     field[3],
                       'end':       field[4],
                       'score':     field[5],
                       'strand':    field[6],
                       'frame':     field[7],
                       'attribute': field[8],
                       }
        return True

    def attribute2gff(self):
        """-----------------------------------------------------------------------------------------
        Convert
            "jg2008" -> ID=jg2008;                      : naked tag only seen for genes
            gene_id "jg2008" -> Parent=jg2008
            transcript_id "jg7073.t1"; -> ID=jg7073.t1
        :return:
        -----------------------------------------------------------------------------------------"""
        new = ''

        feature = self.parsed['feature']
        attr = self.parsed['attribute'].strip()

        if attr.endswith(';'):
            # remove trailing ;
            # there should't be any , but there are
            attr = attr[:-1]

        if feature == 'gene':
            # only one attribute, the gene ID
            new = 'ID={}'.format(attr)

        else:
            # remove all quotes and split into tag-field pairs
            note = attr.replace('"', '').split(';')

            for n in note:
                n = n.strip()
                tag, field = n.split(' ', maxsplit=1)
                if tag == 'gene_id':
                    gid = field
                else:
                    tid = field

            if feature == 'transcript':
                new = 'ID={};Parent={}'.format(tid, gid)
                self.exon = 0
            elif feature == 'CDS':
                eid = tid + '.e{}'.format(self.exon)
                self.exon += 1
                new = 'ID={};Parent={}'.format(eid, tid)
            elif feature == 'start_codon':
                new = 'ID={}.start;Parent={}'.format(tid, tid)
            elif feature == 'stop_codon':
                new = 'ID={}.stop;Parent={}'.format(tid, tid)

        self.parsed['attribute'] = new.rstrip(';')

        return True


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    gtf = Gtf('C:/Users/michael/Desktop/augustus.hints.gtf')
    gff = open('Hp.augustus.hints.gff', 'w')

    target = ['gene', 'transcript', 'CDS', 'stop_codon', 'start_codon']

    n = 0
    while gtf.next():
        gtf.parse()
        gtf.attribute2gff()
        # print(gtf.line)
        # print(gtf.parsed)

        if gtf.parsed['feature'] in target:
            out = ''
            for k in gtf.parsed:
                if k == 'seqname':
                    underline = gtf.parsed[k].find('_')
                    gtf.parsed[k] = gtf.parsed[k][:underline]
                    scaffold = int(gtf.parsed[k].replace('scaffold',''))
                    if scaffold < 8:
                        continue
                out += gtf.parsed[k] + '\t'
            if scaffold >= 8:
                out = out.rstrip('\t')
                gff.write('{}\n'.format(out))

        # n += 1
        # if n > 25:
        #     break

    exit(0)
