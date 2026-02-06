"""#################################################################################################
Another GFF/GTF (GXF) parser with the individual records and collections of records split into two
classes. This should make it easier reuse.

Michael Gribskov     21 November 2025
#################################################################################################"""
import sys
version = '2.0.0'


class GxfRecord:
    """=============================================================================================
    The nine standard columns appear as attributes of the GxfRecord object.
    The standard columns are available in the class variable GxfRecord.column

    The GFF standard defines columns 0-7, column 8 is used for any attributes the project defines
    For details see https://useast.ensembl.org/info/website/upload/gff.html or
    https://agat.readthedocs.io/en/latest/gxf.html

    Column 0: "seqid"
      The ID of the landmark used to establish the coordinate system for the current feature,
      typically the ID of the chromosome or scaffold. IDs may contain any characters, but must
      escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
      contain unescaped whitespace and must not begin with an unescaped ">".
    Column 1: "source"
      Free text describing the origin of the feature annotation. Typically, this is the name of
      program, such as "Genescan" or a database name, such as "Genbank."
    Column 2: "type"
      The feature type. This is constrained to be either a term from the Sequence Ontology. Common
      features include gene, transcript, mRNA, 5'UTR, 3'UTR, exon, CDS, start_codon, stop_codon, etc
    Columns 3 & 4: "start" and "end"
      The start and end coordinates of the feature in 1-based integer coordinates, relative to the
      sequence in column 1. Start is always less than or equal to end. For features that cross the
      origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes),
      end end + the length of sequence, i.e., the end will be past the end of the sequence.
      For zero-length features, such as insertion sites, start equals end and the implied site is to
      the right of the indicated base in the direction of the landmark.
    Column 5: "score"
      The score of the feature, a floating point number. As in earlier versions of the format,
    Column 6: "strand"
      The strand of the feature. + for positive strand , - for minus strand, and . for features that
      are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.
    Column 7: "phase"
      For features of type "CDS", the phase indicates where the next codon begins relative to the 5'
      end of the current CDS feature. The 5' end for CDS features on the plus strand is the feature
      start and the 5' end for CDS features on the minus strand is the feature's end. The phase is
      0, 1, or 2, indicating the number of bases forward from the start of the current CDS feature
      the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide of
      the CDS feature (i.e. 0 bases forward). Note that ‘Phase’ in the context of a GFF3 CDS feature
      should not be confused with the similar concept of frame that is also a common concept in
      bioinformatics. Frame is generally calculated as a value for a given base relative to the
      start of the complete open reading frame (ORF).

      The phase is REQUIRED for all CDS features.
    Column 8: GTF/GFF2
     A semicolon-separated list of tag-value pairs, providing additional information about each
     feature. Except for the terminal semicolon, semicolons should be followed by spaces

     Each attribute is a key:value pair. Keys are expected to be single tokens; Values must be
     enclosed in parentheses so that spaces may be present. Attributes may contain leading/trailing
     spaces. gene_id and transcript_id are required, but may be empty strings.

    Column 8: GFF3
      A semicolon-separated list of feature attributes in the format tag=value. Multiple tag=value
      pairs are separated by semicolons. URL escaping rules are used for tags or values containing
      the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced
      with the %09 URL escape. Attribute values do not need to be and should not be quoted unless
      they are actually part of the value. Quotes should be included as part of the value by parsers
      and not stripped.

      For attributes with predefined meanings see
      https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    ============================================================================================="""
    # Predefined keys for columns 0-8, agrees with GFF3
    column = {'seqid': '', 'source': '', 'type': '', 'start': '0', 'end': '0',
              'score': '', 'strand': '', 'phase': '', 'attribute': None}

    def __init__(self, row, fmt="gtf"):
        """-----------------------------------------------------------------------------------------
        format is not a class variable because you may want to have both gff3 and gtf formats
        active at the same time

        :param row: string  a string that can be parsed as gff3 or gtf
        :param fmt: string  GTF is the default, use gff3 for GFF3
        -----------------------------------------------------------------------------------------"""
        for key, value in GxfRecord.column.items():
            # set up the standard columns as attributes
            setattr(self, key, value)

        # attr_sep separates the key/value pairs in attribute column
        self.format = None
        self.attr_sep = None
        self.format_set(fmt)

        if row:
            self.feature_parse(row)

    def format_set(self, fmt):
        """-----------------------------------------------------------------------------------------
        Set up GFF3 or GTF format. GFF2 is the same as GTF.

        :param fmt: str      format for attributes gff3|gtf|gff2
        :return: str         format
        -----------------------------------------------------------------------------------------"""
        if fmt == 'gff3' or fmt == 'gff':
            self.format = 'gff3'
            self.attr_sep = '='
        else:
            self.format = 'gtf'
            self.attr_sep = ' '  # space

        return self.format

    def feature_parse(self, row):
        """-----------------------------------------------------------------------------------------
        parse a feature line. Attributes are split on semicolons, and key/value pairs split on
        attr_sep character.

        :param row: str     string to parse
        :return: True
        -----------------------------------------------------------------------------------------"""
        field = row.split(maxsplit=8)
        col_n = 0
        for key in GxfRecord.column:
            # extract the 9 defined columns
            if col_n in [3, 4]:
                # change begin, end to int
                # should score be float?
                setattr(self, key, int(field[col_n]))

            else:
                try:
                    setattr(self, key, field[col_n])
                except:
                    print(f'GxfRecord: error in feature_parse(); line {sys._getframe().f_lineno}')

            # split the attributes on ; and save as a hash in self.attributes
            attrs = {}
            if key == 'attribute':
                field = field[col_n].rstrip().split(';')
                # attribute may end in; so last field may be blank
                if not field[-1]:
                    field.pop()

                for f in field:
                    try:
                        (key, value) = f.strip().split(self.attr_sep, maxsplit=1)
                    except:
                        print(f'GxfRecord: error in feature_parse(); line {sys._getframe().f_lineno}')
                    # TODO check this works for GFF3
                    attrs[key] = value.strip('"')

                setattr(self, 'attribute', attrs)

            col_n += 1

        return True

    @staticmethod
    def attribute_format(attr_dict, fmt):
        """-----------------------------------------------------------------------------------------
        Return a string containing the attributes (column 8) formatted as GTF or GFF3 according to

        :param attr_dict: dict      attributes to be formatted
        :param fmt: str             gtf|gff3
        :return: str                formatted string ready to write in column 8
        -----------------------------------------------------------------------------------------"""
        attr_str = ''
        if fmt == 'gff3':
            for attr in sorted(attr_dict, key=lambda x: (x != 'transcript_id', x != 'gene_id', x)):
                attr_str += f'{attr}={attr_dict[attr]}; '

        elif fmt == 'gtf' or fmt == 'gff2':
            for attr in sorted(attr_dict, key=lambda x: (x != 'transcript_id', x != 'gene_id', x)):
                attr_str += f'{attr} "{attr_dict[attr]}"; '

        else:
            print(f'GxfRecord:attribute_format - unknown format ({fmt})')

        return attr_str


class GxfSet:
    """=============================================================================================
    A container for GxfRecords with tools for selecting sets of features from files or GxfSets

    ============================================================================================="""

    def __init__(self, file='', fmt='gtf'):
        """-----------------------------------------------------------------------------------------
        file        string          path to data file
        fh          filehandle      open file (read or write)
        features    list            list of GxfRecord
        -----------------------------------------------------------------------------------------"""
        self.fh = None
        self.features = []
        self.format = fmt

        if file:
            self.fh = self.opensafe(file)
        else:
            self.file = 'in.gtf'

    def opensafe(self, file, mode='r'):
        """-----------------------------------------------------------------------------------------
        Open file with check for whether the file exists. If the file cannot be opened, exit with
        status=1

        :param file: string     path to file
        :param mode: string     r|w (or any valid mode)
        :return: filehandle     opened file, also stored in self.fh
        -----------------------------------------------------------------------------------------"""
        if file:
            self.file = file

        try:
            fh = open(file, mode)
        except OSError:
            print(f'GxfSet::opensafe - Error opening file ({file})')
            exit(1)

        self.fh = fh
        return fh

    def feature_get(self, choice):
        """-----------------------------------------------------------------------------------------
        Read all matching features from self.fh. choice is a list of desired features

        :param choice: list     features to select
        :return: int            number of features read
        -----------------------------------------------------------------------------------------"""
        fh = self.fh
        for line in fh:
            line = line.rstrip()
            if not line:
                # blank lines (really shouldn't be there, but sometimes they are)
                continue
            if line.startswith('#') or line.startswith(chr(0)):
                # skip comments
                continue

            # print(line)
            record = GxfRecord(line, fmt=self.format)
            if record.type in choice:
                self.features.append(record)

        return len(self.features)


# ##################################################################################################
# Testing
# #################################################################################################
if __name__ == '__main__':
    gtfin = 'data/stringtie_merged_jm.gtf'
    gtf = GxfSet(file=gtfin)
    feature_n = gtf.feature_get(['transcript'])
    print(f'{feature_n} features read from {gtfin}')


    exit(0)
