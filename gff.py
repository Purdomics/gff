"""#################################################################################################
gff.py

contains Gff class
#################################################################################################"""
import sys
import re


class Dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Gff:
    """=============================================================================================
    Class for manipulating both GFF and GTF files
    Note that the contents of data are not explicitly defined so they can be extended without
    subclassing

     9 October 2019    Michael Gribskov
    ============================================================================================="""

    column = ['sequence', 'method', 'feature', 'begin', 'end', 'score', 'strand', 'frame',
              'attribute']

    def __init__(self, file="", mode="GFF"):
        """-----------------------------------------------------------------------------------------
        Data is initially stored as a list of dicts with the keys in the class variable column.
        After reading the attribute column is split into key-value pairs and the keys used to
        store each attribute separately in the dict.

        data: list of dicts
        gff_in: filehandle for gff file
        attr_sep: the inner separator for attributes. key-value pairs are separated by ; and the
            keys and values separated by attr_sep

        :param file: string     path to file with GFF/GTF data
        :param mode: string     should be 'GFF' or 'GTF'
        -----------------------------------------------------------------------------------------"""
        self.data = []
        self.gff_in = None
        self.mode = 'GFF'
        self.attr_sep = '='
        self.setmode(mode)

        if file:
            self.open(file)

    def open(self, file):
        """-----------------------------------------------------------------------------------------
        safely open a file for reading

        :param self: gff
        :param file: str, path to a GFF file
        :return: fh or False if unsuccessful
        -----------------------------------------------------------------------------------------"""
        try:
            fh = open(file, 'r')

        except IOError:
            sys.stderr.write("gff.open unable to open GFF file ({})".format(file))
            exit(1)

        self.gff_in = fh
        return fh

    def setmode(self, mode):
        """-----------------------------------------------------------------------------------------
        changing from GFF to GTF mode requires changing the attribute separator self.attr_set
        :param mode: string     'GFF' or 'GTF'
        :return: string         new value for self.mode
        -----------------------------------------------------------------------------------------"""
        if self.mode == mode:
            # correct mode is already set
            return
        elif mode == 'GFF':
            self.mode == 'GFF'
            self.attr_sep = '='
        elif mode == 'GTF':
            self.mode == 'GTF'
            self.attr_sep = ' '
        else:
            sys.stderr.write(f'gff.setmode - unknown mode ({mode}), mode is {self.mode}')

        return self.mode

    def attribute_add(self, attr, value, begin=0, end=0):
        """-----------------------------------------------------------------------------------------
        Add an attribute to the data list with name attr and value value. if begin and end are not
        provided, the attribute is added to every row

        :param attr: string     key for new attribute
        :param value: any       initial value for new attribute
        :param begin: int       first row
        :param end: int         last row to modify + 1
        :return: int            number of rows modified
        -----------------------------------------------------------------------------------------"""
        if end == 0:
            end = len(self.data)

        # check that attr is unique
        if attr in self.data[begin]:
            sys.stderr.write(f'Gff.attribute_add - attribute ({attr}) already exists in Gff.data')
            return 0

        for row in self.data[begin:end]:
            row[attr] = value

        return end - begin

    def read(self):
        """-----------------------------------------------------------------------------------------
        Read a line from gff_in and delegate to the proper parsing function

        :return:
        -----------------------------------------------------------------------------------------"""
        line = self.gff_in.readline().rstrip()

        if line:
            if line.startswith('#'):
                self.comment_parse()
            else:
                self.data.append(self.feature_parse(line))

            return True

        else:
            # EOF
            return False

    def read_all(self):
        """-----------------------------------------------------------------------------------------
        read the entire file, all features, into self.data. lines are parsed.

        :return: int, number of lines read
        -----------------------------------------------------------------------------------------"""
        nline = 0
        while self.read():
            nline += 1

        return nline

    def read_feature(self, feature_list):
        """-----------------------------------------------------------------------------------------
        Read the whole file and store only features in feature_list in self.data

        :param feature_list:
        :return: int, number of features read
        -----------------------------------------------------------------------------------------"""
        count = 0
        for line in self.gff_in:
            if line:
                if line.startswith('#'):
                    self.comment_parse()

                else:
                    parsed = self.feature_parse(line)
                    if parsed['feature'] in feature_list:
                        self.data.append(parsed)
                        count += 1

        return count

    def feature_parse(self, line):
        """-----------------------------------------------------------------------------------------
        parse a feature line. the final field holds attributes in key value format. for GFF,
        the format is
        Parent=maker-Ctg0001-augustus-gene-0.4;ID=maker-Ctg0001-augustus-gene-0.4.mRNA1;Name=maker-Ctg0001-augustus-gene-0.4.mRNA1;Alias=maker-Ctg0001-augustus-gene-0.4,maker-Ctg0001-augustus-gene-0.4-mRNA-1;mRNA=maker-Ctg0001-augustus-gene-0.4.mRNA1;coge_fid=936743213
        for GTF is is
        gene_id "MSTRG.13"; transcript_id "MSTRG.13.3"; exon_number "1";

        :return:
        -----------------------------------------------------------------------------------------"""
        field = line.rstrip().split(maxsplit=8)
        parsed = {}
        for i in range(len(field)):
            # extract the 9 defined columns
            parsed[Gff.column[i]] = field[i]

        # set numerical fields to int
        for col in ('begin', 'end', 'frame'):
            if parsed[col] != '.':
                parsed[col] = int(parsed[col])
        # parsed['score'] = float(parsed['score']) # score is usually '.'

        # split the attributes on ; and restore as a hash

        field = parsed['attribute'].rstrip().split(';')
        # attribute may end in; so last field may be blank
        if not field[-1]:
            field.pop()

        for f in field:
            (key, value) = f.strip().split(self.attr_sep, maxsplit=1)
            parsed[key] = value.replace('"', '')

        # self.data.append(parsed)

        return parsed

    def comment_parse(self):
        """-----------------------------------------------------------------------------------------
        Parse comments, not implemented

        :return: True
        -----------------------------------------------------------------------------------------"""
        pass
        return True

    def get_by_feature(self, key, start=0):
        """-----------------------------------------------------------------------------------------
        Generator for lines that match a feature tag

        :param key: str         string matching a feature in column 3
        :param start: int       row on which to start
        :yield: int, Dotdict    row number, next matching entry
        -----------------------------------------------------------------------------------------"""
        data = self.data
        for n in range(start, len(data)):
            if data[n]['feature'] == key:
                yield n, data[n]

        return

    def get_by_sequence(self, seqid):
        """-----------------------------------------------------------------------------------------
        Generator for entries that match a specific sequence (column 0)

        :param seqid: string        id of sequence
        :yield: Dotdict             next matching entry
        -----------------------------------------------------------------------------------------"""
        data = self.data
        for row in data:
            if row['sequence'] == seqid:
                yield Dotdict(row)

        return

    def get_by_value(self, column, key, start=0, stop=0):
        """-----------------------------------------------------------------------------------------
        A generator that returns rows where the specified column matches the specified value.
        Rows missing columns, e.g., those generated from attributes, are skipped

        :param column: str, predefined or attribute column
        :param key: str, column value to match
        :param start: int, beginning row
        :param stop: int, ending row + 1
        :return: row, content_hash
        -----------------------------------------------------------------------------------------"""
        data = self.data
        if stop == 0:
            stop = len(data)

        for n in range(start, stop):
            if column not in data[n]:
                continue

            if data[n][column] == key:
                yield n, data[n]

        return

    def replace_by_column(self, column, find, replace):
        """-----------------------------------------------------------------------------------------
        Replace all of find by replace in column

        :param column: str, name of a predefined or attribute column
        :param find: str, string to replace
        :param replace: string to substitute for find
        :return: int, rows examined
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0
        for d in data:
            if column in d:
                d[column] = d[column].replace(find, replace)
                n += 1

        return n

    def replace_columns_re(self, column_list, search, replace=''):
        """-----------------------------------------------------------------------------------------
        replace strings in a list of columns using a regular expression
        -----------------------------------------------------------------------------------------"""
        data = self.data
        n = 0

        query = re.compile(search)
        for d in data:
            for column in column_list:
                if column in d:
                    d[column] = query.sub(replace, d[column])
                    n += 1

        return n

    def rename_key(self, old_key, new_key):
        """-----------------------------------------------------------------------------------------
        changes old_key to new_key in data hash.  Does nothing if old_key does not exist.  Useful
        for changing names of attribute columns when comparing GFF files from different sources.

        :param old_key: str, current key
        :param new_key: str, new _key
        :return: True if old_key exists
        -----------------------------------------------------------------------------------------"""
        data = self.data
        if old_key in data[0]:
            for d in data:
                d[new_key] = d.pop(old_key)

            return True

        return False

    def position_to_int(self):
        """-----------------------------------------------------------------------------------------
        convert begin and end positions from str to int
        :return: 
        -----------------------------------------------------------------------------------------"""
        for d in self.data:
            d['begin'] = int(d['begin'])
            d['end'] = int(d['end'])

        return True


# ==================================================================================================
# test
# ==================================================================================================
if __name__ == '__main__':
    genome = None
    test_feature = ''
    test_value = ['', '']
    test = {'gff':          False,
            'gtf':          True,
            'read_all':     True,
            'read_feature': False,
            'get_by_value': True,
            'examples':     True}

    if test['gff']:
        # read gff, example is from CoGe comparative genomics
        genome = Gff(file='data/genome.gff')
        genome.attr_sep = '='
        test_feature = ['mRNA', 'exon']
        test_value = ['sequence', 'Ctg0001']

    if test['gtf']:
        # read gtf file, example is from stringtie, v2.0.3
        genome = Gff(file='data/stringtie.gtf')
        genome.attr_sep = ' '
        test_feature = ['transcript']
        test_value = ['gene_id', 'MSTRG.1']

    if test['read_all']:
        genome.read_all()

    if test['read_feature']:
        nline = transcripts = genome.read_feature(test_feature)
        sys.stdout.write(f'{nline} lines read\n')

    if test['get_by_value']:
        # get_by_value is a generator function
        geneid = 'ID'
        if test['gtf']:
            geneid = 'transcript_id'
        for n, data in genome.get_by_value(test_value[0], test_value[1]):
            print(f"sequence:{data['sequence']}\t{data['begin']}\t{data['end']}\t\t{data[geneid]}")

    if test['examples']:
        # some examples of other functions
        # remove the string 'lcl|' in the sequence names
        target = 'lcl|'
        genome.replace_by_column('sequence', target, '')

        # with regular expression
        query = r'(maker|augustus|masked|processed|gene)(-|_)'
        genome.replace_columns_re(['Parent', 'ID'], query, r'')

        # transcripts is a generator function
        feature = 'mrna'
        geneid = 'ID'
        if test['gtf']:
            feature = 'transcript'
            geneid = 'transcript_id'

        for n, data in genome.get_by_feature(feature):
            print(f"sequence:{data['sequence']}\t{data['begin']}\t{data['end']}\t\t{data[geneid]}")

    exit(0)
