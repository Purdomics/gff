from gff import Gff

class Bundle:

    def __init__(self):
        """

        """
        self.begin = 0
        self.end = 0
        self.sequence = ''
        self.transcript = []

def make_bundle(gff,id_column='transcript_id'):
    """---------------------------------------------------------------------------------------------

    :param gff: Gff object
    :param id_column: name of id column in gff hash
    :return: list of bundles
    ---------------------------------------------------------------------------------------------"""
    bundle_list = []
    b_stop = 0
    sequence_old = ''
    for t in sorted(gff.data, key=lambda k: (k['sequence'], int(k['begin']))):
        id = t[id_column]
        begin = int(t['begin'])
        end = int(t['end'])
        sequence = t['sequence']
        strand = t['strand']

        if begin > b_stop or sequence != sequence_old:
            # new bundle
            b = Bundle()
            b.begin = begin
            b.end = max(b.end,end)
            b.sequence = sequence
            b.transcript.append(t)
            bundle_list.append(b)
            sequence_old = sequence

        else:
            # continue current bundle
            b.stop = end
            b.transcript.append(t)

        b_stop = end

    return bundle_list

def overlap2(t1, t2):
    """-----------------------------------------------------------------------------------------
    left = t2['begin'] -t1['begin']
    right = t2['end'] - t1['end']

    :param t1:
    :param t2:
    :return:
    -----------------------------------------------------------------------------------------"""
    left = int(t2['begin']) - int(t1['begin'])
    right = int(t2['end']) - int(t1['end'])
    x1 = int(t2['end']) - int(t1['begin'])
    x2 = int(t2['begin']) - int(t1['end'])

    test = ((left<0)<<3) + ((right<=0)<<2) + ((x1<0)<<1) + (x2<0)
    if test == 15:
        status = 'None'
    elif test == 13:
        status = 'merge-right'
    elif test == 9:
        status = 'merge-include'
    elif test == 5:
        status = 'merge-extend'
    elif test == 1:
        status = 'merge-left'
    elif test == 0:
        status = 'None'
    else:
        print('Error unknown overlap status ({}'.format(test))

    return left, right, status

def overlap(t1, t2):
    """-----------------------------------------------------------------------------------------
    Check if t2 overlaps t1.  if t2 is exactly the same as t1, the status will be 'merge-t2inside".
    Retur:
        left_extension: t1['begin'] - t2['begin']
        right_extension: t1['end'] - t2['end']
        status: no_overlpa
                merge_t1inside
                merge_t2inside
                merge_extleft
                merge_extright

    :param t1:
    :param t2:
    :return:int, int, str, left extension, right extension, type
    -----------------------------------------------------------------------------------------"""
    lextend = 0
    rextend = 0
    type = 'no_overlap'

    left = t1
    right = t2
    if left['begin'] > right['begin']:
        # make sure the left wequence begins with the lower coordinate
        left, right = right, left

    if right['begin'] < left['end']:
        # there is an overlap
        if right['end'] <= left['end']:
            # right range is entirely inside left
            if left is t1:
                type = 'merge_t2inside'
            else:
                type = 'merge_t1inside'
        else:
            # right range extends to the right
            if left is t1:
                type = 'merge_extright'
            else:
                type = 'merge_extleft'

        # extensions are given with respect to t1
        lextend = left['begin'] - t1['begin']
        rextend = max(left['end'], right['end']) - t1['end']

    return lextend, rextend, type

def merge(merge, t):
    """
    
    :param merge: 
    :param t: 
    :return: 
    """
    print('C {}\t{}\t{}'.format(t['transcript_id'], t['begin'], t['end']))
    if t['transcript_id'] not in merge['names']:
        if merge['strand'] != t['strand']:
            print('warning different strands {}'.format(t['transcript_id']))
        merge['begin'] = min(merge['begin'], t['begin'])
        merge['end'] = min(merge['end'], t['end'])
        merge['member'].append(t)
        merge['names'].append(t['transcript_id'])
    
if __name__ == '__main__':

    sgff = Gff(file="stranded.merged.stringtie.gff")
    transcript_n = sgff.read_feature('transcript')
    print('{} stranded features read'.format(transcript_n))
    sgff.replace_columns_re(['sequence'], 'lcl\|', r'')
    # sgff.attr_sep = '='
    # transcript_n = sgff.read_feature('mRNA')
    query = r'(maker|augustus|masked|processed|gene|trnascan)(-|_)'
    sgff.replace_columns_re(['gene_id', 'transcript_id'], query, r'')
    sgff.replace_columns_re(['gene_id', 'transcript_id'], r'-tRNA-', r'.tRNA')
    sgff.position_to_int()
    # sgff.rename_key('ID', 'transcript_id')

    sbundle = make_bundle(sgff)

    ugff = Gff(file="unstranded.merged.stringtie.gff")
    transcript_n = ugff.read_feature('transcript')
    print('{} unstranded features read'.format(transcript_n))
    ugff.replace_columns_re(['sequence'], 'lcl\|', r'')
    query = r'(maker|augustus|masked|processed|gene|trnascan)(-|_)'
    ugff.replace_columns_re(['gene_id', 'transcript_id'], query, r'')
    ugff.replace_columns_re(['gene_id', 'transcript_id'], r'-tRNA-', r'.tRNA')
    ugff.position_to_int()
    # ugff.rename_key('ID', 'transcript_id')

    ubundle = make_bundle(ugff)

    set = []
    set_end = 0
    set_seq = ''
    for b in sorted(ubundle+sbundle, key=lambda k: (k.sequence, k.begin)):

        if b.begin > set_end or b.sequence != set_seq:
            # new set
            set_seq = b.sequence
            set_end = b.end
            set.append( {'begin'   : b.begin,
                        'end'     : set_end,
                        'sequence': set_seq,
                        'member'  : [b]} )

        else:
        # continue current set
            set[-1]['end'] = max(set[-1]['end'], b.end)
            set[-1]['member'].append(b)
            set_end = set[-1]['end']

    s_n = 0

    for s in set:
        print('\nset {}\t{}'.format(s_n, s['sequence']))
        s_n += 1
        trans = []
        for b in s['member']:
            # get all members of the set in one list and sort
            for transcript in b.transcript:
                trans.append(transcript)

        trans.sort(key=lambda k: (k['transcript_id'][0]=='C',k['transcript_id'][0]=='S'))

        # if the C (genome) and S (stranded RNA) agree on the strand, the unstranded reads need to
        # be merged on to the known strand

        # The first sorted transcript is the founder of the merge
        names = []
        dir = ''
        t = trans.pop()
        names.append(t['transcript_id'])
        if t['transcript_id'][0] in ('C', 'S'):
            # only set strand for C and S
            dir += t['strand']

        while trans:
            t = trans.pop()
            if t['transcript_id'] in names:
                continue

            if t['transcript_id'][0] in ('C', 'S') and t['strand'] not in dir:
                # only set strand for C and S
                dir += t['strand']
            elif t['strand'] not in dir:
                continue

            names.append(t['transcript_id'])


        print( 'mergetype: {}'.format(dir))
        for name in names:
            print('\t{}'.format(name))


        if s['sequence'] == 'Ctg0003':
            break
        # break

    exit(0)
