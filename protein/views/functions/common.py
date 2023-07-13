import re

def turn2fasta(seq):
    reg_W = re.compile('\W+')
    reg_fasta = re.compile('^>\w.*\n[\w]{3,}', re.M)
    reg_plainSeq = re.compile('^[A-Za-z]+$', re.M)
    reg_genebank = re.compile('^\d+\s+[A-Z]+\s+', re.I)
    reg_genebankRep = re.compile('\d+|\s+')
    reg_blank = re.compile('\d+|\s+')
    # seq = reg_blank.sub('', seq) + '\n'
    seq = seq.upper()
    print(seq)
    if reg_fasta.match(seq):
        print("seq is fasta")
    elif reg_plainSeq.match(seq):
        print('seq is plain')
        seq = reg_blank.sub('', seq)
        seq = f'>id\n' + reg_blank.sub('', seq) + '\n'

    elif reg_genebank.match(seq):
        print('seq is genebank\n', seq)
        seq = f'>id\n' + reg_genebankRep.sub('', seq) + '\n'
    else:
        return False
    return seq