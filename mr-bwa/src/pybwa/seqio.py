from Bio import SeqIO

# Infers filetype based on the first character in file
#
# @param filename
# @returns filetype
def seqfile_type(filename):
    ftype = 'fasta'
    with open(filename) as fh:
        firstline = fh.readline()
        if firstline.startswith('>'):
            ftype = 'fasta'
        elif firstline.startswith('@'):
            ftype = 'fastq'
        elif firstline.startswith('.sff'):
            ftype = 'sff'
        else:
            raise ValueError('{0} not a valid bwa file'.format(filename))
    return ftype

# Checks the number of reads in file
#
# @param filename
# @returns count of reads
def reads_in_file(filename):
    ftype = seqfile_type(filename)
    return sum([1 for seq in SeqIO.parse(filename, ftype)])
