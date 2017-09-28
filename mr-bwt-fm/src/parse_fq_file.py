#!/usr/bin/env python

import os
import sys
from os.path import isfile


# Extracts reads from fastq files by returning 2nd row of every 4 lines
#
# @param input fastq file
# @returns the reads from the fastq file
def extractReadsFromFASTQ(filename):
    f = open(filename)

    initial_line = 1
    counter = 0

    for line in f:
        if len(line) > 0:
            # fastq files have 4 rows per read. read is found on 2nd row of every 4 lines
            if counter == initial_line or counter == initial_line + 4:
                initial_line = counter
                yield line

            counter += 1

    f.close()


def main():
    if len(sys.argv) < 2:
        print 'Usage: parse_fq_file.py <input_file.fastq>'
        os.abort()

    input_fastq_filename = str(sys.argv[1])

    if not isfile(input_fastq_filename):
        print '{0} is not a valid input file'.format(input_fastq_filename)
        os.abort()

    reads = extractReadsFromFASTQ(input_fastq_filename)
    output_reads_file = open('output.fq.reads', 'w')
    for r in reads:
        output_reads_file.write(r)

    output_reads_file.close()


if __name__ == '__main__':
    main()
