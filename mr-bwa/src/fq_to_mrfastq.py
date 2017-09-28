#!/usr/bin/env python

import sys
import os
from os.path import isfile
from pybwa import MR_FASTQ_LINE_SEPARATOR


# This program takes a fastq file and puts each read in a single line with its fields separated by <sep>.
def main():
    if len(sys.argv) < 2:
        print 'Usage: fq_to_mrfastq.py <input_file.fastq>'
        os.abort()

    input_fastq_filename = str(sys.argv[1])
    if not isfile(input_fastq_filename):
        print '{0} is not a valid input file'.format(input_fastq_filename)
        os.abort()

    input_fastq_file = open(input_fastq_filename, 'r')
    output_fastq_file = open('output.mr.fastq', 'w')

    i = 1
    joined_line = ''

    for line in input_fastq_file:
        # fastq files have 4 lines per read and we want to join these 4 lines separated by <sep>
        # and write them in 1 line. This is crucial since MapReduce streams one line at a time, however in our case,
        # every 4 lines mark a read.
        # Eventually these are split again in mapper and re-assemble the original fastq file
        if i % 4 == 0:
            joined_line += line.replace('\n', '')
            output_fastq_file.write(joined_line + '\n')
            joined_line = ''
        else:
            joined_line += line.replace('\n', MR_FASTQ_LINE_SEPARATOR)

        i += 1

    input_fastq_file.close()
    output_fastq_file.close()


if __name__ == '__main__':
    main()
