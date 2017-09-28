#!/usr/bin/python
import os
from os.path import isfile
import cPickle as pickle
import sys
import time
import bwt_fmindex


# Serializes the input index to binary using cPickle and saves
#
# @param filename to wirte to
# @param index data
def save(filename, idx):
    f = open(filename, 'w')
    # files created with protocols >= 1 are in binary mode (https://docs.python.org/2/library/pickle.html)
    pickle.dump(idx, f, pickle.HIGHEST_PROTOCOL)
    f.close()


def main():
    if not len(sys.argv) in [3]:
        print 'Usage: %s input_file output_index_file' % sys.argv[0]
        os.abort()
    else:
        if not isfile(sys.argv[1]):
            print 'Input file does not exist'
            os.abort()

        inp = open(sys.argv[1])
        # read input
        reference_data = inp.read()
        inp.close()

        start_time = time.time()
        print 'Started indexing at {0}'.format(start_time)

        # create index
        fm_idx = bwt_fmindex.make_index(reference_data)

        end_time = time.time()
        print 'Finished indexing at {0}'.format(end_time)
        print 'Total Indexing Time: {0} seconds'.format(round(end_time - start_time, 2))

        # save index to file
        save(sys.argv[2], fm_idx)

if __name__ == '__main__':
    main()
