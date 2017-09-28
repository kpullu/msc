#!/usr/bin/python
import pickle
import time
import sys
import bwt_fmindex
from utils import query_bps_count_index, logger


# Loads the binary serialized index file
#
# @param index filename
# @returns the binary index
def load_serialized_file(filename):
    f = open(filename)
    idx = pickle.load(f)
    return idx


# Loads reference genome file
#
# @param reference genome filename
# @returns the reference genome text
def read_reference_genome(filename):
    text_file = open(filename, 'r')
    genome = text_file.read()
    text_file.close()
    return genome


# Reads the input data
#
# @param input reads
# @returns one read at a time
def read_input(data_input):
    for line in data_input:
        # Go through each line
        yield line.rstrip()


def main(main_separator='\t', tuple_separator=';'):

    start_time = time.time()
    logger.info('Mapper Start Time: {0}'.format(start_time))

    # load index
    bwt_fm_idx = load_serialized_file('/data/index/hg38_idx')
    if bwt_fmindex is None:
        logger.error('Error while loading reference FM-Index')
        raise RuntimeError('Error while loading reference FM-Index')

    # load ref genome
    ref_gen = read_reference_genome('/data/index/hg38.fa')
    if bwt_fmindex is None:
        logger.error('Error while loading reference genome')
        raise RuntimeError('Error while loading reference genome')

    # load reads
    input_reads = read_input(sys.stdin)

    # this map should contain all distinct matched ref_index together with the count of matched nucleotides
    index_alignments_map = dict()

    for read in input_reads:
        first_occurrence = bwt_fmindex.first_occurrence(read, bwt_fmindex=bwt_fm_idx, mismatches=2)
        if first_occurrence != -1:
            ref_index = first_occurrence
            read_counter = 0
            for i in range(len(read)):
                ref_bp = ref_gen[ref_index]
                query_bp = read[read_counter]

                if ref_index in index_alignments_map:
                    current_tuple = index_alignments_map[ref_index]
                    bps_counts = current_tuple[1]

                    bp_list_index = query_bps_count_index(query_bp)
                    bps_counts[bp_list_index] += 1
                    index_alignments_map[ref_index] = current_tuple[0], bps_counts
                else:
                    # initialize bps counts to 0.
                    # IT'A MUST USING A LIST SINCE ORDER IS IMPORTANT
                    bps_counts = [0] * 6
                    bp_list_index = query_bps_count_index(query_bp)
                    bps_counts[bp_list_index] += 1
                    index_alignments_map[ref_index] = ref_bp, bps_counts

                ref_index += 1
                read_counter += 1

    for ref_index in index_alignments_map.keys():
        # output format ref_index, ref_char;csv_list of combined counts of ACGTDN
        ref_char_name_and_read_counts_tuple = index_alignments_map[ref_index]
        counts_csv = ','.join(str(i) for i in ref_char_name_and_read_counts_tuple[1])

        print '%d%s%s%s%s' % (
            ref_index, main_separator, ref_char_name_and_read_counts_tuple[0], tuple_separator, counts_csv)

    logger.info('Total Mapper Time: {0} seconds'.format(round(time.time() - start_time, 2)))


if __name__ == '__main__':
    main()
