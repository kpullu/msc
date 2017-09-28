#!/usr/bin/env python

# import modules
from itertools import groupby
from operator import itemgetter
import sys
import time
from utils import logger, query_bps_count_index

# Reads the mapper output and returns one line at a time
def read_mapper_output(data_input, main_separator='\t'):
    # Go through each line
    for line in data_input:
        # Strip out the separator character
        yield line.rstrip().split(main_separator, 1)


def main(main_separator='\t', tuple_separator=';', list_separator=','):
    start_time = time.time()
    logger.info('Combiner Start Time: {0}'.format(start_time))

    # Read the data using read_mapper_output
    data = read_mapper_output(sys.stdin, main_separator)

    for ref_index, group in groupby(data, itemgetter(0)):
        # group is a list of [ref_index, ref_char;csv_list of combined counts of ACGTDN]
        try:
            total_count_A = total_count_C = total_count_G = total_count_T = total_count_D = total_count_N = 0
            ref_char = None

            for ref_index, ref_char_and_read_counts in group:
                # this represents ref_char;csv_list of combined counts of ACGTDN
                ref_char_name_and_read_counts_tuple = str(ref_char_and_read_counts).split(tuple_separator)

                ref_char = ref_char_name_and_read_counts_tuple[0]
                combined_counts = ref_char_name_and_read_counts_tuple[1].split(list_separator)

                total_count_A += int(combined_counts[query_bps_count_index('A')])
                total_count_C += int(combined_counts[query_bps_count_index('C')])
                total_count_G += int(combined_counts[query_bps_count_index('G')])
                total_count_T += int(combined_counts[query_bps_count_index('T')])
                total_count_D += int(combined_counts[query_bps_count_index('D')])
                total_count_N += int(combined_counts[query_bps_count_index('N')])

            print '%s%s%s%s%s' % (
                ref_index, main_separator, ref_char, tuple_separator,
                ','.join(
                    [str(total_count_A), str(total_count_C), str(total_count_G), str(total_count_T), str(total_count_D),
                     str(total_count_N)]))

        except ValueError as err:
            logger.error('Error: {0}'.format(err))
            pass

    logger.info('Total Combiner Time: {0} seconds'.format(round(time.time() - start_time, 2)))

if __name__ == '__main__':
    main()