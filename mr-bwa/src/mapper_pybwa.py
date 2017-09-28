#!/usr/bin/env python
import uuid
from pybwa import logger, query_bps_count_index, MR_FASTQ_LINE_SEPARATOR

import pybwa.bwa
import pysam
import sys
import time


# The input fastq is adapted to MapReduce whereby every 4 lines are joined into a single line, since MapReduce
# streams 1 line at a time and in fastq every 4 lines refer to a single read.
# Here we separate each line into 4 entries again
#
# @param the pre-processed input file joined by '<sep>'
# @returns the input in valid fastq format
def read_input(data_input):
    reads_data = ''

    for line in data_input:
        whole_fastq_entry = line.split(MR_FASTQ_LINE_SEPARATOR)
        reads_data += '\n'.join(sep_line for sep_line in whole_fastq_entry)

    return reads_data


# Write to temp file to be used by bwa as input
#
# @param the valid fastq format data
# @param tempfile to write to
def write_to_file(reads_data, filename):
    text_file = open(filename, 'w')
    text_file.write(reads_data)
    text_file.close()


def main(main_separator='\t', tuple_separator=';'):
    random_id = uuid.uuid4()

    start_time = time.time()
    logger.info('Mapper Start Time: {0}'.format(start_time))

    reads_file = '/data/input_reads_{0}.fq'.format(random_id)
    reference_file = '/data/index/hg38.fa'
    output_file = '/data/aln_input_reads_mem_{0}.sam'.format(random_id)

    # Ensure reference is indexed
    if not pybwa.bwa.index_ref(reference_file):
        logger.error('Error while checking reference index')
        raise RuntimeError('Error while checking reference index')

    # load reads from stdin
    reads_data = read_input(sys.stdin)
    write_to_file(reads_data, reads_file)

    # Setup and run bwa mem
    mem = pybwa.bwa.BWAMem(reference_file, reads_file)
    return_code = mem.run(output_file)

    # Check return status
    if return_code != 0:
        logger.error('Error running bwa')
        raise RuntimeError('General error when running bwa')
    elif return_code == 0:
        # parse sam file
        sam_alignment_file = pysam.AlignmentFile(output_file, 'r')
        sam_iterator = sam_alignment_file.fetch()

        # mapper output <ref_index, <ref_char; ref_chromosone_name; [combined counts of ACGTDN]>> e.g. <531, <G;chr17;0,2,1,0,0,0>>
        # meaning at ref index 531, there is bp 'G' and alignment returned  [0A, 2C, 1G, 0T, 0D, 0N]
        # where D stands for deletion in person's genome and N for no-call in query read

        # this map should contain all distinct matched ref_index together with the count of matched nucleotides
        index_alignments_map = dict()

        # each item in iter is an instance of http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.
        # hence we can get details for each read
        for aligned_segment in sam_iterator:
            if not aligned_segment.is_unmapped:
                # check if aligned segment has insertions, if yes get indices
                insertions_query_indexes = []
                if 'I' in aligned_segment.cigarstring.upper():
                    counter = 0
                    for cigartuple in aligned_segment.cigartuples:
                        # check insertions tuples. 1 means insert
                        if cigartuple[0] == 1:
                            start_index = None
                            if counter == 0:
                                start_index = 0
                            elif aligned_segment.cigartuples[counter - 1][0] == 0:
                                # get start_index from previous cigar tuple, if previous tuple was a match/mismatch
                                # i.e. = 0
                                start_index = aligned_segment.cigartuples[counter - 1][1]

                            if start_index is not None:
                                for i in range(start_index, start_index + cigartuple[1]):
                                    insertions_query_indexes.append(i)

                        counter += 1

                # each aligned pair is a tuple <query_index, reference_index, aligned_ref_bp>
                # A None ref_index refers to insertions in person's genome
                # A None query_index refers to deletions in person's genome
                tuple_counter = 0
                for aligned_tuple in aligned_segment.get_aligned_pairs(with_seq=True):
                    query_index = aligned_tuple[0]
                    ref_index = aligned_tuple[1]
                    aligned_ref_bp = aligned_tuple[2]

                    # get aligned ref sequence
                    # aligned_segment.get_reference_sequence()

                    # get aligned query bp
                    if query_index is None:
                        # we have a deletion in person's genome, hence we mark it with 'D'
                        query_bp = 'D'
                    else:
                        query_bp = aligned_segment.query_sequence[query_index]

                    if ref_index is None and query_index in insertions_query_indexes:
                        # handle insertions by updating ref_index to 'ref_index.n', hence we know that at nth index after ref_index, there are insertions
                        # set insert_ref_index to last ref_index before insertions.
                        insert_ref_index = None
                        if query_index - 1 not in insertions_query_indexes and \
                                        aligned_segment.get_aligned_pairs()[tuple_counter - 1][1] is not None:
                            insert_ref_index = aligned_segment.get_aligned_pairs()[tuple_counter - 1][1]

                        # append '.n' where n is the index in insertions_query_indexes +1.
                        # we use . to ensure indexes are sorted numerically
                        # cater up to 99 consecutive insertions
                        appended_index = insertions_query_indexes.index(query_index) + 1;
                        appended_index_str = str(appended_index)
                        if appended_index < 10:
                            appended_index_str = '0' + appended_index_str
                        if insert_ref_index is not None:
                            ref_index = float(str(insert_ref_index) + '.' + appended_index_str)

                    if ref_index in index_alignments_map:
                        current_tuple = index_alignments_map[ref_index]
                        bps_counts = current_tuple[2]

                        bp_list_index = query_bps_count_index(query_bp)
                        bps_counts[bp_list_index] += 1
                        index_alignments_map[ref_index] = current_tuple[0], current_tuple[1], bps_counts

                    elif ref_index is not None:
                        # initialize bps counts to 0.
                        # IT'A MUST USING A LIST SINCE ORDER IS IMPORTANT
                        bps_counts = [0] * 6
                        bp_list_index = query_bps_count_index(query_bp)
                        bps_counts[bp_list_index] += 1
                        index_alignments_map[
                            ref_index] = None if aligned_ref_bp is None else aligned_ref_bp.upper(), aligned_segment.reference_name, bps_counts

                    tuple_counter += 1

        for ref_index in index_alignments_map.keys():
            # output format ref_index, ref_char;ref_chromosone_name;csv_list of combined counts of ACGTDN
            ref_char_name_and_read_counts_tuple = index_alignments_map[ref_index]
            counts_csv = ','.join(str(i) for i in ref_char_name_and_read_counts_tuple[2])
            if float(ref_index).is_integer():
                print '%d%s%s%s%s%s%s' % (
                    ref_index, main_separator, ref_char_name_and_read_counts_tuple[0], tuple_separator,
                    ref_char_name_and_read_counts_tuple[1], tuple_separator, counts_csv)
            else:
                print '%.2f%s%s%s%s%s%s' % (
                    ref_index, main_separator, ref_char_name_and_read_counts_tuple[0], tuple_separator,
                    ref_char_name_and_read_counts_tuple[1], tuple_separator, counts_csv)

        logger.info('Total Mapper Time: {0} seconds'.format(round(time.time() - start_time, 2)))


if __name__ == '__main__':
    main()