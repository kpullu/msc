# This script is an adaptation from a library used to interact with BWA mapper using Python
# https://github.com/VDBWRAIR/pyBWA

import re
from subprocess import Popen, PIPE
import os
import os.path
import glob
import tempfile
import sh
from pybwa import logger
import seqio


# Checks to see if a given reference is indexed already
#
# @param - Reference file name
# @return True if ref is indexed, False if not
def is_indexed(ref):
    ref_ext = set(['.amb', '.ann', '.bwt', '.pac', '.sa'])
    ref_indexes = glob.glob(ref + '.*')
    logger.debug('Indexes found for {0}: {1}'.format(ref, ref_indexes))

    if len(ref_indexes) == 0:
        return False

    ext = set([os.path.splitext(f)[1] for f in ref_indexes])
    intersec = ref_ext & ext

    # Return true only if the intersection of the found extensions is equal to
    # the expected ref_ext set
    return intersec == ref_ext


# Indexes a given reference
#
# @param ref - Reference file path to index
def index_ref(ref):
    # Don't reindex an already indexed ref
    if is_indexed(ref):
        logger.debug('{0} is already indexed'.format(ref))
        return True

    if not os.path.exists(ref):
        logger.critical('Reference path {0} cannot be read'.format(ref))
        return False

    logger.info('Indexing {0}'.format(ref))

    try:
        ret = BWAIndex(ref).run()
    except ValueError as e:
        logger.error(e)

    if ret != 0:
        logger.error('Error running bwa index on {0}'.format(ref))
        return False
    else:
        logger.info('bwa index ran on {0}'.format(ref))
        return True


# Returns full path of installed version of bwa using the shell command 'which bwa'
def which_bwa():
    return str(sh.which('bwa')).strip()


# Parent class exposing common functions for BWA commands

class BWA(object):
    # Options that are required for an
    REQUIRED_OPTIONS = ['bwa_path', 'command']
    # regex to detect usage output
    USAGE_REGEX = re.compile('Usage:\s*bwa')

    # args represent the required options(handled in subclasses) the first of these should be the main command e.g. mem
    # kwargs represent any option that has a dash before e.g. -t=4
    def __init__(self, *args, **kwargs):
        # Save args, kwargs for parsing
        self.kwargs = kwargs
        self.args = list(args)
        # Initialize required options values used for main commands e.g. bwa path
        self.required_options_values = []
        # Initialize options list used for optional params e.g. -t=4
        self.options = []
        # Parse kwargs and populates self.required_options_values
        self.required_options()
        # Parse rest of kwargs and populates self.options
        self.compile_bwa_options()
        # This needs to be implemented in subclass
        self.required_args()

    def required_args(self):
        raise NotImplementedError('This function needs to be implemented by subclasses')

    # Parse out REQUIRED_OPTIONS from kwargs and set them in self.required_options_values
    def required_options(self):
        # bwa_path required for all operations
        self.kwargs['bwa_path'] = which_bwa()

        try:
            # Build up the values in order they appear in REQUIRED_OPTIONS
            for op in self.REQUIRED_OPTIONS:
                self.required_options_values.append(self.kwargs[op])
                # Remove from kwargs
                del self.kwargs[op]
        except KeyError as e:
            # Detects if a parameter is missing
            raise ValueError('{0} is a required parameter'.format(op))

    # Parses remaining kwargs and adds them to options list
    # Assumes REQUIRED_OPTIONS since they have been previously processed by required_options()
    # @returns list of options aka [k1, v1, k2, v2...]
    def compile_bwa_options(self):
        # Build up self.options from kwargs
        for op, val in self.kwargs.items():
            # None type values should be ignored
            if val is None or val == '':
                continue
            # Append dash to option
            self.options.append('-' + op)
            # Options should all be strings(just being passed to command line anyways)
            val = str(val)
            # True false values only have option
            if val.lower() not in ('true', 'false'):
                self.options.append(val)

    # Parse stderr output to find if command was executed without errors
    # Since it seems that bwa does not set return codes we have to parse
    # stderr output instead
    #
    # If the following regex is found then the usage statement was printed
    # which indicates a failure of one of the options:
    #     ^Usage:\s+bwa
    #
    # Subclasses need to implement this as well and call this but they need
    # to parse the rest of the output if this returns success in order to
    # tell if the algorithm ran correctly or not
    #
    # @returns 0 if no usage was found, 1 if usage was found
    def bwa_return_code(self, output):
        # Search the output
        m = self.USAGE_REGEX.search(output)

        # If there is a match return 1
        if m:
            logger.warning('Error running BWA command')
            return 2
        # Otherwise return 0
        return 0

    # Wrapper function to make running bwa easier to avoid supplying many arguments
    #
    # @param output_file - The file path to write the sam output to
    # @returns output of self.run_bwa
    def run(self, output_file='bwa_output.sam'):
        return self.run_bwa(self.required_options_values, self.options,
                            self.args, output_file)

    # Executes bwa based on the passed options
    #
    # @param required_options - Should correspond to self.REQUIRED_OPTIONS
    # @param options_list - Full options for bwa as a list (ex. ['-t', '2'])
    # @param args_list - Required arguments that come after options
    # @param output_file - Output location for bwa process output
    #
    # @returns 0 for success, 2 if incorrect options
    #
    # Subclass implementation should return 1 for any other failures
    def run_bwa(self, required_options, options_list, args_list, output_file='bwa_output.sam'):
        # check that bwa path is valid
        if not os.path.exists(required_options[0]):
            raise ValueError('{0} is not a valid bwa path'.format(required_options[0]))

        # Run bwa
        with open(output_file, 'wb') as fh:
            cmd = required_options + options_list + args_list
            logger.debug('Running {0}'.format(cmd))
            p = Popen(cmd, shell=False, stdout=fh, stderr=PIPE)

            # Get the output
            _, process_stderr = p.communicate()

        logger.debug('BWA Output: {0}'.format(process_stderr))

        # Parse the status
        return self.bwa_return_code(process_stderr)

    # Make sure fastapath is a valid path and already has an index
    #
    # @param fastapath - Path to fasta file
    def validate_indexed_fasta(self, fastapath):
        if not is_indexed(fastapath):
            raise ValueError('{0} does not have an index'.format(fastapath))
        if not os.path.exists(fastapath):
            raise ValueError('{0} does not exist'.format(fastapath))

    # Make sure input fastq file is a valid path
    #
    # @param inputpath - Path to fastq file
    def validate_input(self, inputpath):
        if os.path.exists(inputpath):
            try:
                seqio.seqfile_type(inputpath)
            except ValueError:
                logger.warning('{0} is not valid'.format(inputpath))
        else:
            raise ValueError('{0} does not exist'.format(inputpath))


# Subclass exposing the 'bwa index' command
class BWAIndex(BWA):
    def __init__(self, *args, **kwargs):
        # Injects index command and runs super
        kwargs['command'] = 'index'
        super(BWAIndex, self).__init__(*args, **kwargs)

    # 'index' command only requires a valid input fasta file
    def required_args(self):
        if len(self.args) != 1:
            raise ValueError('BWAIndex needs 1 parameter: the input fasta file')
        self.validate_input(self.args[0])
        if seqio.reads_in_file(self.args[0]) == 0:
            raise ValueError('{0} is not a valid file to index'.format(self.args[0]))

    # 'bwa index' return code
    #
    # @param bwa index command output
    # @returns 0 for success, 1 if index file could not be opened
    def bwa_return_code(self, output):
        if '[bwa_index] fail to open file' in output:
            return 1
        return 0

    # Calls super and then removes output file
    def run(self):
        fd, tmpf = tempfile.mkstemp()
        ret = super(BWAIndex, self).run(tmpf)
        os.unlink(tmpf)
        return ret


# Subclass exposing the 'bwa mem' command
class BWAMem(BWA):
    def __init__(self, *args, **kwargs):
        # Injects mem command and runs super
        kwargs['command'] = 'mem'
        # use 4 threads
        # kwargs['t'] = '4'
        super(BWAMem, self).__init__(*args, **kwargs)

    # 'mem' command requires the indexed reference genome and the input fastq file
    def required_args(self):
        # First argument has to be a valid indexed fasta
        self.validate_args(self.args)

    # Validates both the indexed reference genome and the input fastq file
    def validate_args(self, args):
        if len(args) != 2:
            raise ValueError('BWAMem needs 2 paramters: 1. index 2. input fastq file')
        else:
            self.validate_indexed_fasta(self.args[0])
            self.validate_input(self.args[1])

    # 'bwa mem' return code.
    # Output is validated by using regex. 'bwa mem' output has the format:
    # [M::process] read 100 sequences (111350 bp)...
    # [main] Version: 0.7.4-r385
    #
    # @param bwa index command output
    # @returns 0 for success, 1 if index file could not be opened
    def bwa_return_code(self, output):
        read_line_pat = '\[M::process\] read (\d+) sequences \((\d+) bp\)...'
        cpat = re.compile(read_line_pat)

        total_reads = 0
        total_bp = 0
        counts = cpat.findall(output)
        if len(counts) == 0:
            logger.error('BWA Failed as no reads have been processed')
            return 1

        for reads, bps in counts:
            total_reads += int(reads)
            total_bp += int(bps)

        if total_reads == 0 or total_bp == 0:
            logger.error('BWA Failed as no reads and basepairs have been processed')
            return 1
        return super(BWAMem, self).bwa_return_code(output)
