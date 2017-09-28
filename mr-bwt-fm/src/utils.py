import logging

# setup console logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
# by default sys.stderr stream is used, hence this won't interfere with the print to console used by mapper.
# specifying StreamHandler(sys.stdout) will mix up mapper's output with logging.
ch = logging.StreamHandler()
# For persistent logging, we can use file handler and log to file
# ch = logging.FileHandler('bwa.log', 'w', encoding=None, delay='true')
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('[%(asctime)s] %(filename)s:%(lineno)d %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


# returns the index per bp in the bps count per ref index
def query_bps_count_index(query_bp):
    query_bp = query_bp.upper()
    switcher = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
        'D': 4,
        'N': 5
    }
    return switcher.get(query_bp)