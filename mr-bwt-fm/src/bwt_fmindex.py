from itertools import islice, izip_longest

dollar_initial_char = '$'


# This is a crucial part of the optimized SA which uses Timsort to extract integer keys from the given string.
# This works fine on repetitive data such as DNA and works best on data which is already partially sorted, like
# this case
#
# @param input text to create suffix array for
# @returns the character index for each distinct character in input text
def text_to_int_keys(input_text):
    seen = set()
    ls = []
    for e in input_text:
        if e not in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in input_text]


# This is an optimized version of SA which skips the last of the suffix matrix it contains repetitive elements
#
# @param input text to create suffix array for
# @returns suffix array
def suffix_array(input_text):
    n = len(input_text)
    k = 1
    line = text_to_int_keys(input_text)
    while max(line) < n - 1:
        line = text_to_int_keys(
            [a * (n + 1) + b + 1
             for (a, b) in
             izip_longest(line, islice(line, k, None), fillvalue=-1)])
        k <<= 1

    # obtain sorted indices of all suffixes by doing an inverse permutation of sa
    n = len(line)
    ans = [0] * n
    for i in range(n):
        ans[line[i]] = i
    return ans


# Creates bwt for given inout text based on suffix array
#
# @param input text to create bwt for
# @param suffix array
# @returns bwt
def bwt_from_sa(input_text, sa=None):
    bw = []
    # ensure sa exists
    if sa is None:
        sa = suffix_array(input_text)
    for si in sa:
        if si == 0:
            bw.append(dollar_initial_char)
        else:
            bw.append(input_text[si - 1])

    # return bwt list format
    return ''.join(bw)


# Base class managing rank checkpoints and handling rank queries
class FmCheckpoints(object):
    # Scan BWT, creating periodic checkpoints on the fly
    def __init__(self, bw, cpIval=128):
        # checkpoints
        self.cps = {}
        # spacing between checkpoints
        self.cpIval = cpIval
        # checkpoints tally
        tally = {}
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Build the checkpoints
        for i in xrange(0, len(bw)):
            tally[bw[i]] += 1
            if (i % cpIval) == 0:
                for c in tally.iterkeys():
                    self.cps[c].append(tally[c])

    # Returns the number of chars there are in bw up to and including row
    def rank(self, bw, c, row):
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc


# Downsamples the input suffix array by taking only the suffix-array entries for every nth suffix.
#
# @param original suffix array
# @returns downsampled suffix array
def downsample_suffix_array(sa, n=32):
    ssa = {}
    for i in xrange(0, len(sa)):
        if sa[i] % n == 0:
            ssa[i] = sa[i]
    return ssa


# Creates the custom FM-Index for the input text
#
# @param input text
# @returns a tuple containing 4 main data structures making up fm-index:
# -- bwt
# -- downsampled suffix array
# -- rank checkpoints
# -- the first column containing number of occurrences of each character
def make_index(input_text, cpIval=128, ssaIval=32):
    if input_text[-1] != dollar_initial_char:
        input_text += dollar_initial_char  # add dollar if not there already

    # create suffix array
    sa = suffix_array(input_text)
    # create bwt
    bwt = bwt_from_sa(input_text, sa)
    # downsample suffix array
    ssa = downsample_suffix_array(sa, ssaIval)
    # create rank checkpoints
    checkpoints = FmCheckpoints(bwt, cpIval)
    # Calculate no occurrences of each character
    tots = dict()
    for c in bwt:
        tots[c] = tots.get(c, 0) + 1
    # Calculate concise representation of first column
    first_col = {}
    totc = 0
    for c, count in sorted(tots.iteritems()):
        first_col[c] = totc
        totc += count

    return bwt, ssa, checkpoints, first_col


# Return number of occurrences of characters < c
#
# @param fm-index first column
# @param character
# @returns number of occurrences as specified in fm-index first column
def count_occurrences(first_col, character):
    if character not in first_col:
        # (Unusual) case where does not occur in text
        for cc in sorted(first_col.iterkeys()):
            if character < cc: return first_col[cc]
        return first_col[cc]
    else:
        return first_col[character]

# Returns the range of BWM rows having query as a prefix
#
# @param the fm-index
# @param input query text
# @returns the bwm rows where input query text is found in fm-index
def bwm_range(bwt_fmindex, query, mismatches=1):
    bwt, ssa, checkpoints, first_col = bwt_fmindex
    l, r = 0, len(bwt) - 1
    for i in xrange(len(query) - 1, -1, -1):  # from right to left
        l = checkpoints.rank(bwt, query[i], l - 1) + count_occurrences(first_col, query[i])
        r = checkpoints.rank(bwt, query[i], r) + count_occurrences(first_col, query[i]) - 1
        if r < l:
            break
    return l, r + 1

# Returns the offset of a BWM row wrt to original text t
#
# @param the fm-index
# @param bwm row
# @returns the offset of the BWM row in fm-index
def resolve(bwt_fmindex, row):
    bwt, ssa, checkpoints, first_col = bwt_fmindex

    # move left according to character in given BWT row
    def stepLeft(row):
        c = bwt[row]
        return checkpoints.rank(bwt, c, row - 1) + count_occurrences(first_col, c)

    nsteps = 0
    while row not in ssa:
        row = stepLeft(row)
        nsteps += 1
    return ssa[row] + nsteps


# Returns offsets for all occurrences of query based on the supplied BWT FM-Index
#
# @param query to search for
# @param the fm-index
# @param number of mismatches allowed
# @returns the offsets for all occurrences of query searching for
def all_occurrences(query, bwt_fmindex, mismatches=1):
    l, r = bwm_range(bwt_fmindex, query, mismatches)
    return [resolve(bwt_fmindex, x) for x in xrange(l, r)]

# Returns the offset for first occurrence of query based on the supplied BWT FM-Index.
#
# @param query to search for
# @param the fm-index
# @param number of mismatches allowed
# @returns the offsets for first occurrence of query searching for
def first_occurrence(query, bwt_fmindex, mismatches=1):
    l, r = bwm_range(bwt_fmindex, query, mismatches)
    if l >= r:
        return -1  # no occurrence
    else:
        # just first occurrence
        return resolve(bwt_fmindex, l)
