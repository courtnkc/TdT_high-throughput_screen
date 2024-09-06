# Utility functions for analysis of TdT NGS data for TURTLES system
# - See Jupyter notebooks in ./Notebooks directory for example usage and
# - generation of figures in the paper.
#
# Author: Jonathan Strutz (jonathanstrutz2021@u.northwestern.edu)
#
# Note: I copied two functions, closure and clr, from skbio.stats.composition
# to this file since I had issues installing skbio on Windows.

import math
import os
import random
import re

import numpy as np
import pandas as pd
from Bio import SeqIO

def get_pct_len_base_pair_counts(data_dir, num_bins=1000, seq_len_req=None,
                                 filename_end='trimmed.fq', degen=0,
                                 cutoff=0.0, p_remove=0.0):
    """Read R1 fastq files in subdirectories in data_dir. In each file, count
    the number of A, C, G, and T at each position. Then normalize counts by
    the sequence's length by binning.

    Parameters
    ----------
    data_dir : str
        Filepath to folder which contains subfolders which each have an R1
        fastq file.
    num_bins : int (default: 1000)
        The number of bins to divide each sequence up into.
    seq_len_req : int (default: 1000)
        If provided, only sequences of length seq_len_req will be processed.
    filename_end : str (default: 'trimmed.fq')
        Suffix of the fastq filenames to read. Should be the same for every
        condition.
    degen : int (default: 0)
        Number of degenerate bases at the 3' end of the 5' primer. These bases
        get cut out when the sequences are read.
    cutoff : numeric (default: 0.0)
        How many bases to cut off at the end of the sequence. Also gets rid of
        sequences that are less than cutoff length. Can be a decimal value, in
        which case, on average that many bases will be cut off. For example,
        if cutoff=2.5, 2 bases will be cut off 50% of the time and 3 bases will
        be cut off 50% of the time. All sequences less than 4 bases long would
        be thrown out.
    p_remove : float (default: 0.0)
        Randomly discard this proportion of sequences (0 for none, 1 for all).
        Useful for investigating how much data is required for good prediction.

    Returns
    -------
    counts_dict : dict
        Contains an array of counts at each bin (indexed) for each condition.
        Structure is {'directory_name': {'Base': [#, #, ..., #]}}.
    """

    counts_dict = {}

    bins = np.arange(0, 1 + 1 / num_bins, 1 / num_bins)

    for directory in os.listdir(data_dir):
        if os.path.isfile(directory):
            continue
        print(directory)
        for filename in os.listdir(data_dir + directory):
            print(filename)
            if filename.endswith(filename_end) and 'R1' in filename:
                data_file = filename

        filepath = data_dir + directory + '/' + data_file
        if 'R1_Trimmed_Reads' not in filepath and 'No_TdT' not in filepath:

            counts_dict[directory] = {base: np.zeros(num_bins)
                                      for base in ['A', 'C', 'G', 'N', 'T']}

            print('Loading', directory, '\n')
            with open(filepath) as infile:
                fasta_sequences = SeqIO.parse(infile, 'fastq')

                n = 0
                for fasta in fasta_sequences:
                    seq = list(str(fasta.seq))
                    n += 1

                    if p_remove:
                        remove_seq = (random.random() < p_remove)
                        if remove_seq:
                            continue

                    if degen:
                        if len(seq) > degen:
                            seq = seq[degen:]
                        else:
                            continue

                    if cutoff:
                        if cutoff.is_integer and len(seq) > cutoff:
                            seq = seq[:-int(cutoff)]
                        elif not cutoff.is_integer and len(seq) > (cutoff + 1):
                            floor = int(math.floor(cutoff))
                            ceil = int(math.ceil(cutoff))
                            p_floor = cutoff - floor

                            use_floor = (random.random() < p_floor)

                            if use_floor:
                                seq = seq[:-floor]
                            else:
                                seq = seq[:-ceil]
                        else:
                            continue

                    seq_fracs = [i / len(seq) for i in range(len(seq))]

                    if not seq_len_req or len(seq) == seq_len_req:

                        bin_ind = np.digitize(seq_fracs, bins) - 1
                        bin_ind = np.append(bin_ind, num_bins - 1)

                        inds = {base: [] for base in ['A', 'C', 'G', 'T', 'N']}
                        prev_base = seq[0]
                        i_init = 0
                        if len(seq) > 1:
                            for i, base in enumerate(seq[1:]):
                                i = i + 1
                                if base != prev_base:
                                    inds[prev_base].append((bin_ind[i_init],
                                                            bin_ind[i]))
                                    i_init = i

                                if i == len(seq) - 1:
                                    inds[base].append((bin_ind[i_init],
                                                       num_bins))

                                prev_base = seq[i]

                        else:
                            inds[seq[0]].append((0, num_bins))

                        for base in ['A', 'C', 'G', 'T', 'N']:
                            for (i1, i2) in inds[base]:
                                counts_dict[directory][base][i1:i2] += 1

                print('Read', str(n), 'sequences...\n')

    return counts_dict



def read_seqs(data_dir, filename_end='trimmed.fq', degen=0, cutoff=0.0,
              seq_len_req=None, p_discard=0.0, cond_text=None, exact_cond_text=False):
    """Read trimmed R1 fastq files in subdirectories in data_dir. Options for
    removing degenerate bases in primer, cutting off bases at the end of each
    sequence, only reading sequences of a specified length, and for randomly
    discarding a given proportion of sequences are available.

    Parameters
    ----------
    data_dir : str
        Filepath to folder which contains subfolders which each have an R1
        fastq file. read_seqs loops through these subfolders to find these
        files.
    filename_end : str (default: 'trimmed.fq')
        Suffix of the fastq filenames to read. Should be the same for every
        condition.
    degen : int (default: 0)
        Number of degenerate bases at the 3' end of the 5' initiator. This
        number of bases will get cut out from the beginning of each sequence.
    cutoff : float (default: 0.0)
        How many bases to cut off at the end of the sequence. Also gets rid of
        sequences that are less than cutoff length. Can be a decimal value, in
        which case, on average that many bases will be cut off. For example,
        if cutoff=2.5, 2 bases will be cut off 50% of the time and 3 bases will
        be cut off 50% of the time. All sequences less than 3 bases long would
        be thrown out. Any sequences of length 0 after this step are thrown
        out as well.
    seq_len_req : int (default: None)
        If provided, only sequences of length seq_len_req will be processed.
        This length is the sequence length after applying degen and cutoff.
    p_discard : float (default: 0.0)
        Randomly discard this proportion of sequences (0 for none, 1 for all).
        Useful for investigating how much data is required for good prediction.
    cond_text : list of str (default: None)
        If given, searches for these strings within each directory name. If
        none are present, ignores that directory completely. Useful for
        only reading sequences for specific conditions.
    exact_cond_text : bool(default:False)
        If True, matches directory name exactly (using cond_text). If false, searches
        within directory name for cond_text. Discards that condition if not found.
    
    Returns
    -------
    seqs : dict
        Key is directory (condition) name, value is a list of further
        processed sequences given input options.

    See Also
    --------
    cutoff_float
    get_norm_len_base_counts

    """

    if degen:
        if isinstance(degen, float):
            if degen.is_integer():
                degen = int(degen)
            else:
                raise ValueError('degen should be a whole number.')

        if degen < 0:
            raise ValueError('degen should be non-negative.')

    if cutoff:
        if isinstance(cutoff, int):
            cutoff = float(cutoff)

        if cutoff < 0:
            raise ValueError('cutoff should be non-negative.')

    if p_discard:
        if isinstance(p_discard, int):
            p_discard = float(p_discard)

        if p_discard < 0 or p_discard > 1:
            raise ValueError('p_discard should be between 0.0 and 1.0, \
                              inclusive.')

    if cond_text:
        if isinstance(cond_text, str):
            cond_text = [cond_text]

    seqs = {}

    for directory in os.listdir(data_dir):
        if os.path.isfile(directory):
            # Just in case there are files in the data_dir, skip them
            continue

        data_file = None
        for filename in os.listdir(data_dir + directory):
            if filename.endswith(filename_end) and 'R1' in filename:
                data_file = filename

        if not data_file:
            raise FileNotFoundError('File with {} and \'R1\' in filename not '
                                    'found'.format(filename_end))

        if cond_text:
            discard_cond = True
            for cond in cond_text:
                if exact_cond_text and cond == directory:
                    discard_cond = False
                elif not exact_cond_text and cond in directory:
                    discard_cond = False
            if discard_cond:
                continue

        filepath = data_dir + directory + '/' + data_file

        n_seqs = 0  # Need counter since we are using a generator
        seqs[directory] = []
        print('Loading', directory)
        for fasta in SeqIO.parse(filepath, 'fastq'):
            seq = list(str(fasta.seq))
            n_seqs += 1

            if p_discard:
                discard_seq = (random.random() < p_discard)
                if discard_seq:
                    continue

            if degen:
                if len(seq) >= degen:
                    seq = seq[degen:]
                else:
                    continue

            if cutoff:
                if len(seq) < cutoff:
                    continue
                elif cutoff.is_integer():
                    seq = seq[:-int(cutoff)]
                else:
                    # Cutoff is not a whole number
                    seq = cutoff_float(seq, cutoff)
                    if len(seq) == 0:
                        continue

            if seq_len_req and len(seq) != seq_len_req:
                continue

            # If at this point, seq has been processed and passed all checks
            seqs[directory].append(seq)

        print('Read', str(n_seqs), 'sequences...\n')

    return seqs


def cutoff_float(seq, cutoff):
    """Probabilistically cuts off a whole number of bases from the end of a
    sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence to cut bases off of.
    cutoff : float
        Non-whole number of bases to cut off, on average.

    Returns
    -------
    seq : str
        Processed sequence with bases cut off of the end.

    """

    floor = int(math.floor(cutoff))
    ceil = int(math.ceil(cutoff))
    p_floor = cutoff - floor

    use_floor = (random.random() < p_floor)

    if use_floor:
        seq = seq[:-floor]
    else:
        seq = seq[:-ceil]

    return seq


def get_seq_intervals(seqs, positions):
    """Get bases at specified locations within sequences.

    Parameters
    ----------
    seqs : dict
        Key is directory (condition) name, value is a list of further
        processed sequences given input options.
    positions : list of int tuples and/or negative ints
        Base positions to parse each seq. For example, [(1, 5), (11, 15)]
        would return a list of two dicts, one for seqs from bp 1-5, the other
        for seqs from bp 11-15. If just an int, goes from that position to the
        end of the sequence, e.g. -5 would get the last 5 bases. 1-indexed and
        bounds are inclusive.

    Returns
    -------
    seq_intervals : list of dicts
        Key is condition, value is list of new seqs.

    See Also
    --------
    read_seqs

    """

    seq_intervals = [{} for _ in range(len(positions))]

    for directory in seqs:
        for i, pos in enumerate(positions):
            if isinstance(pos, tuple):
                seq_intervals[i][directory] = [seq[(pos[0] - 1): pos[1]] for
                                               seq in seqs[directory]]
            elif isinstance(pos, int):
                if pos > 0:
                    pos -= 1
                seq_intervals[i][directory] = [seq[(pos):] for seq in
                                               seqs[directory]]
            else:
                raise Exception('Positions should be ints and/or tuples of '
                                'positive integers.')

    return seq_intervals


def get_norm_len_base_counts(seqs, num_bins=1000):
    """For each condition in seqs dict, count the number of A, C, G, and T at
    each position. Then normalize counts by the sequence's length by binning.

    Parameters
    ----------
    seqs : dict
        Key is directory (condition) name, value is a list of further
        processed sequences given input options.
    num_bins : int (default: 1000)
        The number of bins to divide each sequence up into.

    Returns
    -------
    counts_dict : dict
        Contains an array of counts at each bin (indexed) for each condition.
        Structure is {'directory_name': {'Base': [#, #, ..., #]}}.

    See Also
    --------
    read_seqs

    """

    counts_dict = {}

    bases = ['A', 'C', 'G', 'N', 'T']

    bins = np.arange(0, 1 + 1 / num_bins, 1 / num_bins)

    for directory in seqs:

        counts_dict[directory] = {base: np.zeros(num_bins) for base in bases}

        for seq in seqs[directory]:

            inds = bin_seq(seq, bins)

            # Given bin interval (e.g. bins 10 to 20) for a given base in this
            # sequence, add 1 to those bins for that base in counts_dict
            for base in bases:
                for (i1, i2) in inds[base]:
                    counts_dict[directory][base][i1:i2] += 1

        print(directory, 'processed\n')

    return counts_dict


def bin_seq(seq, bins):
    """Divides the bases of a sequence into the number of bins specified.
    Returns the bin indices for each base.

    Parameters
    ----------
    seq : list
        Ordered list of bases in sequence.
    bins : np.arange
        Array of bin indices (e.g. array([0, 1, 2, ..., 997, 998, 999])).

    See Also
    --------
    get_norm_len_base_counts

    Notes
    -----
    This function has been optimized for binning which is why we precalculate
    intervals rather than counting bases in each individual bin on the fly.

    """

    num_bins = len(bins)

    seq_fracs = [i / len(seq) for i in range(len(seq))]

    bin_ind = np.digitize(seq_fracs, bins, right=True)

    inds = {base: [] for base in ['A', 'C', 'G', 'T', 'N']}
    prev_base = seq[0]
    i_init = 0

    if len(seq) > 1:
        for i, base in enumerate(seq[1:]):
            i = i + 1
            if base != prev_base:
                inds[prev_base].append((bin_ind[i_init], bin_ind[i]))
                i_init = i

            if i == len(seq) - 1:
                inds[base].append((bin_ind[i_init], num_bins))

            prev_base = seq[i]

    else:
        inds[seq[0]].append((0, num_bins))

    return inds


def calc_norm_len_base_pcts(counts_dict, exclude_n=True):
    """Takes a dict of counts and converts them into percents of A, C, G, and
    T at each bin. Percents reported as fractions (0 --> 1 space).

    Parameters
    ----------
    counts_dict : dict
        Contains an array of counts at each bin (indexed) for each condition.
        Structure is {'directory_name': {'Base': [#, #, ..., #]}}.
    exclude_n : bool (default: True)
        If True, calculate A, C, G, and T percent relative to total A, C, G,
        and T. Othwerise, calculate A, C, G, T, and N percent relative to total
        number of all base calls.

    Returns
    -------
    pcts_dict : dict
        Contains an array of percents at each bin (indexed) for each condition.
        Structure is {'directory_name': {'Base': [#, #, ..., #]}}.

    """

    pcts_dict = {}

    if exclude_n:
        bases = ['A', 'C', 'G', 'T']
    else:
        bases = ['A', 'C', 'G', 'T', 'N']

    for condition in counts_dict:
        total_counts = sum([np.array(counts_dict[condition][base])
                            for base in bases])
        total_counts = np.array(total_counts)

        pcts_dict[condition] = {}
        for base in bases:
            counts = np.array(counts_dict[condition][base])
            pcts_dict[condition][base] = np.divide(counts, total_counts)

    return pcts_dict


def get_total_base_pcts(seqs, exclude_n=True):
    """Read R1 fastq files in subdirectories in data_dir. In each file, count
    the total numbers of A, C, G, and T added across all sequences. Then,
    normalize to percent.

    Parameters
    ----------
    seqs : dict
        Key is directory (condition) name, value is a list of further
        processed sequences given input options.
    exclude_n : bool (default: True)
        If True, calculate A, C, G, and T percent relative to total A, C, G,
        and T. Othwerise, calculate A, C, G, T, and N percent relative to total
        number of all base calls.

    Returns
    -------
    pcts_dict : dict
        Contains a dict of percents for A, C, G, T, and N (unless exclude_n is
        True). {'directory_name': {'A': #, 'C': #, 'G': #, 'T': #}}.

    See Also
    --------
    read_seqs

    """

    pcts_dict = {}

    if exclude_n:
        bases = ['A', 'C', 'G', 'T']
    else:
        bases = ['A', 'C', 'G', 'T', 'N']

    for directory in seqs:

        base_n = {base: 0 for base in bases}

        for seq in seqs[directory]:

            for base in seq:
                if base in bases:
                    base_n[base] += 1

        total = sum([base_n[base] for base in bases])
        pcts_n = {base: base_n[base] / total for base in bases}

        pcts_dict[directory] = pcts_n

    print('Processed all sequences\n\n')

    return pcts_dict


def calc_aitchison_distance(pcts_dict):
    """Transforms compositional data into Aitchison space.

    Parameters
    ----------
    pcts_dict : dict
        Key is condition (directory), value is dict {base: percent}. Percent is
        in fractional form (0 to 1).

    Returns
    -------
    clr_data : dict
        Key is condition (directory), value is dict {base: Aitchison distance}.

    See Also
    --------
    get_total_base_percents

    """

    total_counts = {}
    bases = ['A', 'C', 'G', 'T']

    for condition in pcts_dict:
        # CLR transformation requires %A, C, G, T to be a column vector
        total_counts[condition] = \
            np.column_stack([pcts_dict[condition][base] for base in bases])

    clr_data = {}
    for condition in total_counts:
        clr_data[condition] = {}
        clr_matrix = clr(total_counts[condition])
        for i, base in enumerate(bases):
            clr_data[condition][base] = clr_matrix[i] #[i]?

    return clr_data


def closure(mat):
    """Copied from skbio.stats.composition.closure since I had issues
    installing the skbio package on Windows. Performs closure to ensure that
    all elements add up to 1.

    Parameters
    ----------
    mat : array_like
       a matrix of proportions where
       rows = compositions
       columns = components

    Returns
    -------
    array_like, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1

    Raises
    ------
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import closure
    >>> X = np.array([[2, 2, 6], [4, 4, 2]])
    >>> closure(X)
    array([[ 0.2,  0.2,  0.6],
           [ 0.4,  0.4,  0.2]])

    """

    mat = np.atleast_2d(mat)

    if np.any(mat < 0):
        raise ValueError("Cannot have negative proportions")

    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")

    if np.all(mat == 0, axis=1).sum() > 0:
        raise ValueError("Input matrix cannot have rows with all zeros")

    mat = mat / mat.sum(axis=1, keepdims=True)

    return mat.squeeze()


def clr(mat):
    r"""Copied from skbio.stats.composition.clr since I had issues installing
    the skbio package on Windows. Performs centre log ratio transformation.
    This function transforms compositions from Aitchison geometry to the real
    space. The :math:`clr` transform is both an isometry and an isomorphism
    defined on the following spaces

    :math:`clr: S^D \rightarrow U`

    where :math:`U=\{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in
        \mathbb{R}^D\}`

    It is defined for a composition :math:`x` as follows:
    .. math::
        clr(x) = \ln\left[\frac{x_1}{g_m(x)}, \ldots, \frac{x_D}{g_m(x)}\right]

    where :math:`g_m(x) = (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric
        mean of :math:`x`.

    Parameters
    ----------
    mat : array_like, float
       a matrix of proportions where
       rows = compositions and
       columns = components
    Returns
    -------
    numpy.ndarray
         clr transformed matrix
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])

    """

    mat = closure(mat)
    lmat = np.log(mat)
    gm = lmat.mean(axis=-1, keepdims=True)

    return (lmat - gm).squeeze()


def generate_aitch_df(pcts_dict, clr_data, condition_dict, rep_dict,
                      zero_control_conds, one_control_conds):
    """Take dNTP frequency data and Aitchison distance data and then calculate
    aitchison distance. Format results in a long-form dataframe for plotting.

    Parameters
    ----------
    pcts_dict : dict
        Key is condition (directory), value is dict {base: [#, #, ..., #, #]}
        where each # is a percent and index is bin. Percents are in fractional
        form (0 to 1).
    clr_data : dict
        Key is condition (directory), value is dict {base: [#, #, ..., #, #]}
        where each # is a percent transformed into Aitchison space and index is
        bin.
    condition_dict : dict
        Maps directory name to more readable condition name. Format is
        {directory: condition_name}.
    rep_dict : dict
        Maps directory name to replicate number. Format is {directory: #}.
    zero_control_conds : list
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the zero-signal control. Used for signal calculation.
    one_control_conds : list
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the one-signal control. Used for signal calculation.

    Returns
    -------
    data : pd.DataFrame
        Contains information on dNTP frequency as well as the transformed
        Aitchison value at each bin for each condition, replicate pair.
        Contains the following columns: 'Directory', 'Condition', 'Replicate',
        'Bin Number', 'Aitch Dist (from 0)', 'Aitch Dist (from 1)', 'A %
        Aitch', 'C % Aitch', 'G % Aitch', 'T % Aitch', 'A % Aitch Diff from 0',
        'C % Aitch Diff from 0' 'G % Aitch Diff from 0', 'T % Aitch Diff from
        0', 'A % Aitch Diff from 1, 'C % Aitch Diff from 1', 'G % Aitch Diff
        from 1', 'T % Aitch Diff from 1', 'A %', 'C %', 'G %', 'T %'.

    See Also
    --------
    calc_norm_len_base_pcts : can generate pcts_dict
    calc_aitchison_distance : can generate clr_data
    calc_signal : adds signal column to df

    """

    data = pd.DataFrame()
    bases = ['A', 'C', 'G', 'T']

    directory_col = []
    condition_col = []
    replicate_col = []
    bin_num_col = []
    aitch_dist0_col = []
    aitch_dist1_col = []
    a_ad_col = []
    c_ad_col = []
    g_ad_col = []
    t_ad_col = []
    a_diff0_col = []
    c_diff0_col = []
    g_diff0_col = []
    t_diff0_col = []
    a_diff1_col = []
    c_diff1_col = []
    g_diff1_col = []
    t_diff1_col = []
    a_col = []
    c_col = []
    g_col = []
    t_col = []

    zero_cond_clr_means = {base: np.mean([clr_data[cond][base] for cond in
                           zero_control_conds], axis=0) for base in bases}

    one_cond_clr_means = {base: np.mean([clr_data[cond][base] for cond in
                          one_control_conds], axis=0) for base in bases}

    for directory in pcts_dict:

        a_ad_col += list(clr_data[directory]['A'])
        c_ad_col += list(clr_data[directory]['C'])
        g_ad_col += list(clr_data[directory]['G'])
        t_ad_col += list(clr_data[directory]['T'])

        diffs0 = [abs(clr_data[directory][base] - zero_cond_clr_means[base])
                  for base in bases]
        diffs1 = [abs(clr_data[directory][base] - one_cond_clr_means[base])
                  for base in bases]

        a_diff0_col += list(diffs0[0])
        c_diff0_col += list(diffs0[1])
        g_diff0_col += list(diffs0[2])
        t_diff0_col += list(diffs0[3])

        a_diff1_col += list(diffs1[0])
        c_diff1_col += list(diffs1[1])
        g_diff1_col += list(diffs1[2])
        t_diff1_col += list(diffs1[3])

        aitch_dists0 = np.linalg.norm(diffs0, axis=0)
        aitch_dists1 = np.linalg.norm(diffs1, axis=0)

        condition_name = condition_dict[directory]
        condition_rep = rep_dict[directory]

        for i, (x0, x1) in enumerate(zip(aitch_dists0, aitch_dists1)):
            bin_num = i + 1
            bin_num_col.append(bin_num)
            aitch_dist0_col.append(x0)
            aitch_dist1_col.append(x1)

            directory_col.append(directory)
            condition_col.append(condition_name)
            replicate_col.append(condition_rep)

        a_col += list(pcts_dict[directory]['A'])
        c_col += list(pcts_dict[directory]['C'])
        g_col += list(pcts_dict[directory]['G'])
        t_col += list(pcts_dict[directory]['T'])

    data['Directory'] = directory_col
    data['Condition'] = condition_col
    data['Replicate'] = replicate_col
    data['Bin Number'] = bin_num_col
    data['Aitch Dist (from 0)'] = aitch_dist0_col
    data['Aitch Dist (from 1)'] = aitch_dist1_col
    data['A % Aitch'] = a_ad_col
    data['C % Aitch'] = c_ad_col
    data['G % Aitch'] = g_ad_col
    data['T % Aitch'] = t_ad_col
    data['A % Aitch Diff from 0'] = a_diff0_col
    data['C % Aitch Diff from 0'] = c_diff0_col
    data['G % Aitch Diff from 0'] = g_diff0_col
    data['T % Aitch Diff from 0'] = t_diff0_col
    data['A % Aitch Diff from 1'] = a_diff1_col
    data['C % Aitch Diff from 1'] = c_diff1_col
    data['G % Aitch Diff from 1'] = g_diff1_col
    data['T % Aitch Diff from 1'] = t_diff1_col
    data['A %'] = a_col
    data['C %'] = c_col
    data['G %'] = g_col
    data['T %'] = t_col

    return data


def calc_signal(data, zero_control_conds, one_control_conds,
                zero_control_name='0 Control', one_control_name='1 Control'):
    """Takes in long-form dataframe with columns 'Aitch Dist (from 0)', 'Aitch
    Dist (from 1)', 'Directory', 'Condition', and 'Replicate'. Calculates
    signal and creates new 'Signal' column in dataframe. See paper for details.

    Parameters
    ----------
    data : pd.DataFrame
        Contains information on dNTP frequency as well as the transformed
        Aitchison value at each bin for each condition, replicate pair.
        Contains the following columns: 'Directory', 'Condition', 'Replicate',
        'Bin Number', 'Aitch Dist (from 0)', 'Aitch Dist (from 1)', 'A %
        Aitch', 'C % Aitch', 'G % Aitch', 'T % Aitch', 'A % Aitch Diff from 0',
        'C % Aitch Diff from 0' 'G % Aitch Diff from 0', 'T % Aitch Diff from
        0', 'A % Aitch Diff from 1, 'C % Aitch Diff from 1', 'G % Aitch Diff
        from 1', 'T % Aitch Diff from 1', 'A %', 'C %', 'G %', 'T %'.
    zero_control_conds : list
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the zero-signal control. Used for signal calculation.
    one_control_conds : list
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the one-signal control. Used for signal calculation.
    zero_control_name : str (default: '0 Control')
        The name of the 0 control in the 'Condition' column of data.
    one_control_name : str (default: '0 Control')
        The name of the 1 control in the 'Condition' column of data.

    Returns
    -------
    data : pd.DataFrame
        Updated input data DataFrame with 'Aitch Fraction' and 'Signal' column.

    See Also
    --------
    generate_aitch_df

    """

    n_conditions = len(data.Directory.unique())

    # Calculate relative distance between 0 and 1 control averages
    aitch_dists0 = data['Aitch Dist (from 0)']
    aitch_dists1 = data['Aitch Dist (from 1)']
    data['Aitch Fraction'] = aitch_dists0 / (aitch_dists0 + aitch_dists1)

    # Now, renormalize such that 0 and 1 control means are 0 and 1,
    # respectively. Need to do this because 0 and 1 controls are individually
    # compared to the 0 and 1 control means in previous step, resulting in
    # the across-replicate mean relative distance being greater than 0 for 0
    # controls and less than 1 for 1 controls.
    data0 = data.loc[data.Condition == zero_control_name, :]
    data1 = data.loc[data.Condition == one_control_name, :]

    reps0 = data0.Replicate.unique()
    reps1 = data1.Replicate.unique()

    aitch_fracs0_vals = []
    for rep0 in reps0:
        data0_rep0 = data0.loc[data0.Replicate == rep0, :]
        aitch_fracs0_vals.append(data0_rep0['Aitch Fraction'])

    aitch_fracs0_means = np.mean(aitch_fracs0_vals, axis=0)
    aitch_fracs0_col = np.array(list(aitch_fracs0_means) * n_conditions)

    aitch_fracs1_vals = []
    for rep1 in reps1:
        data1_rep1 = data1.loc[data1.Replicate == rep1, :]
        aitch_fracs1_vals.append(data1_rep1['Aitch Fraction'])

    aitch_fracs1_means = np.mean(aitch_fracs1_vals, axis=0)
    aitch_fracs1_col = np.array(list(aitch_fracs1_means) * n_conditions)

    # Calculate final signal
    signal_col = ((data['Aitch Fraction'] - aitch_fracs0_col) /
                  (aitch_fracs1_col - aitch_fracs0_col))

    data['Signal'] = signal_col

    return data


def get_length_dists(seqs, max_len=200):
    """Read R1 fastq files in subdirectories in data_dir. In each file, count
    the number of sequences of each length.

    Parameters
    ----------
    seqs : dict
        Key is directory (condition) name, value is a list of further
        processed sequences given input options.
    max_len : int
        Maximum sequence length you expect to see in all your data.

    Returns
    -------
    len_dists : dict
        Contains the number of sequences for a given length for each condition.
        Structure is {'directory_name': [num_0_bp_seqs, num_1_bp_seqs, etc.]}.'

    See Also
    --------
    read_seqs
    parse_fastq_lengths
    calc_cum_length_dists

    """

    len_dists = {}

    for directory in seqs:

        lens = [0 for _ in range(max_len + 1)]

        for seq in seqs[directory]:
            lens[len(seq)] += 1

        len_dists[directory] = lens

    return len_dists


def calc_cum_length_dists(len_dists):
    """Takes length distribution data and calculates the cumulative length
    distribution from the regular length distribution. Does this for each
    condition, replicate pair (directory).

    Parameters
    ----------
    len_dists : dict
        Contains the number of sequences for a given length for each condition.
        Structure is {'directory_name': [num_1_bp_seqs, num_2_bp_seqs, etc.]}.'

    Returns
    -------
    cum_len_dists: dict
        A dict that stores a list of cumulative percents for each directory.

    See Also
    --------
    get_length_dists

    """

    cum_len_dists = {}

    for directory in len_dists:
        len_dist = np.array(len_dists[directory])
        pct_counts = len_dist / sum(len_dist)

        cum_counts = []
        cum_val = 0
        for pct in pct_counts:
            cum_val += pct
            cum_counts.append(cum_val)

        cum_len_dists[directory] = cum_counts

    return cum_len_dists


def generate_length_df(len_dists, condition_dict, rep_dict,
                       zero_control_conds=None, one_control_conds=None):
    """Normalize length counts (in len_dists) to percents and format into a
    DataFrame for plotting and further analysis.

    Parameters
    ----------
    len_dists : dict
        Contains the number of sequences for a given length for each condition.
        Structure is {'directory_name': [num_1_bp_seqs, num_2_bp_seqs, etc.]}.'
    condition_dict : dict
        Maps directory name to more readable condition name. Format is
        {directory: condition_name}.
    rep_dict : dict
        Maps directory name to replicate number. Format is {directory: #}.
    zero_control_conds : list
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the zero-signal control. Used for signal calculation.
    one_control_conds : list (default: None)
        List of condition (directory) prefixes of those conditions that are to
        be labeled as the one-signal control. Used for signal calculation.

    Returns
    -------
    len_data : pandas.DataFrame
        Contains length distribution data. Columns are: 'Condition',
        'Replicate', 'Signal', 'Length', 'Count', and '% Count'.

    See Also
    --------
    get_len_dists
    generate_aitch_df

    """

    len_data = pd.DataFrame()

    lengths = range(200)  # Assume longest seq is no longer than 200 bases

    lens_col = []
    len_counts_col = []
    norm_counts_col = []
    cond_col = []
    rep_col = []
    control_signal_val_col = []

    for directory in len_dists:
        lens_col += lengths

        len_counts = len_dists[directory]
        len_counts_col += len_counts

        norm_counts = list(np.array(len_counts) / sum(len_counts))
        norm_counts_col += norm_counts

        cond = condition_dict[directory]
        rep = rep_dict[directory]

        cond_col += [cond] * len(lengths)
        rep_col += [rep] * len(lengths)

        if zero_control_conds and directory in zero_control_conds:
            control_signal_val = 0
        elif one_control_conds and directory in one_control_conds:
            control_signal_val = 1
        else:
            control_signal_val = np.nan

        control_signal_val_col += [control_signal_val] * len(lengths)

    len_data['Condition'] = cond_col
    len_data['Replicate'] = rep_col
    len_data['Signal'] = control_signal_val_col
    len_data['Length'] = lens_col
    len_data['Count'] = len_counts_col
    len_data['% Count'] = norm_counts_col

    return len_data


def parse_fastq_lengths(len_dists, expt_time=60, shift=0):
    """For each condition in the overall data_dir, look at each R1 fastq file
    (all sequences, not deduplicated or trimmed based on length), and tally
    the lengths of all the sequences. Then calculate mean, std length.

    Parameters
    ----------
    len_dists : dict
        Contains the number of sequences for a given length for each condition.
        Structure is {'directory_name': [num_1_bp_seqs, num_2_bp_seqs, etc.]}.
    expt_time : numeric (default: 60)
        Length of the experiment (in minutes).
    shift : int (default: 0)
        See enumerate_indices_from_counts().

    Returns
    -------
    averages : pd.DataFrame
        Pandas DataFrame with fastq length information for each condition.
        Columns are 'Directory', 'Mean' (length), 'Std Devs' (length), and
        'Rate (nt/min)'.

    See Also
    --------
    get_length_dists
    calc_switch_bins

    """

    conds = []
    means = []
    stds = []

    for cond in len_dists:
        len_counts = len_dists[cond]
        lengths = enumerate_indices_from_counts(len_counts, shift=shift)
        conds.append(cond)
        means.append(np.mean(lengths))
        stds.append(np.std(lengths))

    averages = pd.DataFrame()
    averages['Directory'] = conds
    averages['Mean'] = means
    averages['Std Devs'] = stds

    averages['Rate (nt/min)'] = np.array(averages['Mean']) / expt_time

    return averages


def enumerate_indices_from_counts(counts, shift=0):
    """Take a list, where result is specified by index and number of those
    results is specified by the list value, and enumerate all values. For
    example, the list [2, 3, 1] would become [0, 0, 1, 1, 1, 2].

    Parameters
    ----------
    counts : list
        List of counts, where index specifies the result being counted.
    shift : int (default: 0)
        Added to the index value to determine result. Most useful when counts
        are 1-indexed rather than 0-indexed (i.e. by setting shift to 1).

    Returns
    -------
    enumerated_list : list
        List of enumerated values.

    """

    enumerated_list = []
    for index, count in enumerate(counts):
        result = index + shift
        results = [result] * count
        enumerated_list += results

    return enumerated_list


def calc_switch_bins(averages, data_df, mode='01'):
    """Assumes that we want to use the normalized 0 --> 1 signal in data_df
    rather than the norm % difference. Assumes 0 --> 1 signal is precalculated
    and tabulated in data_df in 'Signal' column.

    Parameters
    ----------
    averages : pd.DataFrame
        Pandas DataFrame with fastq length information for each condition.
        Columns are 'Directory', 'Mean' (length), 'Std Devs' (length), and
        'Rate (nt/min)'.
    data_df : pd.DataFrame
        Should contain columns 'Directory' and 'Signal'.
    mode : str
        Determines which way we are switching (from 0 to 1 or from 1 to 0).
        Also, will look for this string to mark switch conditions. For example,
        if mode == '01', then all rows in data_df with '01' in Condition will
        be assumed to be switching experiments.

    Returns
    -------
    averages : pd.DataFrame
        Updated input averages DataFrame with 'Switch Bin' column.

    See Also
    --------
    parse_fastq_lengths
    calc_switch_bin
    calc_switch_times

    """

    directories = list(averages.Directory)
    switch_bin_col = []

    for directory in directories:
        sub_data_df = data_df[data_df.Directory == directory]
        signals = np.array(sub_data_df['Signal'])

        switch_bin = calc_switch_bin(signals, mode=mode)
        switch_bin_col.append(switch_bin)

    averages['Switch Bin'] = switch_bin_col

    return averages


def calc_switch_bin(signals, mode='01'):
    """Find two bin interval which contains 0.50 signal. Then interpolate
    between those two bins to find the exact switch bin. For example, if
    signals array is [0, 0.3, 0.7, 1] (four bins), then switch bin would be
    calculated as 1.5.

    Parameters
    ----------
    signals : np.array
        1-D array of signal values.
    mode : str
        Determines which way we are switching (from 0 to 1 or from 1 to 0).

    Returns
    -------
    switch_bin : float
        Index (bin) in signals at which signal crosses 0.50 threshold. This
        value is interpolated between the two bins around 0.50.

    See Also
    --------
    calc_switch_bins

    """

    switch_bin = None

    for index, signal in enumerate(signals):
        if mode == '01' and signal > 0.5:
            # Switch bins found - interpolate between them
            bin_upper = index
            bin_lower = index - 1
            m = signals[bin_upper] - signals[bin_lower]
            x = bin_lower + (0.5 - signals[bin_lower]) / m
            switch_bin = x + 1
            break
        elif mode == '10' and signal < 0.5:
            # Switch bins found - interpolate between them
            bin_upper = index - 1
            bin_lower = index
            m = signals[bin_lower] - signals[bin_upper]
            x = bin_upper + (0.5 - signals[bin_upper]) / m
            switch_bin = x + 1
            break

    if not switch_bin or switch_bin < 0:
        switch_bin = np.nan

    return switch_bin


def calc_switch_times(averages, num_bins=1000, start_control_conds=None,
                      end_control_conds=None, t_expt=60):
    """Calculate switch times for 01 or 10 (does not support more than one
    switch). If no controls for start and end are specified, assumes constant
    rate throughout. However, if controls are specified, assumes constant 0
    rate and constant (but potentially different) 1 rate (see methods in paper
    for details).

    Parameters
    ----------
    averages : pd.DataFrame
        Pandas DataFrame with fastq length information for each condition.
        Columns are 'Directory', 'Mean' (length), 'Std Devs' (length),
        'Rate (nt/min)', and 'Switch Bin'.
    num_bins : int (default: 1000)
        Number of bins sequences were binned into.
    start_control_conds : list (default: None)
        Name of the conditions in averages.Directory that are the 1 control if
        switching from 1 --> 0 or the 0 control if switching from 0 --> 1.
    end_control_conds : list (default: None)
        Name of the conditions in averages.Directory that are the 0 control if
        switching from 1 --> 0 or the 1 control if switching from 0 --> 1.
    t_expt : numeric (default: 60)
        Length of the experiment (in minutes).

    Returns
    -------
    averages : pd.DataFrame
        Updated input averages DataFrame with 'Switch Time' column.

    See Also
    --------
    calc_switch_bins

    """

    if start_control_conds and end_control_conds:
        r_start_df = averages[averages['Directory'].isin(start_control_conds)]
        r_start = r_start_df['Rate (nt/min)'].mean()

        r_end_df = averages[averages['Directory'].isin(end_control_conds)]
        r_end = r_end_df['Rate (nt/min)'].mean()

        switch_times = calc_switch_time_w_rates(averages['Switch Bin'],
                                                num_bins, t_expt, r_start,
                                                r_end)

        averages['Switch Time'] = switch_times

    else:
        averages['Switch Time'] = (averages['Switch Bin'] / num_bins * t_expt)

    return averages


def calc_switch_time_w_rates(switch_bins, num_bins, t_expt, r_start, r_end):
    """Uses Equations 5 and 6 in paper to calculate switch time. Accounts for
    difference in polymerization rates between 0 and 1 conditions.

    Parameters
    ----------
    switch_bins : pd.Series
        Series of switch bin indices.
    t_expt : numeric
        Length of experiment (in minutes).
    r_start : numeric
        Average rate at the beginning of the reaction.
    r_end : numeric
        Average rate at the end of the reaction.

    Returns
    -------
    switch_times : pd.Series
        Series of switch times (minutes).

    """

    alpha = r_start / r_end

    num = alpha * t_expt
    denom = ((num_bins / switch_bins) + alpha - 1)

    switch_times = num / denom

    return switch_times
