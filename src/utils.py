# unctions with no imports within SQANTI3
import bisect
from collections.abc import Iterable
import hashlib
import lzma
import math
import os
import re

import numpy as np

def is_expressed(count, min_expression):
    """Whether a count is expression in a sample, given --min_expression.

    With min_expression 0 (the default) any count above 0 is expression; otherwise
    the count must reach min_expression. Works on numbers and on pandas objects.
    """
    if min_expression == 0:
        return count > 0
    return count >= min_expression

def mergeDict(dict1, dict2):
    """ Merge dictionaries to collect info from several files"""
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
                dict3[key] = [value , dict1[key]]
    return dict3

def flatten(lis):
     """ Recursively flattens a nested iterable"""
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:
             yield item

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    mean = sum(data)*1. / n  # mean
    var = sum(pow(x - mean, 2) for x in data) / n  # variance
    return math.sqrt(var)  # standard deviation

def find_polyA_motif(genome_seq, polyA_motif_list):
    """    
    Searches for the first occurrence of any polyA motif from a ranked list within a given genomic sequence.

    Args:
        genome_seq (str): The genomic sequence to search for polyA motifs. The sequence must already be oriented.
        polyA_motif_list (list of str): A ranked list of motifs to search for. The function will report the first motif found.

    Returns:
        tuple: A tuple containing:
            - polyA_motif (str): The first polyA motif found in the sequence. If no motif is found, returns 'NA'.
            - polyA_dist (int or str): The distance (in bases) upstream from the end of the sequence where the motif is found. If no motif is found, returns 'NA'.
            - found (str): 'TRUE' if a motif is found, otherwise 'FALSE'.
    """
    for motif in polyA_motif_list:
        i = genome_seq.find(motif)
        if i >= 0:
            return motif, -(len(genome_seq)-i-len(motif)+1), 'TRUE'
    return 'NA', 'NA', 'FALSE'


def get_files_from_dir(directory, extension):
    """ Get all files with a given extension from a directory or a file"""
    if os.path.isfile(directory):
        with open(directory) as f:
            return [line.strip() for line in f]  # Corrected strip method call
    else:
        return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(extension)]
    

def find_closest_in_list(lst, pos):
    """
    Finds the closest value in a sorted list to a given position.

    Args:
        lst (list of int/float): A sorted list of numbers.
        pos (int/float): The position to find the closest value to.

    Returns:
        int/float: The difference between the closest value in the list and the given position.

    Example:
        >>> find_closest_in_list([1, 3, 5, 7], 4)
        -1
    """
    i = bisect.bisect_left(lst, pos)
    if i == 0:
        return lst[0]-pos
    elif i == len(lst):
        return lst[-1]-pos
    else:
        a, b = lst[i-1]-pos, lst[i]-pos
        if abs(a) < abs(b): return a
        else: return b
    
def alphanum_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

    
def calculate_tss(strand, start0, end1):
    """
    Strand aware calculation of the middle of the peak in the bed file
    If the cage peak length is of 1 nucleotide, the average is not calculcated
    """
    if end1 - start0 > 1:
        tss0 = int((start0 + end1) / 2)
        if strand == '+':
            return tss0
        else:
            return tss0 + 1
    else:
        if strand == '+':
            return start0
        else:
            return end1


def hash_intron_chains(keys: list[str], bits: int) -> np.ndarray:
    """Hash intron chain keys to integers of a given number of bits.

    The hash is the leading bits of the MD5 digest of the key, so it is
    the same on every platform and Python version (unlike hash()).

    Args:
        keys: Intron chain keys (see intron_chain_key()).
        bits: Number of bits to keep, from 1 to 64.

    Returns:
        uint64 array with one hash per key, in the order of keys.

    Raises:
        ValueError: If bits is outside 1-64.
    """
    if not 1 <= bits <= 64:
        raise ValueError(f"Hash bits must be between 1 and 64, got {bits}.")
    shift = 64 - bits
    result = np.fromiter(
        (int.from_bytes(hashlib.md5(key.encode()).digest()[:8], "big")
         >> shift for key in keys),
        dtype=np.uint64, count=len(keys))

    return result


def pack_array(values: np.ndarray) -> np.ndarray:
    """Compress a numeric array for the pantranscriptome catalog.

    The bytes of the little-endian values are regrouped by position
    (all the first bytes, then all the second bytes...) before xz
    compression, because the high bytes of small numbers are mostly
    zeros and compress far better together.

    Args:
        values: One-dimensional array of unsigned integers.

    Returns:
        uint8 array with the compressed bytes.
    """
    little = values.astype(values.dtype.newbyteorder("<"), copy=False)
    planes = little.view(np.uint8).reshape(-1, values.dtype.itemsize).T
    compressed = lzma.compress(planes.tobytes(),
                               preset=9 | lzma.PRESET_EXTREME)
    result = np.frombuffer(compressed, dtype=np.uint8)

    return result


def unpack_array(packed: np.ndarray, dtype: str) -> np.ndarray:
    """Undo pack_array().

    Args:
        packed: uint8 array returned by pack_array().
        dtype: Unsigned integer type of the original values, e.g.
            "uint32".

    Returns:
        The original values with the native byte order.
    """
    little = np.dtype(dtype).newbyteorder("<")
    planes = np.frombuffer(lzma.decompress(packed.tobytes()),
                           dtype=np.uint8)
    values = planes.reshape(little.itemsize, -1).T.copy().view(little)
    result = values.ravel().astype(dtype)

    return result
