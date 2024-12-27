# unctions with no imports within SQANTI3
import bisect
from collections.abc import Iterable
import math
import os

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
    print(i)
    if i == 0:
        return lst[0]-pos
    elif i == len(lst):
        return lst[-1]-pos
    else:
        a, b = lst[i-1]-pos, lst[i]-pos
        if abs(a) < abs(b): return a
        else: return b