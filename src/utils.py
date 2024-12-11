# unctions with no imports within SQANTI3
from collections.abc import Iterable
import math

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
    :param genome_seq: genomic sequence to search polyA motifs from, must already be oriented
    :param polyA_motif_list: ranked list of motifs to find, report the top one found
    :return: polyA_motif, polyA_dist (how many bases upstream is this found)
    """
    for motif in polyA_motif_list:
        i = genome_seq.find(motif)
        if i >= 0:
            return motif, -(len(genome_seq)-i-len(motif)+1), 'TRUE'
    return 'NA', 'NA', 'FALSE'

