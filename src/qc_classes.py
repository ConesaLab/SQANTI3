#!/usr/bin/env python
import os
from bx.intervals import Interval, IntervalTree
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional, Union, get_args, get_origin
from src.module_logging import qc_logger
from src.utils import calculate_tss

class genePredReader(object):
    """
    A class to read gene prediction records from a file.

    Attributes:
    -----------
    f : file object
        The file object for the input file containing gene prediction records.

    Methods:
    --------
    __iter__():
        Returns the iterator object itself.

    __next__():
        Reads the next line from the file, converts it to a genePredRecord object, and returns it.
        Raises StopIteration when the end of the file is reached.
    """
    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)

class genePredRecord(object):
    """
    A class to represent a gene prediction record.

    Attributes:
    -----------
    id : str
        Identifier for the gene prediction record.
    chrom : str
        Chromosome name.
    strand : str
        Strand information ('+' or '-').
    txStart : int
        0-based transcription start position.
    txEnd : int
        0-based transcription end position.
    cdsStart : int
        0-based coding sequence start position.
    cdsEnd : int
        0-based coding sequence end position.
    exonCount : int
        Number of exons.
    exonStarts : list of int
        List of 0-based exon start positions.
    exonEnds : list of int
        List of 1-based exon end positions.
    gene : str, optional
        Gene name (default is None).
    length : int
        Total length of the exons.
    exons : list of Interval
        List of exon intervals.
    junctions : list of tuple
        List of junctions represented as (1-based last base of previous exon, 1-based first base of next exon).

    Methods:
    --------
    segments:
        Returns the list of exon intervals.

    from_line(cls, line):
        Creates an instance of genePredRecord from a tab-separated line.

    get_splice_site(genome_dict, i):
        Returns the donor-acceptor splice site pattern for the i-th junction.
        """
    def __init__(self, id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene=None):
        # Validate inputs
        self._validate_inputs(txStart, txEnd, cdsStart, cdsEnd, strand, exonCount, exonStarts, exonEnds)


        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart         # 1-based start
        self.txEnd = txEnd             # 1-based end
        self.cdsStart = cdsStart       # 1-based start
        self.cdsEnd = cdsEnd           # 1-based end
        self.exonCount = exonCount
        self.exonStarts = exonStarts   # 0-based starts
        self.exonEnds = exonEnds       # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s,e in zip(exonStarts, exonEnds):
            self.length += e-s
            self.exons.append(Interval(s, e))

        # junctions are stored (1-based last base of prev exon, 1-based first base of next exon)
        self.junctions = [(self.exonEnds[i],self.exonStarts[i+1]) for i in range(self.exonCount-1)]

    def _validate_inputs(self, txStart, txEnd, cdsStart, cdsEnd,strand, exonCount, exonStarts, exonEnds):

        if txStart < 0 or txEnd < 0 or cdsStart < 0 or cdsEnd < 0:
            raise ValueError("Transcription and coding start/end positions must be non-negative.")

        if txStart >= txEnd:
            raise ValueError("Transcription start must be less than transcription end.")

        if cdsStart > cdsEnd:
            raise ValueError("CDS start must be less than CDS end.")

        if exonCount <= 0:
            raise ValueError("Exon count must be a positive integer.")

        if len(exonStarts) != exonCount or len(exonEnds) != exonCount:
            raise ValueError("Exon starts and ends must match the exon count.")

        for start, end in zip(exonStarts, exonEnds):
            if start >= end:
                raise ValueError("Exon start positions must be less than exon end positions.")

        if txStart > min(exonStarts) or txEnd < max(exonEnds):
            raise ValueError("Transcript boundaries must encompass all exons.")
    @property
    def segments(self):
        return self.exons


    @classmethod
    def from_line(cls, line):
        raw = line.strip().split('\t')
        return cls(id=raw[0],
                  chrom=raw[1],
                  strand=raw[2],
                  txStart=int(raw[3]),
                  txEnd=int(raw[4]),
                  cdsStart=int(raw[5]),
                  cdsEnd=int(raw[6]),
                  exonCount=int(raw[7]),
                  exonStarts=[int(x) for x in raw[8][:-1].split(',')],  #exonStarts string has extra , at end
                  exonEnds=[int(x) for x in raw[9][:-1].split(',')],     #exonEnds string has extra , at end
                  gene=raw[11] if len(raw)>=12 else None,
                  )

    def get_splice_site(self, genome_dict, i):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        assert 0 <= i < self.exonCount-1

        d = self.exonEnds[i]
        a = self.exonStarts[i+1]

        seq_d = genome_dict[self.chrom].seq[d:d+2]
        seq_a = genome_dict[self.chrom].seq[a-2:a]

        if self.strand == '+':
            return (str(seq_d)+str(seq_a)).upper()
        else:
            return (str(seq_a.reverse_complement())+str(seq_d.reverse_complement())).upper()

@dataclass
class myQueryTranscripts:
    isoform: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[str] = None
    length: Optional[int] = None
    exons: Optional[int] = None

    structural_category: str = ""
    subcategory: str = "no_subcategory"
    FSM_class: Optional[str] = None
    associated_gene: List[str] = field(default_factory=list)
    associated_transcript: List[str] = field(default_factory=list)

    ref_length: Optional[int] = None
    ref_exons: Optional[int] = None

    q_splicesite_hit: Optional[int] = None
    q_exon_overlap: Optional[int] = None

    diff_to_TSS: Optional[int] = None
    diff_to_TTS: Optional[int] = None
    diff_to_gene_TSS: Optional[float] = None
    diff_to_gene_TTS: Optional[float] = None

    RTS_stage: Optional[str] = None
    all_canonical: Optional[str] = None

    min_sample_cov: Optional[float] = None
    min_cov: Optional[float] = None
    min_cov_pos: Optional[float] = None
    sd_cov: Optional[float] = None

    FL: Optional[int] = None

    n_indels: Optional[int] = None
    n_indels_junc: Optional[int] = None
    bite: Optional[str] = None

    iso_exp: Optional[int] = None
    gene_exp: Optional[int] = None
    ratio_exp: Optional[float] = None  # Computed field, may be updated dynamically

    coding: str = "non_coding"
    CDS_length: Optional[int] = None
    protein_length: Optional[int] = None
    CDS_start: Optional[int] = None
    CDS_end: Optional[int] = None
    CDS_genomic_start: Optional[int] = None
    CDS_genomic_end: Optional[int] = None
    psauron_score: Optional[float] = None
    CDS_type: Optional[str] = None
    predicted_NMD: Optional[bool] = None
    
    perc_A_downstream_TTS: Optional[float] = None
    seq_A_downstream_TTS: Optional[str] = None

    dist_to_CAGE_peak: Optional[float] = None
    within_CAGE_peak: Optional[bool] = None
    dist_to_polyA_site: Optional[float] = None
    within_polyA_site: Optional[bool] = None

    polyA_motif: Optional[str] = None
    polyA_dist: Optional[int] = None
    polyA_motif_found: Optional[bool] = None

    protein_seq: Optional[str] = None
    ratio_TSS: Optional[float] = None

    # Extra fields not in FIELDS_CLASS TODO: take them out of the dictionary
    FL_dict: Dict[str, int] = field(default_factory=dict)
    AS_genes: Set[str] = field(default_factory=set)
    genes: Optional[list] = None  # List of genes associated with the isoform
    transcripts: Optional[list] = None  # List of transcripts associated with the isoform
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None
    ref_strand: Optional[str] = None

    def __post_init__(self):
        self._validate_input()
        self._validate_types()
    
    def _validate_input(self):
        if self.isoform == "":
            raise ValueError("Isoform identifier cannot be empty.")
        if self.CDS_start is not None and self.CDS_end is not None:
            if self.CDS_length is None:
                self.CDS_length = self.CDS_end - self.CDS_start + 1
            if self.CDS_start > self.CDS_end:
                raise ValueError("CDS start must be less than CDS end.")
        
    def _validate_types(self):
        for field_name, expected_type in self.__annotations__.items():
            value = getattr(self, field_name)

            # Skip None values unless the field is not Optional
            if value is None:
                origin = get_origin(expected_type)
                if origin is Union and type(None) in get_args(expected_type):
                    continue  # It's Optional[...] and value is None, which is fine
                if expected_type is Optional:
                    continue  # Fallback
                if expected_type is type(None):
                    continue
                # If it's not optional and is None, raise
                raise TypeError(f"'{field_name}' must be of type {expected_type}, but got None")

            if not self._is_instance_of_type(value, expected_type):
                raise TypeError(f"'{field_name}' must be of type {expected_type}, got {type(value).__name__}")
    @staticmethod
    def _is_instance_of_type(value, expected_type):
        origin = get_origin(expected_type)
        args = get_args(expected_type)

        if origin is None:
            return isinstance(value, expected_type)
        if origin in (list, List):
            return isinstance(value, list) and all(isinstance(v, args[0]) for v in value)
        if origin in (set, Set):
            return isinstance(value, set) and all(isinstance(v, args[0]) for v in value)
        if origin in (dict, Dict):
            return isinstance(value, dict) and all(
                isinstance(k, args[0]) and isinstance(v, args[1]) for k, v in value.items()
            )
        # Fallback for unsupported generic types
        return isinstance(value, expected_type)

    def ratioExp(self):
        if self.gene_exp in (None, 0) or self.iso_exp is None:
            return None
        return float(self.iso_exp) / float(self.gene_exp)

    def geneName(self):
        return "_".join(sorted(set(self.genes)))

    def add_gene(self, gene):
        """
        Add a gene to the associated genes list.
        If the gene is already present, it will not be added again.
        """
        if self.genes is None:
            self.genes = [gene]
            return
        if gene not in self.genes:
            self.genes.append(gene)
            return
        
    def get_total_diff(self):
        if self.diff_to_TSS is None or self.diff_to_TTS is None:
            return None
        return abs(self.diff_to_TSS) + abs(self.diff_to_TTS)

    def get_orf_size(self):
        if self.coding == 'coding' and self.CDS_genomic_end is not None and self.CDS_genomic_start is not None:
            try:
                return abs(int(self.CDS_genomic_end) - int(self.CDS_genomic_start)) + 1
            except:
                return None
        return None

    def update(self, attrs: dict):
        for key, val in attrs.items():
            setattr(self, key, val)
    
    # Output methods
    def __str__(self):
        return str([{k: getattr(self, k)} for k in vars(self)])
        # return f"{self.isoform}: {self.geneName()} ({self.structural_category})"

    def as_dict(self):
        base = self.__dict__.copy()
        base["ratio_exp"] = self.ratioExp()
        base["ORF_length"] = self.get_orf_size()
        base["associated_gene"] = self.geneName()
        base["associated_transcript"] = '_'.join(set(self.transcripts))

        # Replace None with "NA"
        for k, v in base.items():
            if v is None:
                base[k] = "NA"

        for sample, count in self.FL_dict.items():
            base[f"FL.{sample}"] = count

        # Eliminate non-report attributes
        non_report_attrs = ['AS_genes','FL_dict','genes','transcripts', 'ref_start', 'ref_end', 'ref_strand']
        for attr in non_report_attrs:
            if attr in base:
                del base[attr]
        return base
    
@dataclass
class myQueryProteins:
    cds_start: int
    cds_end: int
    protein_length: int

    cds_type: Optional[str] = None
    cds_length: Optional[int] = field(init=False)

    protein_seq: Optional[str] = None
    proteinID: Optional[str] = None
    psauron_score: Optional[float] = None

    def __post_init__(self):
        self._validate_input()

        # Automatically compute CDS length
        self.cds_length = self.cds_end - self.cds_start + 1

    def _validate_input(self):
        if self.cds_start > self.cds_end:
            raise ValueError("CDS start must be less than CDS end.")
        if self.cds_start < 0 or self.cds_end < 0:
            raise ValueError("CDS start and end must be non-negative.")
        if self.protein_length < 0:
            raise ValueError("Protein length must be non-negative.")
        if self.protein_seq is not None and not isinstance(self.protein_seq, str):
            raise ValueError("Protein sequence must be a string.")
        if self.proteinID is not None and not isinstance(self.proteinID, str):
            raise ValueError("Protein ID must be a string.")


class CAGEPeak:
    """
    A class to represent and query CAGE (Cap Analysis of Gene Expression) peaks from a BED file.

    Attributes
    ----------
    cage_bed_filename : str
        The filename of the BED file containing CAGE peak data.
    cage_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects containing intervals of peaks.

    Methods
    -------
    __init__(cage_bed_filename):
        Initializes the CAGEPeak object with the given BED filename and reads the BED file to populate the peaks.
    read_bed():
        Reads the BED file and populates the cage_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=10000):
        Queries the CAGE peaks to determine if a given position falls within a peak and calculates the distance to the nearest TSS.
    """
    def __init__(self, cage_bed_filename):
        self._validate_input(cage_bed_filename)
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def _validate_input(self, cage_bed_filename):
        if not cage_bed_filename.endswith('.bed'):
            raise ValueError("CAGE peak file must be in BED format.")
        if not os.path.exists(cage_bed_filename):
            raise FileNotFoundError("CAGE peak file does not exist.")

    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split('\t')
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = calculate_tss(strand, start0, end1)
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (tss0, start0, end1))
    
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
    def find(self, chrom, strand, query, search_window=10000):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within a cage peak>, <nearest distance to TSS>
        If the distance is negative, the query is upstream of the TSS.
        If the query is outside of the peak upstream of it, the distance is NA
        """
        within_peak, dist_peak = 'FALSE', float('inf')
        peaks = self.cage_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for tss0, start0, end1 in peaks:
            # Checks if the TSS is upstream of a peak
            qc_logger.debug(f"Query {query} is within peak {start0}-{end1} on {chrom} strand {strand}")
            if (strand == '+' and start0 > query and end1 > query) or \
            (strand == '-' and start0 < query and end1 < query):
                qc_logger.debug(f"Query {query} is upstream of peak {start0}-{end1} on {chrom} strand {strand}")
                continue
            # Checks if the query is within the peak and the distance to the TSS
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            d = (tss0 - query) * (-1 if strand == '-' else 1)
            w = 'TRUE' if within_out else 'FALSE'

            if not within_peak=='TRUE':
                within_peak, d = w, (tss0 - query) * (-1 if strand=='-' else +1)
                if within_peak == 'TRUE' or abs(d) < abs(dist_peak):
                    dist_peak = d

            else:
                d = (tss0 - query) * (-1 if strand=='-' else +1)
                if abs(d) < abs(dist_peak) and not(w == 'FALSE' and within_peak == 'TRUE'):
                    within_peak, dist_peak = w, d
        if dist_peak == float('inf'):
            dist_peak = 'NA'
        return within_peak, dist_peak

class PolyAPeak:
    """
    A class to represent and query polyA peaks from a BED file.

    Attributes
    ----------
    polya_bed_filename : str
        The filename of the BED file containing polyA peak information.
    polya_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects representing intervals of peaks.

    Methods
    -------
    __init__(polya_bed_filename)
        Initializes the PolyAPeak object with the given BED filename and reads the BED file to populate polyA peaks.
    read_bed()
        Reads the BED file and populates the polya_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=100)
        Queries the polyA peaks to determine if a given position falls within a specified search window of any peak.
    """
    def __init__(self, polya_bed_filename):
        self.polya_bed_filename = polya_bed_filename
        self.polya_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.polya_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            self.polya_peaks[(chrom,strand)].insert(start0, end1, (start0, end1))

    def find(self, chrom, strand, query, search_window=100):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within some distance to polyA>, <distance to closest>
        - if downstream, + if upstream (watch for strand!!!)
        """
        assert strand in ('+', '-')
        within_polyA, dist_polyA = 'FALSE', 'NA'
        hits = self.polya_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for start0, end1 in hits:
            # Checks if the query is within the tail and the distance to the 5'
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            distance = start0 - query if strand == '+' else query - end1

            if within_out:
                within_polyA = 'TRUE'
            if dist_polyA == 'NA' or abs(distance) < abs(dist_polyA):
                dist_polyA = distance

        return within_polyA, dist_polyA
