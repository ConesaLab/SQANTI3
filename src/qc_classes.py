#!/usr/bin/env python
import os
from bx.intervals import Interval, IntervalTree
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Union, Optional
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

# TODO: Make this class attributes to directly create the header of the classification file
@dataclass
class myQueryTranscripts:
    id: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[str] = None
    tss_diff: Union[int, str] = "NA"
    tts_diff: Union[int, str] = "NA"
    num_exons: Optional[int] = None
    length: Optional[int] = None
    str_class: str = ""
    subtype: str = "no_subcategory"

    genes: List[str] = field(default_factory=list)
    transcripts: List[str] = field(default_factory=list)

    # Expression & structure
    bite: str = "NA"
    RT_switching: str = "????"
    canonical: str = "NA"
    min_cov: Union[str, float] = "NA"
    min_cov_pos: Union[str, float] = "NA"
    min_samp_cov: Union[str, float] = "NA"
    sd: Union[str, float] = "NA"

    FL: Union[str, int] = "NA"
    FL_dict: Dict[str, int] = field(default_factory=dict)

    nIndels: Union[str, int] = "NA"
    nIndelsJunc: Union[str, int] = "NA"

    proteinID: Optional[str] = None
    protein_length: Union[str,int] = "NA"
    protein_seq: str = "NA"
    psauron_score: Union[str,float] = "NA"
    CDS_start: Union[str, int] = "NA"
    CDS_end: Union[str, int] = "NA"
    CDS_len: Union[str, int] = "NA" #TODO: Perhaps get this from the already present data
    CDS_genomic_start: Union[str, int] = "NA"
    CDS_genomic_end: Union[str, int] = "NA"
    CDS_type: str = "NA"
    is_NMD: Union[str, bool] = "NA"
    coding: str = "non_coding"

    isoExp: Union[str, float] = "NA"
    geneExp: Union[str, float] = "NA"

    refLen: Union[str, int] = "NA"
    refExons: Union[str, int] = "NA"
    refStart: Union[str, int] = "NA"
    refEnd: Union[str, int] = "NA"

    q_splicesite_hit: int = 0
    q_exon_overlap: int = 0
    FSM_class: Optional[str] = None

    percAdownTTS: Union[str, float] = "NA"
    seqAdownTTS: Optional[str] = None

    dist_CAGE: Union[str, float] = "NA"
    within_CAGE: Union[str, bool] = "NA"
    dist_polyA_site: Union[str, float] = "NA"
    within_polyA_site: Union[str, bool] = "NA"
    polyA_motif: str = "NA"
    polyA_dist: Union[str, int] = "NA"
    polyA_motif_found: Union[str, bool] = "NA"
    ratio_TSS: Union[str, float] = "NA"

    tss_gene_diff: Union[str, float] = "NA"
    tts_gene_diff: Union[str, float] = "NA"
    AS_genes: set = field(default_factory=set)

    def ratioExp(self):
        if self.geneExp in ("NA", None, 0):
            return "NA"
        return float(self.isoExp) / float(self.geneExp)

    def geneName(self):
        return "_".join(sorted(set(self.genes)))

    def get_total_diff(self):
        if self.tss_diff == "NA" or self.tts_diff == "NA":
            return "NA"
        return abs(int(self.tss_diff)) + abs(int(self.tts_diff))

    def get_orf_size(self):
        if self.coding == 'coding':
            size = abs(self.CDS_genomic_end - self.CDS_genomic_start) + 1
        else:
            size = "NA"
        return size
    def __str__(self):
        return f"{self.id}: {self.gene_name()} ({self.str_class})"

    def as_dict(self):
        base = self.__dict__.copy()
        base["ratio_exp"] = self.ratioExp()
        base["gene_name"] = self.geneName()
        base["ORF_length"] = self.get_orf_size()
        for sample, count in self.FL_dict.items():
            base[f"FL.{sample}"] = count
        return base

@dataclass
class myQueryProteins:
    def __init__(self, cds_start, cds_end, protein_length, protein_seq=None, proteinID="NA", psauron_score="NA", cds_type="NA"):
        """
        Initialize a myQueryProteins object.

        Parameters
        ----------
        cds_start : int
            The start position of the CDS (coding sequence) on the transcript.
        cds_end : int
            The end position of the CDS on the transcript.
        orf_length : int
            The length of the ORF (open reading frame).
        orf_seq : str, optional
            The amino acid sequence of the ORF (default is None).
        proteinID : str, optional
            The protein ID associated with the ORF (default is "NA").
        psauron_score : float, optional
            The Psauron score for the ORF (default is "NA").
        orf_type : str, optional
            The type of the ORF (default is "NA").
        """
        self._validate_input(cds_start, cds_end,protein_length)

        self.protein_length  = protein_length
        self.cds_start   = cds_start       # 1-based start on transcript
        self.cds_end     = cds_end         # 1-based end on transcript (stop codon), ORF is seq[cds_start-1:cds_end].translate()
        self.cds_genomic_start = None      # 1-based genomic start of ORF, if - strand, is greater than end
        self.cds_genomic_end = None        # 1-based genomic end of ORF
        self.protein_seq     = protein_seq
        self.proteinID   = proteinID
        self.psauron_score = psauron_score # Psauron score for the ORF, if available
        self.cds_type = cds_type           # Type of ORF according to TD2 classification ("complete","5prime_partial",etc...)
        self.cds_length = cds_end - cds_start + 1

    def _validate_input(self, cds_start, cds_end, cds_length, cds_seq=None, proteinID="NA"):
        if cds_start > cds_end:
            raise ValueError("CDS start must be less than CDS end.")
        if cds_length < 0:
            raise ValueError("CDS length must be non-negative.")
        if cds_start < 0 or cds_end < 0:
            raise ValueError("CDS start and end must be non-negative.")
        if cds_seq is not None and not isinstance(cds_seq, str):
            raise ValueError("CDS sequence must be a string.")
        if not isinstance(proteinID, str):
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
