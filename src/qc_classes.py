#!/usr/bin/env python
import os
from bx.intervals import Interval, IntervalTree
from collections import defaultdict
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

class myQueryTranscripts:
    def __init__(self, id, tss_diff, tts_diff, num_exons, length, str_class, subtype=None,
                 genes=None, transcripts=None, chrom=None, strand=None, bite ="NA",
                 RT_switching ="????", canonical="NA", min_cov ="NA",
                 min_cov_pos ="NA", min_samp_cov="NA", sd ="NA", FL ="NA", FL_dict={},
                 nIndels ="NA", nIndelsJunc ="NA", proteinID=None,
                 ORFlen="NA", CDS_start="NA", CDS_end="NA",
                 CDS_genomic_start="NA", CDS_genomic_end="NA",
                 ORFseq="NA",
                 is_NMD="NA",
                 isoExp ="NA", geneExp ="NA", coding ="non_coding",
                 refLen ="NA", refExons ="NA",
                 refStart = "NA", refEnd = "NA",
                 q_splicesite_hit = 0,
                 q_exon_overlap = 0,
                 FSM_class = None, percAdownTTS = None, seqAdownTTS=None,
                 dist_CAGE='NA', within_CAGE='NA',
                 dist_polyA_site='NA', within_polyA_site='NA',
                 polyA_motif='NA', polyA_dist='NA',
                 polyA_motif_found='NA', ratio_TSS='NA'):

        self._validate_inputs(id,strand,CDS_start,CDS_end)

        self.id  = id
        self.tss_diff    = tss_diff   # distance to TSS of best matching ref
        self.tts_diff    = tts_diff   # distance to TTS of best matching ref
        self.tss_gene_diff = 'NA'     # min distance to TSS of all genes matching the ref
        self.tts_gene_diff = 'NA'     # min distance to TTS of all genes matching the ref
        self.genes 		 = genes if genes is not None else []
        self.AS_genes    = set()   # ref genes that are hit on the opposite strand
        self.transcripts = transcripts if transcripts is not None else []
        self.num_exons = num_exons
        self.length      = length
        self.str_class   = str_class  	# structural classification of the isoform
        self.chrom       = chrom
        self.strand 	 = strand
        self.subtype 	 = subtype
        self.RT_switching= RT_switching
        self.canonical   = canonical
        self.min_samp_cov = min_samp_cov
        self.min_cov     = min_cov
        self.min_cov_pos = min_cov_pos
        self.sd 	     = sd
        self.proteinID   = proteinID
        self.ORFlen      = ORFlen
        self.ORFseq      = ORFseq
        self.CDS_start   = CDS_start
        self.CDS_end     = CDS_end
        self.coding      = coding
        self.CDS_genomic_start = CDS_genomic_start  # 1-based genomic coordinate of CDS start - strand aware
        self.CDS_genomic_end = CDS_genomic_end      # 1-based genomic coordinate of CDS end - strand aware
        self.is_NMD      = is_NMD                   # (TRUE,FALSE) for NMD if is coding, otherwise "NA"
        self.FL          = FL                       # count for a single sample
        self.FL_dict     = FL_dict                  # dict of sample -> FL count
        self.nIndels     = nIndels
        self.nIndelsJunc = nIndelsJunc
        self.isoExp      = isoExp
        self.geneExp     = geneExp
        self.refLen      = refLen
        self.refExons    = refExons
        self.refStart    = refStart
        self.refEnd      = refEnd
        self.q_splicesite_hit = q_splicesite_hit
        self.q_exon_overlap = q_exon_overlap
        self.FSM_class   = FSM_class
        self.bite        = bite
        self.percAdownTTS = percAdownTTS
        self.seqAdownTTS  = seqAdownTTS
        self.dist_CAGE   = dist_CAGE
        self.within_CAGE = within_CAGE
        self.within_polyA_site = within_polyA_site
        self.dist_polyA_site   = dist_polyA_site    # distance to the closest polyA site (--polyA_peak, BEF file)
        self.polyA_motif = polyA_motif
        self.polyA_dist  = polyA_dist               # distance to the closest polyA motif (--polyA_motif_list, 6mer motif list)
        self.polyA_motif_found = polyA_motif_found  # boolean output for polyA motif
        self.ratio_TSS = ratio_TSS

    def _validate_inputs(self,id, strand,CDS_start,CDS_end):
        if type(id) != str:
            raise ValueError("ID must be provided in string format.")
        if id == "":
            raise ValueError("ID must be a non-empty string.")
        if (strand == "+" and CDS_start > CDS_end ) or (strand == "-" and CDS_start < CDS_end):
            raise ValueError("CDS start must be less than CDS end in the + strand, and greater in the - strand.")


    def get_total_diff(self):
        return abs(self.tss_diff)+abs(self.tts_diff)
    def add_gene(self, gene):
        self.genes.append(gene)

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons):
        self.transcripts = [ref_transcript]
        self.genes = [ref_gene]
        self.tss_diff = tss_diff
        self.tts_diff = tts_diff
        self.refLen = refLen
        self.refExons = refExons

    def geneName(self):
        geneName = "_".join(sorted(set(self.genes))) # If it is not sorted, the order will be random
        return geneName

    def ratioExp(self):
        if self.geneExp == 0 or self.geneExp == "NA":
            return "NA"
        else:
            ratio = float(self.isoExp)/float(self.geneExp)
        return(ratio)

    def CDSlen(self):
        if self.coding == "coding":
            return(str(int(self.CDS_end) - int(self.CDS_start) + 1))
        else:
            return("NA")

    def __str__(self):
        return '\t'.join([str(v) for v in (
            self.chrom, self.strand,
            self.length, self.num_exons, self.str_class,
            "_".join(set(self.genes)), self.id, self.refLen,
            self.refExons, self.tss_diff, self.tts_diff,
            self.subtype, self.RT_switching, self.canonical,
            self.min_samp_cov, self.min_cov, self.min_cov_pos,
            self.sd, self.FL, self.nIndels, self.nIndelsJunc,
            self.bite, self.isoExp, self.geneExp, self.ratioExp(),
            self.FSM_class, self.coding, self.ORFlen, self.CDSlen(),
            self.CDS_start, self.CDS_end, self.CDS_genomic_start,
            self.CDS_genomic_end, self.is_NMD, self.percAdownTTS,
            self.seqAdownTTS, self.dist_CAGE, self.within_CAGE,
            self.dist_polyA_site, self.within_polyA_site,
            self.polyA_motif, self.polyA_dist,
            self.polyA_motif_found, self.ratio_TSS)])


    def as_dict(self):
        d = {'isoform': self.id,
         'chrom': self.chrom,
         'strand': self.strand,
         'length': self.length,
         'exons': self.num_exons,
         'structural_category': self.str_class,
         'associated_gene': self.geneName(),
         'associated_transcript': "_".join(set(self.transcripts)),
         'ref_length': self.refLen,
         'ref_exons': self.refExons,
         'diff_to_TSS': self.tss_diff,
         'diff_to_TTS': self.tts_diff,
         'diff_to_gene_TSS': self.tss_gene_diff,
         'diff_to_gene_TTS': self.tts_gene_diff,
         'subcategory': self.subtype,
         'RTS_stage': self.RT_switching,
         'all_canonical': self.canonical,
         'min_sample_cov': self.min_samp_cov,
         'min_cov': self.min_cov,
         'min_cov_pos': self.min_cov_pos,
         'sd_cov': self.sd,
         'FL': self.FL,
         'n_indels': self.nIndels,
         'n_indels_junc': self.nIndelsJunc,
         'bite': self.bite,
         'iso_exp': self.isoExp,
         'gene_exp': self.geneExp,
         'ratio_exp': self.ratioExp(),
         'FSM_class': self.FSM_class,
         'coding': self.coding,
         'ORF_length': self.ORFlen,
         'ORF_seq': self.ORFseq,
         'CDS_length': self.CDSlen(),
         'CDS_start': self.CDS_start,
         'CDS_end': self.CDS_end,
         'CDS_genomic_start': self.CDS_genomic_start,
         'CDS_genomic_end': self.CDS_genomic_end,
         'predicted_NMD': self.is_NMD,
         'perc_A_downstream_TTS': self.percAdownTTS,
         'seq_A_downstream_TTS': self.seqAdownTTS,
         'dist_to_CAGE_peak': self.dist_CAGE,
         'within_CAGE_peak': self.within_CAGE,
         'dist_to_polyA_site': self.dist_polyA_site,
         'within_polyA_site': self.within_polyA_site,
         'polyA_motif': self.polyA_motif,
         'polyA_dist': self.polyA_dist,
         'polyA_motif_found':self.polyA_motif_found,
         'ratio_TSS' : self.ratio_TSS
         }
        for sample,count in self.FL_dict.items():
            d["FL."+sample] = count
        return d

class myQueryProteins:

    def __init__(self, cds_start, cds_end, orf_length, orf_seq=None, proteinID="NA"):
        self._validate_input(cds_start, cds_end, orf_length)

        self.orf_length  = orf_length
        self.cds_start   = cds_start       # 1-based start on transcript
        self.cds_end     = cds_end         # 1-based end on transcript (stop codon), ORF is seq[cds_start-1:cds_end].translate()
        self.cds_genomic_start = None      # 1-based genomic start of ORF, if - strand, is greater than end
        self.cds_genomic_end = None        # 1-based genomic end of ORF
        self.orf_seq     = orf_seq
        self.proteinID   = proteinID

    def _validate_input(self, cds_start, cds_end, orf_length, orf_seq=None, proteinID="NA"):
        if cds_start > cds_end:
            raise ValueError("CDS start must be less than CDS end.")
        if orf_length < 0:
            raise ValueError("ORF length must be non-negative.")
        if cds_start < 0 or cds_end < 0:
            raise ValueError("CDS start and end must be non-negative.")



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
