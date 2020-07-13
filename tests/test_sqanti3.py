import logging
from pkg_resources import resource_filename
from os.path import exists, join
import unittest
import sqanti3.sqanti3_qc
import tempfile

try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

logging.basicConfig(level=logging.CRITICAL)


class TestSqanti3_qc(unittest.TestCase):
    def setUp(self):
        self.isoforms = resource_filename("tests", "test_data/test_chr13_seqs.fasta")
        self.annotation = resource_filename(
            "tests", "test_data/Homo_sapiens.GRCh38.86.chr13.gtf"
        )
        self.genome = resource_filename(
            "tests", "test_data/Homo_sapiens.GRCh38.dna.chromosome.13.fa"
        )
        self.cage_peaks = resource_filename(
            "tests", "test_data/hg38.cage_peaks.chr13.bed"
        )
        self.polya_list = resource_filename("tests", "test_data/polyA.list")
        self.expression = resource_filename(
            "tests", "test_data/rsemQuantification.chr13.isoforms.results"
        )
        self.fl_count = resource_filename("tests", "test_data/chr13_FL.abundances.txt")
        self.coverage = resource_filename(
            "tests", "test_data/chr13_SR_support.star.SJ.out.tab"
        )
        self.tappas_gff3 = resource_filename(
            "tests", "test_data/tappAS.Homo_sapiens_GRCh38_Ensembl_86.chr13.gff3"
        )

        self.assertTrue(
            exists(self.annotation),
            msg="Cannot find 'test_data/Homo_sapiens.GRCh38.86.chr13.gtf'",
        )
        self.assertTrue(
            exists(self.genome),
            msg="Cannot find 'test_data/Homo_sapiens.GRCh38.dna.chromosome.13.fa'",
        )
        self.assertTrue(
            exists(self.cage_peaks),
            msg="Cannot find 'test_data/hg38.cage_peaks.chr13.bed'",
        )
        self.assertTrue(
            exists(self.polya_list), msg="Cannot find 'test_data/polyA.list'"
        )
        self.assertTrue(
            exists(self.expression),
            msg="Cannot find 'test_data/rsemQuantification.chr13.isoforms.results'",
        )
        self.assertTrue(
            exists(self.fl_count), msg="Cannot find 'test_data/chr13_FL.abundances.txt'"
        )
        self.assertTrue(
            exists(self.coverage),
            msg="Cannot find 'test_data/chr13_SR_support.star.SJ.out.tab'",
        )
        self.assertTrue(
            exists(self.tappas_gff3),
            msg="Cannot find 'test_data/tappAS.Homo_sapiens_GRCh38_Ensembl_86.chr13.gff3'",
        )

    def test_sqanti3_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sqanti3.sqanti3_qc.sqanti3_qc(
                isoforms=self.isoforms,
                annotation=self.annotation,
                genome=self.genome,
                cage_peak=self.cage_peaks,
                polyA_motif_list=self.polya_list,
                expression=self.expression,
                fl_count=self.fl_count,
                coverage=self.coverage,
                # gff3 = self.tappas_gff3
                directory=tmpdir,
                output="test",
                min_ref_len=200,
                force_id_ignore=False,
                aligner_choice="minimap2",
                polyA_peak=None,
                phyloP_bed=None,
                skipORF=True,
                is_fusion=False,
                gtf=False,
                gmap_index=None,
                cpus=8,
                chunks=1,
                sites="ATAC,GCAG,GTAG",
                window=20,
                genename=False,
                version=False,
                skip_report=True,
                isoAnnotLite=False,
                gff3=None,
                doc=f"{tmpdir}/test.params.txt",
                sense="f",
            )
            assert exists(join(tmpdir, "test_junctions.txt"))


if __name__ == "__main__":
    unittest.main()
