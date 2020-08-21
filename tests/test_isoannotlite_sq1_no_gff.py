import logging
from pkg_resources import resource_filename
from os.path import exists, join
import unittest
from sqanti3.utilities.IsoAnnotLite_SQ1 import isoannot
import tempfile

try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

logging.basicConfig(level=logging.CRITICAL)


class TestIsoAnnotLite_SQ1(unittest.TestCase):
    def setUp(self):
        self.classification = resource_filename(
            "tests", "test_data/test_rep_classification.txt"
        )
        self.annotation = resource_filename("tests", "test_data/test_rep_corrected.gtf")
        self.junctions = resource_filename("tests", "test_data/test_rep_junctions.txt")

        self.assertTrue(
            exists(self.classification),
            msg="Cannot find 'test_data/test_rep_classification.txt'",
        )
        self.assertTrue(
            exists(self.annotation),
            msg="Cannot find 'test_data/test_rep_corrected.gtf'",
        )
        self.assertTrue(
            exists(self.junctions),
            msg="Cannot find 'test_data/test_rep_junctions.txt'",
        )

    def test_isoannot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            isoannot(
                corrected=self.annotation,
                classification=self.classification,
                junctions=self.junctions,
                output=f"{tmpdir}/test_tappas_annot.gtf",
            )
            assert exists(join(tmpdir, "test_tappas_annot.gtf"))


if __name__ == "__main__":
    unittest.main()
