import logging
from pkg_resources import resource_filename
from os.path import exists
import unittest
from sqanti3.utilities.IsoAnnotLite_SQ1 import createGTFFromSqanti
import ujson as json


logging.basicConfig(level=logging.CRITICAL)


class TestCreateGTFFromSqanti(unittest.TestCase):
    def setUp(self):
        self.gtf = resource_filename(
            "tests", "test_data/example_out/melanoma_chr13_corrected.gtf"
        )
        self.classification = resource_filename(
            "tests", "test_data/example_out/melanoma_chr13_classification.txt"
        )
        self.junctions = resource_filename(
            "tests", "test_data/example_out/melanoma_chr13_junctions.txt"
        )
        self.gff3 = resource_filename(
            "tests", "test_data/Homo_sapiens_GRCh38_Ensembl_86.gff3"
        )
        self.filename = "tappAS_annot_from_SQANTI3.gff3"
        self.filenameMod = f"{self.filename[:-5]}_mod{self.filename[-5:]}"

        self.assertTrue(
            exists(
                resource_filename(
                    "tests", "test_data/test_createGTFFromSqanti/dc_SQexons.json"
                )
            )
        )
        self.assertTrue(
            exists(
                resource_filename(
                    "tests", "test_data/test_createGTFFromSqanti/dc_SQcoding.json"
                )
            )
        )
        self.assertTrue(
            exists(
                resource_filename(
                    "tests", "test_data/test_createGTFFromSqanti/dc_SQtransGene.json"
                )
            )
        )
        self.assertTrue(
            exists(
                resource_filename(
                    "tests", "test_data/test_createGTFFromSqanti/dc_SQstrand.json"
                )
            )
        )

        with open(
            resource_filename(
                "tests", "test_data/test_createGTFFromSqanti/dc_SQexons.json"
            ),
            "r",
        ) as f:
            self.dc_SQexons = json.load(f)

        with open(
            resource_filename(
                "tests", "test_data/test_createGTFFromSqanti/dc_SQcoding.json"
            ),
            "r",
        ) as f:
            self.dc_SQcoding = json.load(f)

        with open(
            resource_filename(
                "tests", "test_data/test_createGTFFromSqanti/dc_SQtransGene.json"
            ),
            "r",
        ) as f:
            self.dc_SQtransGene = json.load(f)

        with open(
            resource_filename(
                "tests", "test_data/test_createGTFFromSqanti/dc_SQstrand.json"
            ),
            "r",
        ) as f:
            self.dc_SQstrand = json.load(f)

    def test_createGTFFromSqanti(self):
        (dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQstrand) = createGTFFromSqanti(
            exons_gtf=self.gtf,
            transcript_classification=self.classification,
            junctions=self.junctions,
            filename=self.filename,
        )

        self.assertEqual(dc_SQexons, self.dc_SQexons)
        self.assertEqual(dc_SQcoding, self.dc_SQcoding)
        self.assertEqual(dc_SQtransGene, self.dc_SQtransGene)
        self.assertEqual(dc_SQstrand, self.dc_SQstrand)


if __name__ == "__main__":
    unittest.main()
