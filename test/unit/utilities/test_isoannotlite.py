import pytest,sys,os

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)
from src.utilities.IsoAnnotLite_SQ3 import createGTFFromSqanti

TEST_DATA = os.path.join(main_path, "test/test_data")


@pytest.fixture
def classification(tmp_path):
    # The test classification stores CDS_start and CDS_end as floats
    # ("5384.0"); SQANTI3 6.x writes them as integers.
    with open(os.path.join(TEST_DATA,
                           "isoforms/test_isoforms_classification.tsv")) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    cds_columns = [rows[0].index("CDS_start"), rows[0].index("CDS_end")]
    for row in rows[1:]:
        for i in cds_columns:
            if row[i].endswith(".0"):
                row[i] = row[i][:-2]
    filename = str(tmp_path / "classification.txt")
    with open(filename, "w") as f:
        f.write("".join("\t".join(row) + "\n" for row in rows))
    return filename


@pytest.fixture
def tappas_gtf(tmp_path, classification):
    filename = str(tmp_path / "tappas_annot.gtf")
    createGTFFromSqanti(
        os.path.join(TEST_DATA, "isoforms/test_isoforms.gtf"),
        classification,
        os.path.join(TEST_DATA, "junctions/test_junctions.txt"),
        {},
        filename,
    )
    with open(filename) as f:
        lines = f.read().split("\n")[:-1]
    return lines


def test_no_empty_lines(tappas_gtf):
    # updateGTF() stops reading at the first line without 9 fields, so an
    # empty line drops every exon and feature written after it.
    assert "" not in tappas_gtf


def test_exon_rows(tappas_gtf):
    with open(os.path.join(TEST_DATA, "isoforms/test_isoforms.gtf")) as f:
        gtf_exons = [line.split("\t") for line in f
                     if line.split("\t")[2] == "exon"]
    exons = [line.split("\t") for line in tappas_gtf
             if line.split("\t")[2] == "exon"]
    assert len(exons) == len(gtf_exons)
    assert all(len(fields) == 9 for fields in exons)
    assert [fields[8] for fields in exons] == \
        ["Chr=" + fields[0] for fields in gtf_exons]
