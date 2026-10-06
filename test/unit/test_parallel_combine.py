import os, sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.parallel import concat_tables, concat_sams

HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr22\tLN:50818468\n"

def write(path, content):
    path.write_text(content)
    return str(path)

# concat_sams

def test_concat_sams_keeps_one_header(tmp_path):
    sam1 = write(tmp_path / "0.sam", HEADER + "@PG\tID:minimap2\tCL:chunk0\nPB.1.1\t0\tchr22\t100\t60\t50M\t*\t0\t0\t*\t*\n")
    sam2 = write(tmp_path / "1.sam", HEADER + "@PG\tID:minimap2\tCL:chunk1\nPB.2.1\t0\tchr22\t500\t60\t50M\t*\t0\t0\t*\t*\n"
                                             "PB.2.1\t2048\tchr22\t900\t60\t50H20M\t*\t0\t0\t*\t*\n")
    out = tmp_path / "combined.sam"
    concat_sams([sam1, sam2], str(out))

    lines = out.read_text().splitlines()
    assert [l for l in lines if l.startswith("@")] == HEADER.splitlines() + ["@PG\tID:minimap2\tCL:chunk0"]
    assert [l.split("\t")[0] for l in lines if not l.startswith("@")] == ["PB.1.1", "PB.2.1", "PB.2.1"]

# concat_tables

def test_concat_tables_keeps_one_header(tmp_path):
    t1 = write(tmp_path / "0.txt", "isoform\tvalue\nPB.1.1\t1\n")
    t2 = write(tmp_path / "1.txt", "isoform\tvalue\nPB.2.1\t2\nPB.3.1\t3\n")
    out = tmp_path / "combined.txt"
    assert concat_tables([t1, t2], str(out)) == 3
    assert out.read_text() == "isoform\tvalue\nPB.1.1\t1\nPB.2.1\t2\nPB.3.1\t3\n"

def test_concat_tables_skips_missing_files(tmp_path):
    # e.g. the supplementary report only exists in the chunks that had split alignments
    t2 = write(tmp_path / "1.txt", "isoform\tvalue\nPB.2.1\t2\n")
    out = tmp_path / "combined.txt"
    assert concat_tables([str(tmp_path / "0.txt"), t2], str(out)) == 1
    assert out.read_text() == "isoform\tvalue\nPB.2.1\t2\n"

def test_concat_tables_no_inputs_writes_nothing(tmp_path):
    out = tmp_path / "combined.txt"
    assert concat_tables([str(tmp_path / "0.txt"), str(tmp_path / "1.txt")], str(out)) == 0
    assert not out.exists()
