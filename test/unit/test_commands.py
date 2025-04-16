import os,sys,pytest
import logging
main_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../..")
sys.path.append(main_path)
from src.commands import get_aligner_command


# Mock the global command templates
@pytest.fixture(autouse=True)
def mock_cmd_templates(monkeypatch):
    monkeypatch.setattr('src.commands.GMAP_CMD', "gmap -t {cpus} -d {name} -D {dir} --sense-force {i} > {o}")
    monkeypatch.setattr('src.commands.MINIMAP2_CMD', "minimap2 -ax splice --secondary=no -C5 -uf -t {cpus} {g} {i} > {o}")
    monkeypatch.setattr('src.commands.DESALT_CMD', "deSALT aln -t {cpus} {dir} {i} -o {o}")
    monkeypatch.setattr('src.commands.ULTRA_CMD', "uLTRA map -t {cpus} --prefix {prefix} -g {g} -a {a} -i {i} -o {o_dir}")

@pytest.fixture
def default_args():
    return {
        "genome": "genome.fa",
        "isoforms": "iso.fa",
        "annotation": "anno.gtf",
        "outdir": "/out",
        "corrSAM": "out.sam",
        "n_cpu": 4,
        "gmap_index": "/index/gmap"
    }
# Create a fixture for the tester logger
@pytest.fixture
def tester_logger():
    # Define a logger for testing
    logger = logging.getLogger("tester_logger")

    logger.setLevel(logging.INFO)
    # Return the logger
    return logger

def test_gmap_command(default_args, tester_logger):
    cmd = get_aligner_command("gmap", **default_args, logger=tester_logger)
    expected = "gmap -t 4 -d gmap -D /index --sense-force iso.fa > out.sam"
    assert cmd == expected
    expected = "gmap -t 4 -d gmap -D /index --sense-force iso.fa > out.sam"
    assert cmd == expected

def test_minimap2_command(default_args, tester_logger):
    cmd = get_aligner_command("minimap2", **default_args, logger=tester_logger)
    expected = "minimap2 -ax splice --secondary=no -C5 -uf -t 4 genome.fa iso.fa > out.sam"
    assert cmd == expected

def test_desalt_command(default_args, tester_logger):
    cmd = get_aligner_command("deSALT", **default_args, logger=tester_logger)
    expected = "deSALT aln -t 4 /index/gmap iso.fa -o out.sam"
    assert cmd == expected

def test_ultra_command(default_args, tester_logger):
    cmd = get_aligner_command("uLTRA", **default_args, logger=tester_logger)
    expected = "uLTRA map -t 4 --prefix ../out -g genome.fa -a anno.gtf -i iso.fa -o /out/uLTRA_out/"
    assert cmd == expected

def test_invalid_aligner(default_args,caplog,tester_logger):
    caplog.set_level(logging.ERROR)
    with pytest.raises(ValueError):
        get_aligner_command("invalid_aligner", **default_args,logger=tester_logger)
    assert "Unsupported aligner choice: invalid_aligner" in caplog.text

def test_output_capture(default_args, caplog,tester_logger):
    caplog.set_level(logging.INFO)
    get_aligner_command("gmap", **default_args,logger=tester_logger)
    assert "****Aligning reads with GMAP..." in caplog.text