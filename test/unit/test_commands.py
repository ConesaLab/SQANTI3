import os,sys,pytest

main_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../..")
sys.path.append(main_path)
from src.commands import get_aligner_command


# Mock the global command templates
@pytest.fixture(autouse=True)
def mock_cmd_templates(monkeypatch):
    monkeypatch.setattr('src.commands.GMAP_CMD', "gmap -t {cpus} -d {name} -D {dir} {sense} {i} > {o}")
    monkeypatch.setattr('src.commands.MINIMAP2_CMD', "minimap2 -t {cpus} {sense} {g} {i} > {o}")
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
        "gmap_index": "/index/gmap",
        "sense": "--sense-force"
    }

def test_gmap_command(default_args):
    cmd = get_aligner_command("gmap", **default_args)
    expected = "gmap -t 4 -d gmap -D /index --sense-force iso.fa > out.sam"
    assert cmd == expected

def test_minimap2_command(default_args):
    cmd = get_aligner_command("minimap2", **default_args)
    expected = "minimap2 -t 4 --sense-force genome.fa iso.fa > out.sam"
    assert cmd == expected

def test_desalt_command(default_args):
    cmd = get_aligner_command("deSALT", **default_args)
    expected = "deSALT aln -t 4 /index/gmap iso.fa -o out.sam"
    assert cmd == expected

def test_ultra_command(default_args):
    cmd = get_aligner_command("uLTRA", **default_args)
    expected = "uLTRA map -t 4 --prefix ../out -g genome.fa -a anno.gtf -i iso.fa -o /out/uLTRA_out/"
    assert cmd == expected

def test_invalid_aligner(default_args):
    with pytest.raises(ValueError, match="Unsupported aligner choice: invalid_aligner"):
        get_aligner_command("invalid_aligner", **default_args)

def test_output_capture(default_args, capsys):
    get_aligner_command("gmap", **default_args)
    captured = capsys.readouterr()
    assert "****Aligning reads with GMAP..." in captured.out
