import pytest, os,sys
import warnings
from test.utils import run_command_test
import shutil

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, main_path)
from src.commands import GTF_to_genePred


def test_checkRscript_installation():
    """ Test if Rscript is installed and check its version. """
    command = "Rscript --version"
    result = run_command_test(command)
    assert result.returncode == 0, f"Rscript is not installed or not functioning correctly: {result.stderr}"
    expected_version = "4.3.3"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected Rscript version {expected_version}, but got: {result.stdout}")

def test_kallisto():
    """ Test if Kallisto is installed and check its version. """
    command = "kallisto version"
    result = run_command_test(command)
    assert result.returncode != 127, f"Kallisto command not found: {result.stderr}"
    assert result.returncode == 0, f"Kallisto failed to execute properly: {result.stderr}"
    expected_version = "0.48.0"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected Kallisto version {expected_version}, but got: {result.stdout}")

def test_star():
    """ Test if STAR is installed and check its version. """
    command = "STAR --version"
    result = run_command_test(command)
    assert result.returncode != 127, f"STAR command not found: {result.stderr}"
    assert result.returncode == 0, f"STAR failed to execute properly: {result.stderr}"
    expected_version = "2.7.11b"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected STAR version {expected_version}, but got: {result.stdout}")

def test_gmap():
    """ Test if GMAP is installed and check its version. """
    command = "gmap --version"
    result = run_command_test(command)
    assert result.returncode != 127, f"GMAP command not found: {result.stderr}"
    assert result.returncode == 0, f"GMAP failed to execute properly: {result.stderr}"
    expected_version = "2024-11-20"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected GMAP version {expected_version}, but got: {result.stdout}")

def test_minimap2():
    """ Test if Minimap2 is installed and check its version. """
    command = "minimap2 --version"
    result = run_command_test(command)
    assert result.returncode != 127, f"Minimap2 command not found: {result.stderr}"
    assert result.returncode == 0, f"Minimap2 failed to execute properly: {result.stderr}"
    expected_version = "2.28-r1209"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected Minimap2 version {expected_version}, but got: {result.stdout}")

def test_desalt():
    """ Test if deSALT is installed and check its version. """
    command = "deSALT aln"
    result = run_command_test(command)
    assert result.returncode != 127, f"deSALT command not found: {result.stderr}"
    assert result.returncode == 0, f"deSALT failed to execute properly: {result.stderr}"
    # No version checking for the werid way deSALT is installed

def test_ultra():
    """ Test if uLTRA is installed and check its version. """
    command = "uLTRA --version"
    result = run_command_test(command)
    assert result.returncode != 127, f"uLTRA command not found: {result.stderr}"
    assert result.returncode == 0, f"uLTRA failed to execute properly: {result.stderr}"
    expected_version = "0.1"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected uLTRA version {expected_version}, but got: {result.stdout}")

def test_gffread():
    """ Test if gffread is installed and check its version. """
    command = "gffread --version"
    result = run_command_test(command)
    assert result.returncode != 127, f"gffread command not found: {result.stderr}"
    assert result.returncode == 0, f"gffread failed to execute properly: {result.stderr}"
    expected_version = "0.12.7"  # Modify this as per your required version
    if expected_version not in result.stdout:
        warnings.warn(f"Expected gffread version {expected_version}, but got: {result.stdout}")

## Check programs inside of the utilties directory

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..","src","utilities")

def test_gmst(utiltiesPath=utilitiesPath):
    """ Test if can be found """
    gmst_path = os.path.join(utilitiesPath,"gmst","gmst.pl")
    if shutil.which(gmst_path) is None:
        pytest.fail(f"Cannot find executable {gmst_path}. Abort!")

def test_gtf2genepred(utiltiesPath=utilitiesPath):
    """ Test if can be found """
    gtf2genepred_path = os.path.join(utilitiesPath,"gtfToGenePred")
    if shutil.which(gtf2genepred_path) is None:
        pytest.fail(f"Cannot find executable {gtf2genepred_path}. Abort!")
    # Remove the genePred file and run it again
    try:
        os.remove(os.path.join(main_path,"test","test_data","test_isoforms.genePred"))
    except FileNotFoundError:
        pass
    genePred_file = GTF_to_genePred(os.path.join(main_path,"test","test_data","test_isoforms.gtf"))
    assert os.path.exists(genePred_file), f"Cannot find genePred file {genePred_file}. Abort!"
    #os.remove(genePred_file)

def test_isoannotlite(utiltiesPath=utilitiesPath):
    """ Test if can be found """
    isoannotlite_path = os.path.join(utilitiesPath,"IsoAnnotLite_SQ3.py")
    if shutil.which(isoannotlite_path) is None:
        pytest.fail(f"Cannot find executable {isoannotlite_path}. Abort!")

