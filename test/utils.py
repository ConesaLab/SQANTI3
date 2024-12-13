
import pytest

# Helper function to validate imports
def validate_import(import_path, import_name):
    print(sys.path)
    try:
        __import__(import_path, globals(), locals(), [import_name], 0)
    except ImportError as e:
        pytest.fail(f"ImportError: Failed to import '{import_name}' from '{import_path}'. Error: {e}")

def validate_import_own(import_path, import_name,system_path):
    import sys
    sys.path.insert(0, system_path)
    try:
        __import__(import_path, globals(), locals(), [import_name], 0)
    except ImportError as e:
        pytest.fail(f"ImportError: Failed to import '{import_name}' from '{import_path}'. Error: {e}")

# Helper function to run a command and return the result
def run_command(command):
    import subprocess
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result