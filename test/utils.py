
import pytest

# Helper function to validate imports
def validate_import(import_path, import_name):
    print(sys.path)
    try:
        __import__(import_path, globals(), locals(), [import_name], 0)
    except ImportError as e:
        pytest.fail(f"ImportError: Failed to import '{import_name}' from '{import_path}'. Error: {e}")

import sys
def validate_import_own(import_path, import_name,system_path):
    sys.path.insert(0, system_path)

    print(import_path)
    try:

        __import__(import_path, globals(), locals(), [import_name], 0)
    except ImportError as e:
        pytest.fail(f"ImportError: Failed to import '{import_name}' from '{import_path}'. Error: {e}")
