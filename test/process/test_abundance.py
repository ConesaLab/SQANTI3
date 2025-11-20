from collections import defaultdict
import os, sys, pytest
from typing import Any
from Bio import SeqIO
import pandas as pd

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

# Import functions to test
