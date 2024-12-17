import sys,os,pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)
from src.parsers import reference_parser


