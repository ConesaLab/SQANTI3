import subprocess
import yaml, os
from src.qc_argparse import qc_argparse
from src.filter_argparse import filter_argparse
from src.rescue_argparse import rescue_argparse


def sqanti_path(filename):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),"..",filename)


def create_config(config_path):
    main_args = get_shared_args()
    config = {
        "main" : main_args,
        "qc": get_parser_specific_args_simple(qc_argparse(),main_args),
        "filter": get_parser_specific_args_complex(filter_argparse(),main_args),
        "rescue": get_parser_specific_args_complex(rescue_argparse(),main_args)
    }
    with open(config_path, "w") as f:
        yaml.dump(config, f,sort_keys=False)
    print(f"Config file created at {config_path}")

def get_shared_args():
    return {
        "refGTF": "",
        "refFasta": "",
        "cpus": 4,
        "dir": "sqanti3_results",
        "output": "isoforms"
    }

def get_parser_specific_args_simple(parser,shared_args):
    parser_args = {"enabled": True, "options": {}}
    for action in parser._actions:

        if action.dest not in shared_args:
            if action.default == "==SUPPRESS==":
                continue
            else:
                parser_args["options"][action.dest] = action.default if action.default is not None else "" 
    return parser_args

def get_parser_specific_args_complex(parser,shared_args):
    parser_args = {"enabled": True, "options": {}}
    subparsers = parser._subparsers._group_actions[0].choices
    parser_args["options"]["common"] = {}
    i = 1
    subparsers_names = list(subparsers.keys())
    for subparser_name, subparser in subparsers.items():
        parser_args["options"][subparser_name] = get_parser_specific_args_simple(subparser,shared_args)
        if i != 1:
            parser_args["options"][subparser_name]["enabled"] = False
        i += 1
    options_to_move = []
    
    # Iterate over a copy of the dictionary
    for option, value in parser_args["options"][subparsers_names[0]]["options"].copy().items():
        if option in parser_args["options"][subparsers_names[1]]["options"]:
            parser_args["options"]["common"][option] = value
            options_to_move.append(option)
    
    # Remove the options from subparsers after iteration
    for option in options_to_move:
        del parser_args["options"][subparsers_names[0]]["options"][option]
        del parser_args["options"][subparsers_names[1]]["options"][option]
    
        
    return parser_args  

def format_options(options):
    """Convert a dictionary of options into a command-line argument string."""
    return ' '.join(f'--{key} {value}' for key, value in options.items() if value not in ['',False])

def run_sqanti_module(cmd):
    print(f"Running: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR during SQANTI3 module execution")