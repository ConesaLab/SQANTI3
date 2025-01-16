import yaml
from src.qc_argparse import qc_argparse
from src.filter_argparse import filter_argparse
from src.rescue_argparse import rescue_argparse

def create_config(config_path):
    main_args = get_shared_args()
    config = {
        "main" : main_args,
        "qc": get_parser_specific_args_simple(qc_argparse(),main_args),
        "filter": get_parser_specific_args_complex(filter_argparse(),main_args),
        "rescue": ""# get_parser_specific_args(rescue_argparse(),main_args)
    }
    with open(config_path, "w") as f:
        yaml.dump(config, f,sort_keys=False)
    print(f"Config file created at {config_path}")

def get_shared_args():
    return {
        "refGTF": "",
        "refFASTA": "",
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
    for subparser_name, subparser in subparsers.items():
        parser_args["options"][subparser_name] = get_parser_specific_args_simple(subparser,shared_args)
    return parser_args  