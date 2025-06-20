import subprocess
import sys, shutil
import yaml, os
from src.qc_argparse import qc_argparse
from src.filter_argparse import filter_argparse
from src.rescue_argparse import rescue_argparse
from src.logging_config import main_logger,save_module_logger_info


def sqanti_path(filename):
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..",filename))

def check_conda():
    GFFREAD_PROG = "gffread"

    if shutil.which(GFFREAD_PROG) is None:
        main_logger.error(f"Cannot find executable {GFFREAD_PROG}. Abort!")
        main_logger.error(f"Did you activate SQANTI3's conda environment?")
        raise SystemExit(1)

def create_config(config_path,options,level):
    """
    Create a YAML configuration file. 
    It uses the default values from the parsers, unless the user has specified any of them.
    """
    main_args = get_shared_args(level)
    user_options = None
    config = {
        "main" : main_args,
        "qc": get_parser_specific_args_simple(qc_argparse(),main_args),
        "filter": get_parser_specific_args_complex(filter_argparse(),main_args),
        "rescue": get_parser_specific_args_complex(rescue_argparse(),main_args)
    } 
    if options is not None:
        user_options = get_user_options(options,list(flatten_dict(config).keys()))
        replace_value(config,user_options)
    
    config = set_default_values(config,user_options)

    with open(config_path, "w") as f:
        yaml.dump(config, f,sort_keys=False)
    main_logger.info(f"Config file created at {config_path}")

def get_shared_args(level):
    return {
        "refGTF": "",
        "refFasta": "",
        "cpus": 4,
        "dir": "sqanti3_results",
        "output": "isoforms",
        "log_level": level
    }

def get_user_options(options, sqanti_options):
    """Parse the -a options into a dictionary and convert numeric values."""
    options_dict = {}
    for option in options:
        key, value = option.split('=')
        if key not in sqanti_options:
            main_logger.warning(f"Option '{key}' not found in the default configuration.")
            main_logger.warning("Would you like to continue? (y/n)")
            answer = input().strip().lower()
            if answer != 'y':
                main_logger.error("Aborting...")
                sys.exit(1)
            else:
                continue
            continue
        # Try to convert the value to an integer or float
        if value.isdigit():
            value = int(value)
        else:
            try:
                value = float(value)
            except ValueError:
                pass  # Keep the value as a string if it cannot be converted
        options_dict[key] = value
    return options_dict

def replace_value(default_dict, user_config):
    for key, value in default_dict.items():
        if isinstance(value, dict):
            replace_value(value, user_config)
        else:
            if key in user_config:
                default_dict[key] = user_config[key]
    return default_dict

def generate_default_path(config, filename):
    return f"{config['main']['dir']}/{config['main']['output']}{filename}"

def set_default_values(config,user_options):
    user_options = user_options or {}
    
    default_values = {
        'filter': {
            'sqanti_class': '_classification.txt',
            'filter_isoforms': '_corrected.fasta',
            'filter_gtf': '_corrected.gtf',
            'filter_faa': '_corrected.faa'
        },
        'rescue': {
            'filter_class': '_RulesFilter_result_classification.txt' if config['filter']['options']['rules']['enabled'] else '_ML_result_classification.txt',
            'rescue_isoforms': '_corrected.fasta',
            'rescue_gtf': '.filtered.gtf',
            'random_forest': '_randomforest.RData'
        }
    }

    for section, options in default_values.items():
        for key, default_filename in options.items():
            if key not in user_options:
                if section != 'qc':
                   # Avoids including the protein sequences if the user has set skipORF
                   if key == 'filter_faa' and config['qc']['options']['skipORF']:
                       config[section]['options']['common'][key] = generate_default_path(config, default_filename)
                   else:
                       config[section]['options']['common'][key] = generate_default_path(config, default_filename)
                else:
                     config[section]['options'][key] = generate_default_path(config, default_filename)
    return config

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
    
    # Check with arguments are common in both subparsers
    for option, value in parser_args["options"][subparsers_names[0]]["options"].copy().items():
        if option in parser_args["options"][subparsers_names[1]]["options"]:
            parser_args["options"]["common"][option] = value
            options_to_move.append(option)
    
    # Remove the common options from subparsers after iteration
    for option in options_to_move:
        del parser_args["options"][subparsers_names[0]]["options"][option]
        del parser_args["options"][subparsers_names[1]]["options"][option]
    
    return parser_args  

def flatten_dict(d):
    """Flatten a nested dictionary without changing key names."""
    items = []
    for k, v in d.items():
        if isinstance(v, dict):
            items.extend(flatten_dict(v).items())
        else:
            items.append((k, v))
    return dict(items)

def format_options(options):
    """Convert a dictionary of options into a command-line argument string."""
    return ' '.join(f'--{key}' if value is True or value == "True" or value == "true" else f'--{key} {value}' for key, value in options.items() if value not in ['',False])

def validate_user_options(user_options, valid_keys):
    """Check if any user option is not in the list of valid keys."""
    for key in user_options:
        if key not in valid_keys:
            main_logger.warning(f"Option '{key}' not found in the default configuration.")

def modify_options(options,user_options):
    user_options = get_user_options(user_options, list(options.keys()))
    for key, _ in options.items():
        if key in user_options:
            options[key] = user_options[key]
    return options

def run_step(step,config,dry_run, user_options):
    commands = {
        "qc": f"{sys.executable} {sqanti_path('sqanti3_qc.py')} {{options}}",
        "filter": f"{sys.executable} {sqanti_path('sqanti3_filter.py')} {{type}} {{options}}",
        "rescue": f"{sys.executable} {sqanti_path('sqanti3_rescue.py')} {{type}} {{options}}"
    }
    main_opt = config.get("main", {})
    if step == "qc":
        options = main_opt | config[step].get("options", "")
        if user_options is not None:
            modify_options(options,user_options)

        cmd = commands[step].format(options = format_options(options))
        
    else:
        first_subparser=False
        for subparser, subparser_args in config[step].get("options", {}).items():
            if subparser == "common":
                options = main_opt | subparser_args
                if step == "filter": # In filter we don't need the reference files
                    options.pop("refGTF")
                    options.pop("refFasta")
            else:
                if subparser_args["enabled"]:
                    if first_subparser:
                        main_logger.error("Both rules and machine learning are enabled. Only one can be used at a time.")
                        sys.exit(1)
                    first_subparser = True
                    options = options | subparser_args.get("options", {})
                    if user_options is not None:
                        modify_options(options,user_options)
                    cmd = commands[step].format(type = subparser, options = format_options(options))
    main_logger.info(f"Running {step.upper()}")
    if dry_run:
        main_logger.info(f"{cmd}")
    else:
        #save_module_logger_info(f"{main_opt['dir']}/logs",step,main_opt['log_level'],sqanti_path('src/data/module_logger_config.json'))
        run_sqanti_module(cmd)

def run_sqanti_module(cmd):
    main_logger.info(f"{cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        step = cmd.split(" ")[1].split("/")[-1].split("_")[-1].split(".py")[0]
        main_logger.error(f"There was an unexpected issue during SQANTI3 {step.upper()} module execution")
        sys.exit(1)

def run_step_help(step):
    help = {
        "qc": f"{sys.executable} {sqanti_path('sqanti3_qc.py')} -h",
        "filter": f"{sys.executable} {sqanti_path('sqanti3_filter.py')} {{type}} -h",
        "rescue": f"{sys.executable} {sqanti_path('sqanti3_rescue.py')} {{type}} -h"
    }
    if step != "qc":
        filter_type = input(f"Which {step} option do you want to see the help from (rules/ml)? ")
        if filter_type != "rules" and filter_type != "ml":
            main_logger.error("Invalid option")
            return
        help[step] = help[step].format(type=filter_type)
    subprocess.check_call(help[step], shell=True)