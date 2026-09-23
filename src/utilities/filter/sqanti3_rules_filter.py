import os, sys, json
import pandas as pd
from src.module_logging import message,filter_logger
from typing import Optional

junction_related_columns = [
    "RTS_stage",
    "all_canonical",
    "min_cov",
    "min_cov_pos",
    "sd_cov",
    "n_indels_junc",
    "bite",
    "predicted_NMD"
]


def check_min_prevalence(r, sc: str) -> None:
    """ Checks if min_prevalence is a int value greater than 0
    Args:
        r: value to check
        sc(str): structural category
    
    Exits:
        Calls sys.exit(1) if r is not int or value equal or lower to 0
    """
    if type(r) is int:
        if r < 1:
            filter_logger.error(f"Invalid min_prevalence value {r!r} at {sc}, minimum value is 1")
            sys.exit(1)
    else:
        filter_logger.error(f"Invalid min_prevalence value {r!r} at {sc}, integer value required")
        sys.exit(1)



def read_json_rules(json_file):
    """Parse JSON rules file into structured DataFrame format for filtering.
    
    Processes JSON rules containing filtering criteria for different structural
    categories into a dictionary of DataFrames with standardized rule formats.
    
    Args:
        json_file (str): Path to JSON file containing filtering rules. JSON structure should
            have structural categories as keys with lists of rule dictionaries.
            
    Returns:
        dict: Nested dictionary where keys are structural categories, and values are lists
            of DataFrames containing parsed rules with columns:
            - structural_category (str)
            - column (str): Column name from classification data
            - type (str): Rule type (Category/Min_Threshold/Max_Threshold/min_prevalence)
            - rule (str/float): Threshold value or category requirement

    """
    with open(json_file, 'r') as f:
        json_data = json.load(f)

    names_check(json_data,)

    rules_dict = {} 
    for sc, rules in json_data.items():
        rules_dict[sc] = []
        for rule_set in rules:
            rules_table = []
            for col_name, r in rule_set.items():
                if col_name.lower() in ("min_prevalence","prevalence"):   
                    check_min_prevalence(r, sc)
                    rules_table.append(["prevalence", 'min_prevalence', r])
                elif isinstance(r, list):
                    if all(isinstance(x, (int, float)) for x in r):
                        rules_table.append([col_name, 'Min_Threshold', min(r)])
                        rules_table.append([col_name, 'Max_Threshold', max(r)])
                    else:
                        rules_table.append([col_name, 'Category', [str(value).lower() for value in r]])
                elif isinstance(r, (int, float)):
                    rules_table.append([col_name, 'Min_Threshold', r])
                else:
                    rules_table.append([col_name, 'Category', str(r).lower()])
            rules_dict[sc].append(pd.DataFrame(rules_table, columns=['column', 'type', 'rule']))
    if 'rest' not in rules_dict:
        rules_dict['rest'] = []
        filter_logger.warning("No rules defined for 'rest' structural category. Defaulting to no filtering.")

    return rules_dict

def names_check(input_dict):
    """Check if the input dictionary has the required keys.

    Args:
        input_dict (dict): Dictionary containing structural categories as keys

    Exits:
        Calls sys.exit(1) if any required key is missing in the input dictionary.
    """
    required_keys = [
        'full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog',
        'genic_intron', 'genic', 'antisense', 'fusion', 'intergenic'
    ]

    # Collect all invalid keys
    invalid_keys = [key for key in input_dict if key not in required_keys and key != 'rest']
    for key in invalid_keys: # In case there are any invalid keys
        filter_logger.error(
            f"Invalid structural category '{key}' found in rules file. "
            f"Expected one of {required_keys} or 'rest'."
        )
    if invalid_keys: # We only exit once if there are any invalid keys
        sys.exit(1)

def apply_rules(row, force_multiexon, rules_dict):
    """
    Determine if a transcript should be filtered based on defined rules.
    
    Applies filtering rules to a single transcript row from SQANTI3 classification data.
    Returns "Artifact" if any rule is violated, "Isoform" if all rules pass.
    
    Args:
        row (pd.Series): Single row from SQANTI3 classification dataframe
        force_multiexon (bool): If True, automatically filter mono-exonic transcripts
        rules_dict (dict): Parsed rules from read_json_rules()
        
    Returns:
        str: "Artifact" if transcript fails any filter, "Isoform" otherwise
        
    Note:
        Uses OR evaluation, so in case 
    """
    if force_multiexon and row['exons'] == 1:
        return "Artifact"
    if row['structural_category'] not in rules_dict.keys():
        structural_category = 'rest'
    else:
        structural_category = row['structural_category']
    is_isoform = False
    for rules in rules_dict[structural_category]:
        isoform = True

        for _, rule in rules.iterrows():
            column = rule['column']
            rule_type = rule['type']
            rule_value = rule['rule']
            
            if float(row['exons']) == 1 and column in junction_related_columns:
                continue


            # check if it is nan
            try:
                if pd.isna(row[column]):
                    isoform = False
                else:
                    try:
                        if rule_type == 'Category':
                            if isinstance(rule_value, list):
                                if str(row[column]).lower() not in rule_value:
                                    isoform = False
                            else:
                                if str(row[column]).lower() != rule_value:
                                    isoform = False
                        elif rule_type == 'Min_Threshold':
                            if row[column] < rule_value:
                                isoform = False
                        elif rule_type == 'Max_Threshold':
                            if row[column] > rule_value:
                                isoform = False
                        elif rule_type == "min_prevalence":
                            if row[column] < rule_value:
                                isoform = False
                    except TypeError:
                        filter_logger.error(f"Type error for column {column} with value {row[column]}.")
                        filter_logger.error(f"Check if the column you indicated in the rules file is correct.")
                        sys.exit(1)
            except KeyError:
                filter_logger.error(f"Column {column} not found in SQANTI3 classification data.")
                filter_logger.error(f"Perhaps you misspelled the column name in the rules file?")
                sys.exit(1)
        is_isoform = is_isoform or isoform
        # no need to check next branch if this branch is passed, end early
        if is_isoform:
            break
    if is_isoform:
        return "Isoform"
    else:
        return "Artifact"

def get_reasons(row, force_multiexon, rules_dict):
    """Collect detailed reasons for transcript filtering decisions.
    
    Args:
        row (pd.Series): Single row from SQANTI3 classification dataframe
        force_multiexon (bool): Flag to enforce multi-exon filtering
        rules_dict (dict): Parsed rules from read_json_rules()
        
    Returns:
        pd.Series: Contains three elements:
            - isoform: Transcript ID
            - structural_category: Assigned structural category
            - filter_reason: Semicolon-separated list of failed criteria
            
    Note:
        Uses set to avoid duplicate reasons from multiple rule checks
    """
    reasons = set()
    if force_multiexon and row['exons'] == 1:
        reasons.add("Mono-exonic")
    
    structural_category = row['structural_category'] if row['structural_category'] in rules_dict else 'rest'

    for rules in rules_dict[structural_category]:
        for _, rule in rules.iterrows():
            column = rule['column']
            rule_type = rule['type']
            rule_value = rule['rule']

            if float(row['exons']) == 1 and column in junction_related_columns:
                continue

            if pd.isna(row[column]):
                reasons.add(f"NA value in {column}")
                continue

            if rule_type == 'Category':
                if isinstance(rule_value, list):
                    if str(row[column]).lower() not in rule_value:
                        reasons.add(f"{column}: {row[column]}")
                else:
                    if str(row[column]).lower() != rule_value:
                        reasons.add(f"{column}: {row[column]}")
            elif rule_type == 'Min_Threshold':
                if row[column] < rule_value:
                    reasons.add(f"{column}: {row[column]} < {rule_value}")
            elif rule_type == 'Max_Threshold':
                if row[column] > rule_value:
                    reasons.add(f"{column}: {row[column]} > {rule_value}")
            elif rule_type == "min_prevalence":
                if row[column] < rule_value:
                    reasons.add(f"{column}: {row[column]} < {rule_value} Multisample-artifact")
    
    return pd.Series({
        'isoform': row['isoform'], 
        'structural_category': row['structural_category'], 
        'filter_reason': '; '.join(reasons)
    })


def get_highest_min_prevalence(rules_dict: dict) -> Optional[int]:
    """Return highest min_prevalence threshold present in the parsed rules.
    
    Args:
        rules_dict(dict): input rules

    Return: 
        maximum value of the list of prevalences values
    """
    value = max( 
            [ value
            for rule_sets in rules_dict.values()
            for rules in rule_sets
            for value in rules.loc[rules["type"] == "min_prevalence", "rule"]
            ], default=None)
    
    return value



def drop_prevalence_rules(rules_dict: dict) -> dict:
    """Return a copy of the rules without min_prevalence rules.

    A rule set left empty accepts every isoform..

    Args:
        rules_dict(dict): input rules
    
    Return:
        rules_dict without prevalence rules
    """
    rules_dict = {
        sc: [rules[rules["type"] != "min_prevalence"].reset_index(drop=True)
             for rules in rule_sets]

        for sc, rule_sets in rules_dict.items()
    }

    return rules_dict


def check_prevalence_column(classif, highest_min_prevalence: int)-> None:
    """Check that the classification supports the requested min_prevalence rules.

    Args:
        classif: classification table
        highest_min_prevalence (int): maximum value of the list of prevalences values

    Exits:
        Calls sys.exit(1) if the prevalence column is missing or empty, or if a
        threshold exceeds the number of samples.
    """
    if "prevalence" not in classif.columns or classif["prevalence"].isna().all(): 
            
        filter_logger.error("The rules file contains min_prevalence rules, but the classification file has no prevalence values.")
        filter_logger.error("Run SQANTI3 QC with a multi-sample --fl_count file to generate them or at least use release 6.1")
        sys.exit(1)

    n_samples = sum(col.startswith("FL.") for col in classif.columns)
    if highest_min_prevalence > n_samples: 

        filter_logger.error(f"min_prevalence {highest_min_prevalence} is larger than the number of samples in the classification file ({n_samples}).")
        filter_logger.error(f"Use a lower minimum prevalence value equal or lower to the number of samples({n_samples}).")
        sys.exit(1)

        


def rules_filter(sqanti_class, json_file, force_multi_exon, prefix, logger, ignore_prevalence=False):
    """Main function to execute SQANTI3 filtering workflow.
    
    Args:
        sqanti_class (str): Path to SQANTI3 classification file (TSV format)
        json_file (str): Path to JSON file containing filtering rules
        force_multi_exon (bool): If True, exclude all mono-exonic transcripts
        prefix (str): Output filename prefix
        logger (logging.Logger): Configured logger for progress reporting
        ignore_prevalence (bool): If True, min_prevalence rules are removed before filtering.
            SQANTI3 rescue sets it to filter the reference transcriptome.
        
    Output Files:
        Creates three files in current directory:
        - {prefix}_RulesFilter_classification.txt: Full classification with filter results
        - {prefix}_pass_isoforms.txt: List of passing isoforms
        - {prefix}_filtering_reasons.txt: Detailed filtering reasons for artifacts
        
    Example:
        >>> rules_filter("input.tsv", "rules.json", True, "output", logger)
    """
    message("Reading SQANTI3 classification file",logger)
    
    classif = pd.read_csv(sqanti_class, sep="\t", dtype={'chrom': str})

    message("Reading JSON rules",logger)
    
    rules_dict = read_json_rules(json_file)
    highest_min_prevalence = get_highest_min_prevalence(rules_dict)
    if highest_min_prevalence is not None:
        if ignore_prevalence:
            logger.info("Ignoring min_prevalence rules (--ignore_prevalence).")
            rules_dict = drop_prevalence_rules(rules_dict)
        else:
            check_prevalence_column(classif, highest_min_prevalence)

    message("Applying rules to filter isoforms",logger)

    classif['filter_result'] = classif.apply(lambda row: apply_rules(row, force_multi_exon, rules_dict), axis=1)

    inclusion_list = classif[classif['filter_result'] == "Isoform"]['isoform']

    artifacts_classif = classif[classif['filter_result'] == "Artifact"]
    reasons_df = artifacts_classif.apply(lambda row: get_reasons(row, force_multi_exon, rules_dict), axis=1)

    message("Writing results",logger)
    classif.to_csv(os.path.join(f"{prefix}_RulesFilter_classification.txt"), sep='\t', index=False)
    inclusion_list.to_csv(os.path.join(f"{prefix}_pass_isoforms.txt"), sep='\t', index=False, header=False)
    reasons_df.to_csv(os.path.join(f"{prefix}_filtering_reasons.txt"), sep='\t', index=False)

