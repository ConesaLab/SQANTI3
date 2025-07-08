import os, sys, json
import pandas as pd
from src.module_logging import message,filter_logger


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
            - type (str): Rule type (Category/Min_Threshold/Max_Threshold)
            - rule (str/float): Threshold value or category requirement

    """
    with open(json_file, 'r') as f:
        json_data = json.load(f)

    rules_dict = {} 
    for sc, rules in json_data.items():
        rules_dict[sc] = []
        for rule_set in rules:
            rules_table = []
            for col_name, r in rule_set.items():
                if isinstance(r, list):
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
    
    return pd.Series({
        'isoform': row['isoform'], 
        'structural_category': row['structural_category'], 
        'filter_reason': '; '.join(reasons)
    })

def rules_filter(sqanti_class,json_file,force_multi_exon,prefix,logger):
    """Main function to execute SQANTI3 filtering workflow.
    
    Args:
        sqanti_class (str): Path to SQANTI3 classification file (TSV format)
        json_file (str): Path to JSON file containing filtering rules
        force_multi_exon (bool): If True, exclude all mono-exonic transcripts
        prefix (str): Output filename prefix
        logger (logging.Logger): Configured logger for progress reporting
        
    Output Files:
        Creates three files in current directory:
        - {prefix}_RulesFilter_result_classification.txt: Full classification with filter results
        - {prefix}_inclusion-list.txt: List of passing isoforms
        - {prefix}_filtering_reasons.txt: Detailed filtering reasons for artifacts
        
    Example:
        >>> rules_filter("input.tsv", "rules.json", True, "output", logger)
    """
    message("Reading SQANTI3 classification file",logger)
    
    classif = pd.read_csv(sqanti_class, sep="\t")

    message("Reading JSON rules",logger)
    
    rules_dict = read_json_rules(json_file)

    message("Applying rules to filter isoforms",logger)

    classif['filter_result'] = classif.apply(lambda row: apply_rules(row, force_multi_exon, rules_dict), axis=1)

    inclusion_list = classif[classif['filter_result'] == "Isoform"]['isoform']

    artifacts_classif = classif[classif['filter_result'] == "Artifact"]
    reasons_df = artifacts_classif.apply(lambda row: get_reasons(row, force_multi_exon, rules_dict), axis=1)

    message("Writing results",logger)
    classif.to_csv(os.path.join(f"{prefix}_RulesFilter_result_classification.txt"), sep='\t', index=False)
    inclusion_list.to_csv(os.path.join(f"{prefix}_inclusion-list.txt"), sep='\t', index=False, header=False)
    reasons_df.to_csv(os.path.join(f"{prefix}_filtering_reasons.txt"), sep='\t', index=False)

