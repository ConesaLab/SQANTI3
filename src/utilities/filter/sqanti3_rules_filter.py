import os, sys, json
import pandas as pd
from src.module_logging import message,filter_logger

junction_related_columns = [
    "RTS_stage",
    "all_canonical",
    "min_cov",
    "min_cov_pos",
    "sd_cov",
    "min_intron_length",
    "n_indels_junc",
    "bite",
    "predicted_NMD"
]

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

    names_check(json_data,)

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
    
    return pd.Series({
        'isoform': row['isoform'], 
        'structural_category': row['structural_category'], 
        'filter_reason': '; '.join(reasons)
    })

def rule_failures(classif, rows, column, rule_type, rule_value):
    """Evaluate one rule on a subset of transcripts at once.

    Vectorized counterpart of the per-row checks in apply_rules and
    get_reasons: NA values fail, categories are compared in lower case and
    thresholds are inclusive.

    Args:
        classif (pd.DataFrame): SQANTI3 classification data
        rows (pd.Series): Boolean mask of the transcripts the rule applies to
        column (str): Classification column the rule checks
        rule_type (str): Category, Min_Threshold or Max_Threshold
        rule_value (str/float/list): Accepted value(s) or threshold

    Returns:
        pd.Series: Boolean Series indexed like the selected rows, True where
            the transcript fails the rule

    Exits:
        Calls sys.exit(1) if the column is missing or cannot be compared
        with the rule value.
    """
    if column not in classif.columns:
        filter_logger.error(f"Column {column} not found in SQANTI3 classification data.")
        filter_logger.error(f"Perhaps you misspelled the column name in the rules file?")
        sys.exit(1)
    values = classif.loc[rows, column]
    try:
        if rule_type == 'Category':
            text = values.astype(str).str.lower()
            if isinstance(rule_value, list):
                failed = ~text.isin(rule_value)
            else:
                failed = text != rule_value
        elif rule_type == 'Min_Threshold':
            failed = values < rule_value
        else:
            failed = values > rule_value
    except TypeError:
        filter_logger.error(f"Type error for column {column} with rule {rule_type} {rule_value}.")
        filter_logger.error(f"Check if the column you indicated in the rules file is correct.")
        sys.exit(1)
    result = values.isna() | failed

    return result

def assigned_category(classif, rules_dict):
    """Structural category whose rules apply to each transcript.

    Args:
        classif (pd.DataFrame): SQANTI3 classification data
        rules_dict (dict): Parsed rules from read_json_rules()

    Returns:
        pd.Series: Structural category of each transcript, or 'rest' when
            the rules file has no rules for it
    """
    category = classif['structural_category']
    result = category.where(category.isin(list(rules_dict.keys())), 'rest')

    return result

def apply_rules_table(classif, force_multiexon, rules_dict):
    """Apply the filtering rules to the whole classification at once.

    Gives the same result as apply_rules row by row: a transcript passes
    if it meets every rule of at least one rule set of its structural
    category, junction rules are skipped for mono-exonic transcripts and,
    with force_multiexon, mono-exonic transcripts are always artifacts.

    Args:
        classif (pd.DataFrame): SQANTI3 classification data
        force_multiexon (bool): If True, automatically filter mono-exonic transcripts
        rules_dict (dict): Parsed rules from read_json_rules()

    Returns:
        pd.Series: "Isoform" or "Artifact" for each transcript, indexed like classif
    """
    mono_exonic = classif['exons'].astype(float) == 1
    forced_out = mono_exonic & force_multiexon
    category = assigned_category(classif, rules_dict)
    passed = pd.Series(False, index=classif.index)
    for structural_category, rule_sets in rules_dict.items():
        # As in apply_rules, a rule set is only evaluated on the transcripts
        # that have not passed an earlier set of the same category.
        pending = (category == structural_category) & ~forced_out
        for rules in rule_sets:
            set_passed = pending.copy()
            for column, rule_type, rule_value in rules.itertuples(index=False):
                applies = pending
                if column in junction_related_columns:
                    applies = pending & ~mono_exonic
                if applies.any():
                    failed = rule_failures(classif, applies, column, rule_type, rule_value)
                    set_passed.loc[failed[failed].index] = False
            passed = passed | set_passed
            pending = pending & ~set_passed
    result = passed.map({True: "Isoform", False: "Artifact"})

    return result

def get_reasons_table(artifacts, force_multiexon, rules_dict):
    """Collect the filtering reasons of all artifacts at once.

    Gives the same reasons as get_reasons row by row: every failed rule of
    every rule set of the structural category, without duplicates.

    Args:
        artifacts (pd.DataFrame): Classification rows filtered as artifacts
        force_multiexon (bool): Flag to enforce multi-exon filtering
        rules_dict (dict): Parsed rules from read_json_rules()

    Returns:
        pd.DataFrame: Columns isoform, structural_category and filter_reason
            (semicolon-separated failed criteria), one row per artifact
    """
    mono_exonic = artifacts['exons'].astype(float) == 1
    category = assigned_category(artifacts, rules_dict)
    messages = [pd.Series(dtype=str)]
    if force_multiexon:
        messages.append(pd.Series("Mono-exonic", index=artifacts.index[mono_exonic]))
    for structural_category, rule_sets in rules_dict.items():
        in_category = category == structural_category
        for rules in rule_sets:
            for column, rule_type, rule_value in rules.itertuples(index=False):
                applies = in_category
                if column in junction_related_columns:
                    applies = in_category & ~mono_exonic
                if applies.any():
                    failed = rule_failures(artifacts, applies, column, rule_type, rule_value)
                    values = artifacts.loc[failed[failed].index, column]
                    text = f"{column}: " + values.astype(str)
                    if rule_type == 'Min_Threshold':
                        text = text + f" < {rule_value}"
                    elif rule_type == 'Max_Threshold':
                        text = text + f" > {rule_value}"
                    text[values.isna()] = f"NA value in {column}"
                    messages.append(text)
    reasons = (pd.concat(messages).rename('filter_reason')
               .rename_axis('row').reset_index().drop_duplicates())
    joined = reasons.groupby('row', sort=False)['filter_reason'].agg('; '.join)
    result = pd.DataFrame({
        'isoform': artifacts['isoform'],
        'structural_category': artifacts['structural_category'],
        'filter_reason': joined.reindex(artifacts.index).fillna('')
    })

    return result

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

    message("Applying rules to filter isoforms",logger)

    classif['filter_result'] = apply_rules_table(classif, force_multi_exon, rules_dict)

    inclusion_list = classif[classif['filter_result'] == "Isoform"]['isoform']

    artifacts_classif = classif[classif['filter_result'] == "Artifact"]
    reasons_df = get_reasons_table(artifacts_classif, force_multi_exon, rules_dict)

    message("Writing results",logger)
    classif.to_csv(os.path.join(f"{prefix}_RulesFilter_classification.txt"), sep='\t', index=False)
    inclusion_list.to_csv(os.path.join(f"{prefix}_pass_isoforms.txt"), sep='\t', index=False, header=False)
    reasons_df.to_csv(os.path.join(f"{prefix}_filtering_reasons.txt"), sep='\t', index=False)

