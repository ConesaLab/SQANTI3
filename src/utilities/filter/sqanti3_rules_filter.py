import os, sys, json
import pandas as pd
from src.module_logging import message,filter_logger
from typing import Optional
from src.config import MIN_DETECTION_COUNT

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
    """Check that a min_prevalence value is a valid threshold.

    The value is either one integer, applied to all the samples (or to every
    group when a counts design is given), or an object with one integer per
    experimental group. Group names are checked against the design later, in
    resolve_prevalence_rules(), because the rules are read before it.

    Args:
        r: value to check
        sc(str): structural category

    Exits:
        Calls sys.exit(1) if r is not an integer >= 1, or a non-empty object
        whose values are all integers >= 1.
    """
    if isinstance(r, dict):
        if not r:
            filter_logger.error(f"Empty min_prevalence object at {sc}, at least one group is required")
            sys.exit(1)
        values = r.items()
    else:
        values = [(None, r)]

    for group, value in values:
        where = sc if group is None else f"{sc} (group {group!r})"
        if type(value) is not int:
            filter_logger.error(f"Invalid min_prevalence value {value!r} at {where}, integer value required")
            sys.exit(1)
        if value < 1:
            filter_logger.error(f"Invalid min_prevalence value {value!r} at {where}, minimum value is 1")
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
            - rule (str/float/dict): Threshold value or category requirement.
            min_prevalence keeps the JSON value (int or {group: int}); with a
            counts design, resolve_prevalence_rules() maps it to prevalence columns.

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


def prevalence_thresholds(rule_value) -> dict:
    """Return a min_prevalence value as {prevalence column: threshold}.

    An integer is a threshold on the prevalence column computed by QC over all
    the samples. A dict has already been resolved against a counts design by
    resolve_prevalence_rules() and maps prevalence.<group> columns to thresholds.

    Args:
        rule_value: possible dict of values or int with mininimum prevalence value
    
    Return:
        dict with value or values with minimum prevalence values
    """
        
    if type(rule_value) is dict:
        output = rule_value
    else:
        output = {"prevalence": rule_value}

    return output


def passes_prevalence(row, rule_value) -> bool:
    """Evaluate a min_prevalence rule on one isoform.

    OR between groups: the rule passes if any prevalence column reaches its
    threshold. A NA value means the isoform was detected in no sample (e.g. it
    is missing from the --fl_count file), so that column never reaches it.
    The reference transcriptome, with no expression data, never gets here:
    rescue drops these rules with --ignore_prevalence.

    Args:
        row: single row of the classification table
        rule_value (int or dict): min_prevalence value, either an integer checked
            on the prevalence column, or {prevalence column: threshold} resolved
            by resolve_prevalence_rules()

    Returns:
        bool: True if the isoform passes the rule, False if not
    """
    thresholds = prevalence_thresholds(rule_value)
    
    passes = any(not pd.isna(row[col]) and row[col] >= m
                 for col, m in thresholds.items())

    return passes


def has_prevalence_rules(rules_dict: dict) -> bool:
    """Check whether the parsed rules contain any min_prevalence rule.

    Args:
        rules_dict (dict): parsed rules from read_json_rules(), with structural
            categories as keys and lists of rule DataFrames as values

    Returns:
        bool: True if any rule set contains a min_prevalence rule, False otherwise
    """
    has_rules = any((rules["type"] == "min_prevalence").any()
                    for rule_sets in rules_dict.values()
                    for rules in rule_sets)

    return has_rules


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

            if rule_type == "min_prevalence":
                if not passes_prevalence(row, rule_value):
                    isoform = False
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

            if rule_type == "min_prevalence":
                if not passes_prevalence(row, rule_value):
                    thresholds = prevalence_thresholds(rule_value)
                    na_cols = [col for col in thresholds if pd.isna(row[col])]
                    reasons.update(f"NA value in {col}" for col in na_cols)
                    failed = ", ".join(f"{col}: {row[col]} < {m}"
                                       for col, m in thresholds.items() if col not in na_cols)
                    if failed:
                        reasons.add(f"{failed} Multisample-artifact")
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


def _reject_duplicate_groups(pairs):
    """object_pairs_hook for json.load, which silently keeps the last duplicated key."""
    keys = [key for key, _ in pairs]
    duplicated = sorted({k for k in keys if keys.count(k) > 1})

    if duplicated:
        filter_logger.error(f"Duplicated group names in the counts design file: {duplicated}")
        sys.exit(1)
    return dict(pairs)


def read_counts_design(design_file: str) -> dict:
    """Read the experimental design: which samples belong to each group.

    The file is a JSON object with one key per experimental group and, as
    value, the list of its samples, named as in the header of the --fl_count
    file used in QC (without the "FL." prefix of the classification columns):

        {"K": ["K1", "K2", "K3"], "B": ["B1", "B2", "B3"]}

    Args:
        design_file (str): path to the JSON design file

    Returns:
        dict: group name -> list of sample names, in file order

    Exits:
        Calls sys.exit(1) if the file is not valid JSON, is not an object of
        non-empty lists of sample names, repeats a group, or assigns a sample
        to more than one group.
    """
    try:
        with open(design_file, 'r') as f:
            design = json.load(f, object_pairs_hook=_reject_duplicate_groups)
    except json.JSONDecodeError as e:
        filter_logger.error(f"Counts design file {design_file} is not valid JSON: {e}")
        sys.exit(1)

    if not isinstance(design, dict) or not design:
        filter_logger.error("The counts design file must be a JSON object with one entry per group, "
                            'e.g. {"K": ["K1", "K2"], "B": ["B1", "B2"]}')
        sys.exit(1)

    seen = {}
    for group, samples in design.items():
        if not group.strip():
            filter_logger.error("Empty group name in the counts design file.")
            sys.exit(1)
        if (not isinstance(samples, list) or not samples
                or not all(isinstance(s, str) and s for s in samples)):
            filter_logger.error(f"Group {group!r} in the counts design file must be a non-empty "
                                "list of sample names.")
            sys.exit(1)
        for sample in samples:
            if sample in seen:
                filter_logger.error(f"Sample {sample!r} is assigned more than once in the counts design "
                                    f"file (groups {seen[sample]!r} and {group!r}).")
                sys.exit(1)
            seen[sample] = group

    return design


def check_counts_design(classif: pd.DataFrame, design: dict) -> None:
    """Check the design against the per-sample columns of the classification.

    Every sample of the design must have its FL.<sample> column. Samples of
    the classification left out of the design are allowed: they may be samples
    that should not take part in the filter (for instance, mixtures used as a
    response variable). A warning lists them.

    Exits:
        Calls sys.exit(1) if the classification has no per-sample counts, or if
        a sample of the design has no FL.<sample> column.
    """
    available = [col[3:] for col in classif.columns if col.startswith("FL.")]
    if not available:
        filter_logger.error("A counts design file was given, but the classification file has no "
                            "per-sample counts (FL.<sample> columns).")
        filter_logger.error("Run SQANTI3 QC with a multi-sample --fl_count file to generate them.")
        sys.exit(1)

    in_design = [s for samples in design.values() for s in samples]
    missing = [s for s in in_design if s not in available]
    if missing:
        filter_logger.error(f"Samples in the counts design file not found in the classification file: {missing}")
        filter_logger.error(f"Available samples: {available}")
        sys.exit(1)

    unused = [s for s in available if s not in in_design]
    if unused:
        filter_logger.warning("Samples not assigned to any group in the counts design file, "
                              f"ignored by min_prevalence rules: {unused}")

    for group, samples in design.items():
        if len(samples) == 1:
            filter_logger.warning(f"Group {group!r} has a single sample: its prevalence can only be 0 or 1.")


def resolve_prevalence_rules(rules_dict: dict, design: Optional[dict]) -> dict:
    """Map each min_prevalence value to the prevalence columns it is checked on.

    Without design, integer values are left as they are (checked on the QC
    prevalence column) and per-group values are an error. With design, every
    value becomes {prevalence.<group>: threshold}:

        2                 -> {"prevalence.K": 2, "prevalence.B": 2}
        {"K": 2, "B": 3}  -> {"prevalence.K": 2, "prevalence.B": 3}

    Per-group values must name exactly the groups of the design: a group left
    out would never retain anything, which is what removing it from the design
    already expresses. No threshold can exceed the size of its group.

    Returns:
        dict: a copy of rules_dict with min_prevalence values resolved

    Exits:
        Calls sys.exit(1) on per-group values without design, group names that
        do not match the design, or thresholds larger than their group.
    """


    def resolve(value, sc):
        if design is None:
            if isinstance(value, dict):
                filter_logger.error(f"Per-group min_prevalence at {sc} requires a counts design file "
                                    "(--counts_design).")
                sys.exit(1)
            return value

        if isinstance(value, dict):
            unknown = [g for g in value if g not in design]
            missing = [g for g in design if g not in value]
            if unknown or missing:
                if unknown:
                    filter_logger.error(f"min_prevalence at {sc} uses groups not in the counts design file: {unknown}")
                if missing:
                    filter_logger.error(f"min_prevalence at {sc} has no threshold for groups: {missing}")
                filter_logger.error(f"Groups in the counts design file: {list(design)}")
                sys.exit(1)
        thresholds = {g: value[g] if isinstance(value, dict) else value for g in design}

        for g, m in thresholds.items():
            if m > len(design[g]):
                filter_logger.error(f"min_prevalence {m} at {sc} is larger than the number of samples "
                                    f"in group {g!r} ({len(design[g])}).")
                sys.exit(1)
        return {f"prevalence.{g}": m for g, m in thresholds.items()}

    resolved = {}
    for sc, rule_sets in rules_dict.items():
        resolved[sc] = []
        for rules in rule_sets:
            rules = rules.copy()
            rules["rule"] = rules["rule"].astype(object)  # an int column cannot hold dicts
            is_prev = rules["type"] == "min_prevalence"
            rules.loc[is_prev, "rule"] = pd.Series(
                [resolve(v, sc) for v in rules.loc[is_prev, "rule"]],
                index=rules.index[is_prev], dtype=object)
            resolved[sc].append(rules)
    return resolved


def add_group_prevalence(classif: pd.DataFrame, design: dict) -> pd.DataFrame:
    """Add one prevalence.<group> column per group of the design.

    Each column counts the samples of the group where the transcript reaches
    the detection threshold, the same rule QC uses for the prevalence column.
    Missing counts are not detections.

    Returns:
        pd.DataFrame: classif with the new columns (modified in place)
    """
    for group, samples in design.items():
        counts = classif[[f"FL.{s}" for s in samples]]
        classif[f"prevalence.{group}"] = (counts >= MIN_DETECTION_COUNT).sum(axis=1)

    return classif


def prepare_prevalence_rules(classif, rules_dict, counts_design, ignore_prevalence, logger):
    """Validate and resolve min_prevalence rules before filtering.

    Returns:
        tuple: (classif, rules_dict). With a counts design, classif gains the
        prevalence.<group> columns and min_prevalence values are resolved to them.
    """
    if not has_prevalence_rules(rules_dict):
        if counts_design is not None:
            logger.warning("--counts_design is ignored: the rules file has no min_prevalence rules.")
        return classif, rules_dict

    if ignore_prevalence:
        logger.info("Ignoring min_prevalence rules (--ignore_prevalence).")
        return classif, drop_prevalence_rules(rules_dict)

    if counts_design is None:
        rules_dict = resolve_prevalence_rules(rules_dict, None)
        check_prevalence_column(classif, get_highest_min_prevalence(rules_dict))
        return classif, rules_dict

    design = read_counts_design(counts_design)
    check_counts_design(classif, design)
    rules_dict = resolve_prevalence_rules(rules_dict, design)
    classif = add_group_prevalence(classif, design)
    logger.info("Counts design: " + ", ".join(f"{g} ({len(s)} samples)" for g, s in design.items()))
    return classif, rules_dict
        


def rules_filter(sqanti_class, json_file, force_multi_exon, prefix, logger, ignore_prevalence=False, counts_design=None):
    """Main function to execute SQANTI3 filtering workflow.
    
    Args:
        sqanti_class (str): Path to SQANTI3 classification file (TSV format)
        json_file (str): Path to JSON file containing filtering rules
        force_multi_exon (bool): If True, exclude all mono-exonic transcripts
        prefix (str): Output filename prefix
        logger (logging.Logger): Configured logger for progress reporting
        ignore_prevalence (bool): If True, min_prevalence rules are removed before filtering.
            SQANTI3 rescue sets it to filter the reference transcriptome.
        counts_design (str): Path to a JSON file assigning samples to experimental groups.
        With it, min_prevalence is evaluated per group and an isoform passes if it
    r   eaches the threshold in any group. Without it, over all samples together.
        
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
    classif, rules_dict = prepare_prevalence_rules(classif, rules_dict, counts_design,
                                                   ignore_prevalence, logger)

    message("Applying rules to filter isoforms",logger)

    classif['filter_result'] = classif.apply(lambda row: apply_rules(row, force_multi_exon, rules_dict), axis=1)

    inclusion_list = classif[classif['filter_result'] == "Isoform"]['isoform']

    artifacts_classif = classif[classif['filter_result'] == "Artifact"]
    reasons_df = artifacts_classif.apply(lambda row: get_reasons(row, force_multi_exon, rules_dict), axis=1)

    message("Writing results",logger)
    classif.to_csv(os.path.join(f"{prefix}_RulesFilter_classification.txt"), sep='\t', index=False)
    inclusion_list.to_csv(os.path.join(f"{prefix}_pass_isoforms.txt"), sep='\t', index=False, header=False)
    reasons_df.to_csv(os.path.join(f"{prefix}_filtering_reasons.txt"), sep='\t', index=False)

