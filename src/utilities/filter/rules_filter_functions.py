import pandas as pd

def apply_rules(isoform_info, force_multiexon, json_df, rules_list):
    if force_multiexon and float(isoform_info["exons"]) == 1:
        final_is_isoform = False
    else:
        sc = str(isoform_info["structural_category"])
        final_is_isoform = True
        
        if sc in json_df.keys():
            final_is_isoform = False
            for p in [i for i, name in enumerate(rules_list) if name == sc]:
                is_isoform = True
                rules = rules_list[p]
                for i in range(len(rules['rule'])):
                    parameter = rules.loc[i, "column"]
                    if pd.notna(isoform_info[parameter]):
                        if rules.loc[i, "type"] == "Min_Threshold":
                            if float(isoform_info[parameter]) < float(rules.loc[i, "rule"]):
                                is_isoform = False
                                break
                        elif rules.loc[i, "type"] == "Max_Threshold":
                            if float(isoform_info[parameter]) > float(rules.loc[i, "rule"]):
                                is_isoform = False
                                break
                        elif rules.loc[i, "type"] == "Category":
                            cat_rules = rules[rules['column'] == parameter]
                            if isoform_info[parameter].lower() not in cat_rules["rule"].tolist():
                                is_isoform = False
                                break
                    else:
                        is_isoform = False
                        break
                final_is_isoform = final_is_isoform or is_isoform
        elif "rest" in json_df.keys():
            final_is_isoform = False
            for p in [i for i, name in enumerate(rules_list) if name == "rest"]:
                is_isoform = True
                rules = rules_list[p]
                for i in range(len(rules['rule'])):
                    parameter = rules.loc[i, "column"]
                    if pd.notna(isoform_info[parameter]):
                        if rules.loc[i, "type"] == "Min_Threshold":
                            if float(isoform_info[parameter]) < float(rules.loc[i, "rule"]):
                                is_isoform = False
                                break
                        elif rules.loc[i, "type"] == "Max_Threshold":
                            if float(isoform_info[parameter]) > float(rules.loc[i, "rule"]):
                                is_isoform = False
                                break
                        elif rules.loc[i, "type"] == "Category":
                            cat_rules = rules[rules['column'] == parameter]
                            if isoform_info[parameter].lower() not in cat_rules["rule"].tolist():
                                is_isoform = False
                                break
                    else:
                        is_isoform = False
                        break
                final_is_isoform = final_is_isoform or is_isoform
    
    return "Isoform" if final_is_isoform else "Artifact"

def get_reasons(isoform_info, force_multiexon, json_df, rules_list):
    if force_multiexon and float(isoform_info["exons"]) == 1:
        reasons = ["Mono-exon"]
    else:
        sc = str(isoform_info["structural_category"])
        reasons = []
        
        if sc in json_df.keys():
            for p in [i for i, name in enumerate(rules_list) if name == sc]:
                rules = rules_list[p]
                for i in range(len(rules['rule'])):
                    parameter = rules.loc[i, "column"]
                    if pd.notna(isoform_info[parameter]):
                        if rules.loc[i, "type"] == "Min_Threshold":
                            if float(isoform_info[parameter]) < float(rules.loc[i, "rule"]):
                                reasons.append(f"Low {parameter}")
                        elif rules.loc[i, "type"] == "Max_Threshold":
                            if float(isoform_info[parameter]) > float(rules.loc[i, "rule"]):
                                reasons.append(f"High {parameter}")
                        elif rules.loc[i, "type"] == "Category":
                            cat_rules = rules[rules['column'] == parameter]
                            if isoform_info[parameter].lower() not in cat_rules["rule"].tolist():
                                reasons.append(f"Out {parameter}")
                    else:
                        reasons.append(f"NA value in {parameter}")
        elif "rest" in json_df.keys():
            for p in [i for i, name in enumerate(rules_list) if name == "rest"]:
                rules = rules_list[p]
                for i in range(len(rules['rule'])):
                    parameter = rules.loc[i, "column"]
                    if pd.notna(isoform_info[parameter]):
                        if rules.loc[i, "type"] == "Min_Threshold":
                            if float(isoform_info[parameter]) < float(rules.loc[i, "rule"]):
                                reasons.append(f"Low {parameter}")
                        elif rules.loc[i, "type"] == "Max_Threshold":
                            if float(isoform_info[parameter]) > float(rules.loc[i, "rule"]):
                                reasons.append(f"High {parameter}")
                        elif rules.loc[i, "type"] == "Category":
                            cat_rules = rules[rules['column'] == parameter]
                            if isoform_info[parameter].lower() not in cat_rules["rule"].tolist():
                                reasons.append(f"Out {parameter}")
                    else:
                        reasons.append(f"NA value in {parameter}")
    
    unique_reasons = list(set(reasons))
    final_df = pd.DataFrame({
        'isoform': [isoform_info["isoform"]] * len(unique_reasons),
        'structural_category': [isoform_info["structural_category"]] * len(unique_reasons),
        'reasons': unique_reasons
    })
    return final_df
