import os

import pandas as pd
import numpy as np
from src.module_logging import rescue_logger, message
from src.commands import (
    RSCRIPTPATH, RESCUE_AUTO_PATH, utilitiesPath,run_command 
)
from src.utilities.rescue.automatic_rescue import (
    read_classification, rescue_fsm_monoexons, add_ism_monoexons,
    get_lost_reference_id, 
    rescue_lost_reference, save_automatic_rescue
)

def run_automatic_rescue(args):

    auto_cmd = f"{RSCRIPTPATH} {RESCUE_AUTO_PATH} -c {args.filter_class} -o {args.output} -d {args.dir} \
    -u {utilitiesPath} -g {args.refGTF} -e {args.rescue_mono_exonic} -m {args.mode}"
    logFile = os.path.join(args.dir,"logs","automatic_rescue.log")
    run_command(auto_cmd, rescue_logger,logFile,description="Running automatic rescue",silent=False)
    # print command
    rescue_logger.debug("Automatic rescue run via the following command:")
    rescue_logger.debug(auto_cmd)

    ## load output: transcripts rescued as a result of automatic rescue

    # make file name
    automatic_rescued_list = os.path.join(args.dir,f"{args.output}_automatic_rescued_list.tsv")

    # set object containing rescued list from the output file
    auto_rescue = set(line.strip() for line in open(automatic_rescued_list))

    ## return automatic rescue outputs
    return(auto_rescue)

def run_automatic_rescue_py(classification_file,monoexons,mode,prefix):
    # Load classification
    message("Reading filter classification file",rescue_logger)
    classif_df = read_classification(classification_file)

    message("Performing automatic rescue",rescue_logger)
    # Select the FSM and ISM isoforms with more than one exon
    ism_fsm_classif = classif_df[
        (classif_df['structural_category'].isin(['full-splice_match','incomplete-splice_match'])) & 
        (classif_df['exons'] > 1)
    ]
    # Add the monoexons if indicated
    if monoexons in ['all','fsm']:
        rescue_fsm_me = rescue_fsm_monoexons(classif_df)
        if mode == 'full':
            ism_fsm_classif = add_ism_monoexons(ism_fsm_classif,classif_df)


    # Find the references that are lost and get the ones that are not represented by isoforms
    lost_ref = get_lost_reference_id(ism_fsm_classif)
    rescue_logger.debug(f"Found {len(lost_ref)} lost references")
    rescue = pd.DataFrame()
    for ref_id in lost_ref:
        rescue_df = rescue_lost_reference(ref_id, ism_fsm_classif)
        rescue = pd.concat([rescue,rescue_df])

    # Split into reference transcripts and ISM
    rescue_ism = rescue[rescue['isoform'].isin(ism_fsm_classif['isoform'])]
    rescue_ref = rescue[rescue['isoform'].isin(ism_fsm_classif['associated_transcript'])]

    if monoexons in ['all','fsm']:
        rescue_auto = pd.concat([rescue_ref,rescue_fsm_me])
    else:
        rescue_auto = rescue_ref
    rescue_logger.debug(f"Rescued {rescue_auto.shape[0]} transcripts")

    # Save the automatic rescue
    save_automatic_rescue(rescue_auto,ism_fsm_classif,mode,prefix)
    message("Automatic rescue completed",rescue_logger)
