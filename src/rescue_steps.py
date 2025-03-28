import os
import subprocess
import sys

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
from src.utilities.rescue.rescue_helpers import (
    get_rescue_gene_targets,get_rescue_reference_targets
)

from src.utilities.rescue.candidate_mapping_helpers import (
    filter_transcriptome,
    prepare_fasta_transcriptome
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
    return rescue_ism

def rescue_candidates(classification_file,monoexons,rescue_ism,prefix):
    # Load classification
    classif_df = read_classification(classification_file)
    
    # Get NIC and NNC artifacts
    rescue_novel = classif_df[
        (classif_df['structural_category'].isin(['novel_in_catalog','novel_not_in_catalog'])) &
        (classif_df['filter_result'] == 'Artifact')
    ]
    if monoexons != 'all':
        rescue_novel = rescue_novel[rescue_novel['exons'] > 1]
    
    rescue_candidates = pd.concat([rescue_ism,rescue_novel['isoform']])
    # Write rescue candidates
    rescue_candidates.to_csv(f"{prefix}_rescue_candidates.tsv", 
                            sep="\t",
                            index=False,
                            header=False)

    return rescue_candidates

def rescue_targets(classification_file,rescue_candidates,ref_gtf,prefix):
    # Load classification
    classif_df = read_classification(classification_file)
    
    # Get the genes with associated rescue candidates
    target_genes = get_rescue_gene_targets(classif_df,rescue_candidates)

    rescue_targets_lr = classif_df[
        (classif_df['associated_gene'].isin(target_genes)) &
        (classif_df['filter_result'] == 'Isoform')
    ]['isoform']

    rescue_targets_ref = get_rescue_reference_targets(ref_gtf,target_genes)
    rescue_logger.debug(f"Found {rescue_targets_ref.shape[0]} reference targets")
    rescue_logger.debug(f"Found {rescue_targets_lr.shape[0]} long-read targets")
    # Merge both groups and remove duplicates
    rescue_targets = pd.concat([rescue_targets_lr, rescue_targets_ref]).drop_duplicates().reset_index(drop=True)  
    rescue_logger.debug(f"Rescue targets: {rescue_targets.shape[0]}")
    
    rescue_targets.to_csv(f"{prefix}_rescue_targets.tsv",
                        sep="\t",
                        index=False,
                        header=False)
    return rescue_targets

## Run mapping of rescue candidates (artifacts) to targets
def run_candidate_mapping(args,targets_list,candidates_list):
    rescue_logger.debug(targets_list)
    #### PREPARATION OF FILES FOR MINIMAP2 ####
    rescue_logger.info("**** Preparation of files for artifact mapping:")
    ## Convert reference transcriptome GTF to FASTA

    ref_trans_fasta = prepare_fasta_transcriptome(args.refGTF,args.refFasta,args.dir)

    ## Filter reference transcriptome FASTA to only include target ref transcripts
    rescue_logger.info("Filtering reference transcriptome FASTA to only rescue targets.")

    # make file names
    target_file = f"{args.dir}/{args.output}_rescue_targets.tsv"
    ref_target_fasta = f"{args.dir}/{args.output}_rescue_targets.ref.fasta"

    filter_transcriptome(ref_trans_fasta,targets_list,ref_target_fasta)
    input()
# # make command
# fasta_cmd = f"seqtk subseq {ref_trans_fasta} {target_file} > {ref_target_fasta}"

# # run
# try:
#     subprocess.check_call(fasta_cmd, shell=True)
#     if os.path.isfile(ref_target_fasta):
#         rescue_logger.info(f"Target reference transcript sequences were saved to {ref_target_fasta}")
#         rescue_logger.info("seqtk command used:")
#         rescue_logger.info(fasta_cmd)
#     else:
#         rescue_logger.error("Target reference transcript FASTA was not created - file not found!")
#         sys.exit(1)
# except subprocess.CalledProcessError:
#     rescue_logger.error(f"Error retrieving target reference transcripts from FASTA. Command used: {fasta_cmd}")
#     sys.exit(1)


    ## Filter SQ3 transcriptome FASTA to only include target LR transcripts

    rescue_logger.info("Filtering supplied long read transcriptome FASTA (--isoforms) to only include rescue targets...")

    # make file names
    LR_target_fasta = f"{args.dir}/{args.output}_rescue_targets.LR.fasta"
    filter_transcriptome(args.rescue_isoforms,targets_list,LR_target_fasta)
    # make command
# fasta_cmd = f"seqtk subseq {args.rescue_isoforms} {target_file} > {LR_target_fasta}"

# # run
# try:
#     subprocess.check_call(fasta_cmd, shell=True)
#     if os.path.isfile(LR_target_fasta):
#         rescue_logger.info(f"Target long read transcript sequences were saved to {LR_target_fasta}")
#         rescue_logger.info("seqtk command used:")
#         rescue_logger.info(fasta_cmd)
#     else:
#         rescue_logger.error("Target long read transcript FASTA was not created - file not found!")
#         sys.exit(1)
# except subprocess.CalledProcessError:
#     rescue_logger.error(f"Error retrieving target long-read transcripts from FASTA. Command used: {fasta_cmd}")
#     sys.exit(1)


    ## join both FASTA files
    rescue_logger.info("Joining reference and LR rescue target FASTA files...")
    # TODO: Find a python way to cat both files together. Perhaps not saving the files and returning the object?
    target_fasta = f"{args.dir}/{args.output}_rescue_targets.fasta"
    cat_cmd = f"cat {ref_target_fasta} {LR_target_fasta} > {target_fasta}"

    try:
        subprocess.check_call(cat_cmd, shell=True)
        if os.path.isfile(target_fasta):
            rescue_logger.info(f"Rescue target FASTA was saved to {target_fasta}")
            rescue_logger.info("Command used:")
            rescue_logger.info(cat_cmd)

            # Remove intermediate target FASTA files (LR and ref)
            rescue_logger.info("Removing intermediate target FASTA files...")
            rm_cmd = f"rm {ref_target_fasta} {LR_target_fasta}"
            try:
                subprocess.call(rm_cmd, shell=True)
            except subprocess.CalledProcessError:
                rescue_logger.warning(f"Failed to remove intermediate files. Command used: {rm_cmd}")
        else:
            rescue_logger.error("Target FASTA was not created - file not found!")
            sys.exit(1)
    except subprocess.CalledProcessError:
        rescue_logger.error(f"Error joining target long-read and reference FASTA files. Command used: {cat_cmd}")
        sys.exit(1)


    ## Filter SQ3 FASTA to include rescue candidates
    rescue_logger.info("Creating rescue candidate FASTA from supplied long read transcriptome fasta (--isoforms)...")

    # make file names
    candidate_file = f"{args.dir}/{args.output}_rescue_candidates.tsv"
    candidate_fasta = f"{args.dir}/{args.output}_rescue_candidates.fasta"
    filter_transcriptome(args.rescue_isoforms,candidates_list,candidate_fasta)
# # make command
# fasta_cmd = f"seqtk subseq {args.rescue_isoforms} {candidate_file} > {candidate_fasta}"

# # run
# try:
#     subprocess.check_call(fasta_cmd, shell=True)
#     if os.path.isfile(candidate_fasta):
#         rescue_logger.info(f"Rescue candidate FASTA was saved to {candidate_fasta}")
#         rescue_logger.info("seqtk command used:")
#         rescue_logger.info(fasta_cmd)
#     else:
#         rescue_logger.error("Candidate FASTA was not created - file not found!")
#         sys.exit(1)
# except subprocess.CalledProcessError:
#     rescue_logger.error(f"Error retrieving rescue candidate sequences from FASTA. Command used: {fasta_cmd}")
#     sys.exit(1)


    #### MAPPING ARTIFACTS (CANDIDATES) WITH MINIMAP2 ####
    # TODO: eliminate file logic in the process
    rescue_logger.info("**** Artifact mapping (candidates vs targets)")
    
    # Mapping
    rescue_logger.info("Mapping rescue candidates to rescue targets with minimap2...")

    # make file names
    sam_file = f"{args.dir}/{args.output}_mapped_rescue.sam"

    # make command
    minimap_cmd = f"minimap2 --secondary=yes -ax map-hifi {target_fasta} {candidate_fasta} > {sam_file}"

    # run
    logFile=f"{args.dir}/logs/rescue/minimap2.log"
    run_command(minimap_cmd,rescue_logger,logFile,"Maping rescue candidates to targets")

    if os.path.isfile(sam_file):
        rescue_logger.info(f"Minimap2 results were saved to {sam_file}")
        rescue_logger.info("minimap2 command used:")
        rescue_logger.info(minimap_cmd)

    # Filter mapping results (select SAM columns)
    rescue_logger.info("Building candidate-target table of mapping hits...")

    # remove header from SAM
    sam_tmp_file = f"{args.dir}/{args.output}_mapped_rescue_noheader.sam"
    sam_cmd = f"grep -v '@' {sam_file} > {sam_tmp_file}"

    run_command(sam_cmd,rescue_logger,"log/rescue/sam_noheader.log",description="Remove header from SAM file")
    if os.path.isfile(sam_tmp_file):
      # get cols with candidate-target pairs + alignment type
      hits_file = f"{args.dir}/{args.output}_rescue_mapping_hits.tsv"
      hits_cmd = f"cut -f1-3 {sam_tmp_file} > {hits_file}"

      run_command(hits_cmd,rescue_logger,"log/rescue/hits.log",description="Extract candidate-target pairs from SAM file")

      if os.path.isfile(hits_file):
        rescue_logger.info(f"Mapping hit table was saved to {hits_file}")

        # delete altered SAM file
        os.remove(rm_cmd)