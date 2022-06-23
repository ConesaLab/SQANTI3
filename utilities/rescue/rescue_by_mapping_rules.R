#!/bin/bash Rscript
#
####################################################
#########     SQANTI3 RESCUE BY MAPPING    #########
####################################################
#
# Author: Ángeles Arzalluz-Luque
# Contact: angeles.arzalluz@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# DESCRIPTION:
# This script contains the code for the third module of the SQANTI3 rescue
# pipeline, which integrates rescue candidate alignment results (to their same-gene
# counterparts, both from the reference and from the long read defined transcriptome). 
# Using the rules filter results obtained for all rescue targets (both reference 
# and long read-defined transcripts), it selects the best rescue target from 
# those that were retrieved during alignment.
#