#!/bin/bash
# This script has to be run from the main directory of SQANTI3

# Run SQANTIQC
# sqanti3 qc -c example/config_files/qc_config.yaml
# sqanti3 qc -c example/config_files/qc_config_reference.yaml

# # Run SQANTI filter
# sqanti3 filter -c example/config_files/filter_rules.yaml
# sqanti3 filter -c example/config_files/filter_ML.yaml

# Run SQANTI rescue
# sqanti3 rescue -c example/config_files/rescue_automatic.yaml
# sqanti3 rescue -c example/config_files/rescue_rules.yaml
sqanti3 rescue -c example/config_files/rescue_ml.yaml