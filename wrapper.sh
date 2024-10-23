#! /bin/bash -xe

# Author: Fabián Robledo
# Email: fabian.robledo@csic.es
# Version: 0.9.0

# One SQANTI3 to rule them all
# A wrapper, with a config file, to execute al three steps of a sqanti3 pipeline
# in one execution
# Al three steps (QC, filter and rescue) can be executed at will , as long
# as the conf is told to skip the script in question

function main () {
    
    # We check every optional argument to treat them properly
    # If they are empty (""), then they won't be needed downstream
    # Because bash ignores all empty variables it finds. 
    # Mandatory arguments are not needed to check because they don't use
    # a flag to indicate it, so it can be used as it is in the conf
    
    # SQANTI3 QC parameters
        if [ -n "${cage_peak_bed_file}" ]; then cage_peak_bed_file="--CAGE_peak ${cage_peak_bed_file}"; fi;
        if [ -n "${polyA_motif_list}"  ]; then polyA_motif_list="--polyA_motif_list $polyA_motif_list"; fi; 
        if [ -n "${output_prefix}" ]; then output_prefix="-o  ${output_prefix}";fi;
        if [ -n "${destination_folder_QC}" ]; then mkdir -p ${destination_folder_QC}; destination_folder_QC="-d ${destination_folder_QC}";fi;
        if [ -n "${full_length_pacbio_abundance_tsv}" ]; then full_length_pacbio_abundance_tsv="-fl  ${full_length_pacbio_abundance_tsv}";fi;
        if [ -n "${short_reads_fofn}" ]; then short_reads_fofn="--short_reads ${short_reads_fofn}"; fi;
        if [ -n "${cpus}" ]; then cpus="--cpus $cpus"; else cpus="--cpus 1"; fi;
        if [ -n "${report_file}" ]; then report_file="--report ${report_file}"; fi;
        if [ -n "${chunks}" ]; then chunks="-n ${chunks}"; fi;
        if [ -n "${reference_fasta}"  ]; then reference_fasta_rescue="-f ${reference_fasta}";fi
        
    # SQANTI3 Filter common parameters
        if [ -n "${gtf_to_filter}"  ]; then gtf_to_filter="--gtf ${gtf_to_filter}";fi
        if [ -n "${isoforms}"  ]; then isoforms="--isoforms ${isoforms}";fi
        if [ -n "${isoannotgff3}"  ]; then isoannotgff3="--isoAnnotGFF3 ${isoannotgff3}";fi
        if [ -n "${sam}"  ]; then sam="--sam ${sam}";fi
        if [ -n "${faa}"  ]; then faa="--faa ${faa}";fi
        if [ -n "${monoexonic}"  ]; then monoexonic="-e";fi
        if [ -n "${skip_report}" ] ; then skip_report="--skip_report";fi
        
    # Filter rules parameters
        if [ -n "${ouput_folder_filter_rules}" ]; then mkdir -p ${ouput_folder_filter_rules}; ouput_folder_filter_rules="-d ${ouput_folder_filter_rules}";fi;
        if [ -n "${prefix_filter_rules}"  ]; then prefix_filter_rules="-o ${prefix_filter_rules}";fi
        if [ -n "${json_file}"  ]; then json_file="-j ${json_file}";fi
        
    # Filter ml parameters
        if [ -n "${ouput_folder_filter_ml}" ]; then mkdir -p ${ouput_folder_filter_ml}; ouput_folder_filter_ml="-d ${ouput_folder_filter_ml}";fi;
        echo $ouput_folder_filter_ml
        if [ -n "${prefix_filter_ml}"  ]; then prefix_filter_ml="-o ${prefix_filter_ml}";fi
        if [ -n "${forest}"  ]; then forest="-r ${forest}";fi
        if [ -n "${threshold}"  ]; then threshold="-j ${threshold}";fi
        if [ -n "${percent_training}"  ]; then percent_training="-t ${percent_training}";fi
        
        
     # Sqanti3 rescue rules parameters
        if [ -n "${rescue_rules_reference_classification}"  ]; then rescue_rules_reference_classification="-k ${rescue_rules_reference_classification}";fi
        if [ -n "${rescue_rules_isoforms}"  ]; then rescue_rules_isoforms="--isoforms ${rescue_rules_isoforms}";fi
        if [ -n "${rescue_rules_gtf}"  ]; then rescue_rules_gtf="--gtf ${rescue_rules_gtf}";fi
        if [ -n "${rescue_rules_reference_gtf}"  ]; then rescue_rules_reference_gtf="-g ${rescue_rules_reference_gtf}";fi
        if [ -n "${rescue_rules_monoexons}"  ]; then rescue_rules_monoexons="-e ${rescue_rules_monoexons}";fi
        if [ -n "${rescue_rules_mode}"  ]; then rescue_rules_mode="--mode ${rescue_rules_mode}";fi
        if [ -n "${rescue_rules_output_folder}"  ]; then rescue_rules_output_folder="-d ${rescue_rules_output_folder}";fi
        if [ -n "${rescue_rules_output_prefix}"  ]; then rescue_rules_output_prefix="-o ${rescue_rules_output_prefix}";fi
        if [ -n "${rescue_rules_json_file}"  ]; then rescue_rules_json_file="-j ${rescue_rules_json_file}";fi
        if [ -n "${rescue_rules_reference_genome}"  ]; then rescue_rules_reference_genome="-f ${rescue_rules_reference_genome}";fi
        
    # Sqanti3 Rescue ml parameters
        if [ -n "${rescue_ml_reference_classification}"  ]; then rescue_ml_reference_classification="-k ${rescue_ml_reference_classification}";fi
        if [ -n "${rescue_ml_isoforms}"  ]; then rescue_ml_isoforms="--isoforms ${rescue_ml_isoforms}";fi
        if [ -n "${rescue_ml_gtf}"  ]; then rescue_ml_gtf="--gtf ${rescue_ml_gtf}";fi
        if [ -n "${rescue_ml_reference_gtf}"  ]; then rescue_ml_reference_gtf="-g ${rescue_ml_reference_gtf}";fi
        if [ -n "${rescue_ml_monoexons}"  ]; then rescue_ml_monoexons="-e ${rescue_ml_monoexons}";fi
        if [ -n "${rescue_ml_mode}"  ]; then rescue_ml_mode="--mode ${rescue_ml_mode}";fi
        if [ -n "${rescue_ml_output_folder}"  ]; then rescue_ml_output_folder="-d ${rescue_ml_output_folder}";fi
        if [ -n "${rescue_ml_output_prefix}"  ]; then rescue_ml_output_prefix="-o ${rescue_ml_output_prefix}";fi
        if [ -n "${rescue_ml_reference_genome}"  ]; then rescue_ml_reference_genome="-f ${rescue_ml_reference_genome}";fi
        if [ -n "${rescue_ml_threshold}"  ]; then rescue_ml_threshold="-j ${rescue_ml_threshold}";fi
        if [ -n "${rescue_ml_randomforest_rdata}"  ]; then rescue_ml_randomforest_rdata="-r ${rescue_ml_randomforest_rdata}";fi
    
    # Check if skip_qc exist and/or is false. Possible outcomes
    # skip_qc is not declared, is "" or is false -> QC is performed
    # skip_qc has any other different value -> QC is skipped
    if [ -z "${skip_qc}" ] || [ ${skip_qc} == "false" ]; 
    then
        ${sqanti3_qc} ${input_gtf} ${reference_gtf}  ${reference_fasta}   \
                             ${cage_peak_bed_file}    \
                             ${polyA_motif_list}   \
                             ${output_prefix}  \
                             ${destination_folder_QC} \
                             ${full_length_pacbio_abundance_tsv}   \
                             ${short_reads_fofn} \
                             ${cpus}  \
                             ${report_file} ${chunks}
    else
        echo "Skipping sqanti3_qc"
    fi
    
    # Check if skip_filter exist and/or is false. Possible outcomes
    # skip_filter is not declared, is "" or is false -> Filter is performed
    # skip_filter has any other different value -> Filter is skipped
    if [ -z "${skip_filter}" ] || [ ${skip_filter} == "false" ]; 
    then
            # If filtering is performed, mode is checked
            # There are 4 outcomes according to filter_mode value:
            # "both" or empty string -> rules and ml are performed
            # "rules" -> Only rules is performed
            # "ml" -> Only ml is performed
            # Other values -> Nothing is performed. This may be a source of bugs
            # TODO: check only that valid values are available
            if  [ -z ${filter_mode} ] ||  [ ${filter_mode} == "rules" ] || [ ${filter_mode} == "both" ];
            then
                
                
                ${sqanti3_filter} rules ${ouput_folder_filter_rules} \
                    ${gtf_to_filter} ${prefix_filter_rules} ${isoforms} ${isoannotgff3} \
                    ${sam} ${faa} ${monoexonic} ${skip_report} ${json_file} \
                    ${classification_file_to_filter}
            fi
            
            if [ -z ${filter_mode} ] || [ ${filter_mode} == "ml" ] || [ $filter_mode == "both" ];
            then 
                ${sqanti3_filter} ml ${ouput_folder_filter_ml} \
                ${gtf_to_filter} ${prefix_filter_rules} ${isoforms} ${isoannotgff3} \
                ${sam} ${faa} ${monoexonic} ${skip_report} ${percent_training} ${TP} \
                ${TN} ${threshold} ${forest} ${classification_file_to_filter}
            fi
            
    else
        echo "Skipping sqanti3_filter"
    fi
    
    
    # Check if skip_rescue exist and/or is false. Possible outcomes
    # skip_filter is not declared, is "" or is false -> Filter is performed
    # skip_filter has any other different value -> Filter is skipped
    if [ -z "${skip_rescue}" ] || [ ${skip_rescue} == "false" ]; 
    then
        
        # If rescue is performed, mode is checked
        # There are 4 outcomes according to rescue_mode value:
        # "both" or empty string -> rules and ml are performed
        # "rules" -> Only rules is performed
        # "ml" -> Only ml is performed
        # Other values -> Nothing is performed. This may be a source of bugs
        # TODO: check only that valid values are available
        if [ ${rescue_mode} == "rules" ] || [ ${rescue_mode} == "both" ];
        then
            ${sqanti3_rescue} rules ${rescue_rules_reference_genome}  ${rescue_rules_isoforms} ${rescue_rules_gtf} ${rescue_rules_ref_classif}  \
                ${rescue_rules_reference_gtf} ${rescue_rules_reference_classification} ${rescue_rules_mode} ${rescue_rules_monoexons} \
                ${rescue_rules_output_prefix} ${rescue_rules_output_folder} ${rescue_rules_json_file} \
                ${rescue_rules_filtered_classification}
        fi
        if [ $rescue_mode == "ml" ] || [ $rescue_mode == "both" ];
        then 
            ${sqanti3_rescue} ml ${rescue_ml_isoforms} ${rescue_ml_gtf} ${rescue_ml_reference_classification}  \
                ${rescue_ml_reference_gtf} ${rescue_ml_reference_gtf} ${rescue_ml_mode} \
                ${rescue_ml_threshold} ${rescue_ml_monoexons} ${rescue_ml_reference_genome} \
                ${rescue_ml_output_prefix} ${rescue_ml_output_folder} ${rescue_ml_randomforest_rdata} \
                ${rescue_ml_filtered_classification} 
        fi
    else
        echo "Skipping sqanti3_rescue"
    fi
    
}

# The config file: wrapper.conf
# A list of name=value parameters used in this script
# Check the conf for more details on each parameters
source "$1"

# Execuing the main function
main