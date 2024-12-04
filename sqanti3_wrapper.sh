filter_corrected_gtf#! /bin/bash -xe

# Author: Fabián Robledo
# Email: fabian.robledo@csic.es
# Version: 1.0.0

# One SQANTI3 to run them all
# A wrapper, with a config file, to execute al three steps of a sqanti3 pipeline
# in one execution
# Al three steps (QC, filter and rescue) can be executed at will , as long
# as the conf is told to skip the script in question

# SQANTI3 (QC, Filter and Rescue) combined

function main () {
    
    # We check every optional argument to treat them properly:
    # If optional arguments have no value, they are ignored, and kept as empty strings "".
    # If not, the needed flag is added
    # If they are empty (""), then they won't be needed downstream
    # Because bash ignores all empty variables it finds. 
    # Mandatory arguments are not needed to check because they don't use
    # a flag to indicate it, so it can be used as they are in the conf
    
    # Common parameters for all three steps   
    # SQANTI3 QC parameters
        if [ -n "${QC_min_ref_length}" ]; then QC_min_ref_length="--min-ref-len ${QC_min_ref_length}"; fi;
        if [ -n "${QC_force_id_ignore}" ] &&  [ ${QC_force_id_ignore} == "true" ]; then QC_force_id_ignore="--force_id_ignore"; else QC_force_id_ignore=""; fi;
        if [ -n "${QC_cage_peak_bed_file}" ]; then QC_cage_peak_bed_file="--CAGE_peak ${QC_cage_peak_bed_file}"; fi;
        if [ -n "${QC_aligner_choice}" ]; then QC_aligner_choice="--aligner_choice ${QC_aligner_choice}"; fi;
        if [ -n "${QC_polyA_motif_list}"  ]; then QC_polyA_motif_list="--polyA_motif_list ${QC_polyA_motif_list}"; fi; 
        if [ -n "${QC_phylobed}"  ]; then QC_phylobed="--phylobed ${QC_phylobed}"; fi; 
        if [ -n "${QC_is_fusion}"  ] && [ ${QC_is_fusion} == "true"]; then QC_is_fusion="--is_fusion"; else QC_is_fusion=""; fi; 
        if [ -n "${QC_orf_input}"  ]; then QC_orf_input="--orf_input ${QC_orf_input}"; fi; 
        if [ -n "${QC_is_fastq}"  ] && [ ${QC_is_fastq} == "true" ];  then QC_is_fastq="--fasta"; else QC_is_fastq=""; fi; 
        if [ -n "${QC_expression_matrix}"  ]; then QC_expression_matrix="-e ${QC_expression_matrix}"; fi; 
        if [ -n "${QC_gmap_index}"  ]; then QC_gmap_index="--gmap_index ${QC_gmap_index}"; fi;         if [ -n "${QC_polyA_peak}"  ]; then QC_polyA_peak="--polyA_peak ${QC_polyA_peak}"; fi; 
        if [ -n "${QC_chunks}" ]; then QC_chunks="-n ${QC_chunks}"; fi;
        if [ -n "${QC_output_prefix}" ]; then QC_output_prefix="-o  ${QC_output_prefix}";fi;
        if [ -n "${QC_destination_folder}" ]; then mkdir -p ${QC_destination_folder}; QC_destination_folder="-d ${QC_destination_folder}";fi;
        if [ -n "${QC_coverage}"  ]; then QC_coverage="-c ${QC_coverage}"; fi; 
        if [ -n "${QC_sites}"  ]; then QC_sites="-s ${QC_sites}"; fi; 
        if [ -n "${QC_window}"  ]; then QC_window="-w ${QC_window}"; fi; 
        if [ -n "${QC_genename}"  ]; then QC_genename="--genename ${QC_genename}"; fi; 
        if [ -n "${QC_full_length_pacbio_abundance_tsv}" ]; then QC_full_length_pacbio_abundance_tsv="-fl  ${QC_full_length_pacbio_abundance_tsv}";fi;
        if [ -n "${QC_saturation}"  ] && [${QC_saturation} == "true"]; then QC_saturation="--saturation"; else QC_saturation="";fi; 
        if [ -n "${QC_report_file}" ]; then QC_report_file="--report ${QC_report_file}"; fi;
        if [ -n "${QC_isoAnnotLite}" ]; then QC_isoAnnotLite="--isoAnnotLite ${QC_isoAnnotLite}"; fi;
        if [ -n "${QC_gff3}" ]; then QC_gff3="--report ${QC_gff3}"; fi;
        if [ -n "${QC_short_reads_fofn}" ]; then QC_short_reads_fofn="--short_reads ${QC_short_reads_fofn}"; fi;
        if [ -n "${QC_SR_bam}" ]; then QC_SR_bam="--SR_bam ${QC_SR_bam}"; fi;
        if [ -n "${QC_isoform_hits}" ]; then QC_isoform_hits="--isoform_hits ${QC_isoform_hits}"; fi;
        if [ -n "${QC_ratio_TSS_metric}" ]; then QC_ratio_TSS_metric="--ratio_TSS_metric ${QC_ratio_TSS_metric}"; fi;
        if [ -n "${QC_cpus}" ]; then QC_cpus="--cpus $cpus"; else QC_cpus="--cpus 1"; fi;


    # SQANTI3 Filter common parameters
        if [ -n "${filter_corrected_gtf}"  ]; then filter_corrected_gtf="--gtf ${filter_corrected_gtf}";fi
        if [ -n "${filter_isoforms}"  ]; then filter_isoforms="--isoforms ${filter_isoforms}";fi
        if [ -n "${filter_isoannotgff3}"  ]; then filter_isoannotgff3="--isoAnnotGFF3 ${filter_isoannotgff3}";fi
        if [ -n "${filter_sam}"  ]; then filter_sam="--sam ${filter_sam}";fi
        if [ -n "${filter_faa}"  ]; then filter_faa="--faa ${filter_faa}";fi
        if [ -n "${filter_monoexonic}"  ]; then filter_monoexonic="-e";fi
        if [ -n "${filter_skip_report}" ] ; then filter_skip_report="--skip_report";fi
        
    # Filter rules parameters
        if [ -n "${filter_rules_ouput_folder}" ]; then mkdir -p ${filter_rules_ouput_folder}; filter_rules_ouput_folder="-d ${filter_rules_ouput_folder}";fi;
        if [ -n "${filter_rules_prefix}"  ]; then filter_rules_prefix="-o ${filter_rules_prefix}";fi
        if [ -n "${filter_rules_json_file}"  ]; then filter_rules_json_file="-j ${filter_rules_json_file}";fi
        
    # Filter ml parameters
        if [ -n "${filter_ml_ouput_folder}" ]; then mkdir -p ${filter_ml_ouput_folder}; filter_ml_ouput_folder="-d ${filter_ml_ouput_folder}";fi;
        if [ -n "${filter_ml_prefix}"  ]; then filter_ml_prefix="-o ${filter_ml_prefix}";fi
        if [ -n "${filter_ml_threshold}"  ]; then filter_ml_threshold="-j ${filter_ml_threshold}";fi
        if [ -n "${filter_ml_percent_training}"  ]; then filter_ml_percent_training="-t ${filter_ml_percent_training}";fi
        if [ -n "${filter_ml_remove_columns}"  ]; then filter_ml_remove_columns="-r ${filter_ml_remove_columns}";fi
        if [ -n "${filter_ml_intermediate_files}"  ]; then filter_ml_intermediate_files="--intermediate_files ${filter_ml_intermediate_files}";fi
        if [ -n "${filter_ml_max_class_size}"  ]; then filter_ml_max_class_size="-z ${filter_ml_max_class_size}";fi
        if [ -n "${filter_ml_intrapriming}"  ]; then filter_ml_intrapriming="-i ${filter_ml_intrapriming}";fi
        
    # Rescue parameters
    
        if [ -n "${reference_fasta}"  ]; then reference_fasta="-f ${QC_reference_fasta}";fi
        
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
        ${sqanti3_qc} ${QC_input} ${QC_reference_gtf}  ${QC_reference_fasta}   \
                             ${QC_fastq} ${QC_skip_ORF_param} ${QC_is_fusion_param} \
                             ${QC_cage_peak_bed_file} ${QC_polyA_motif_list} ${QC_output_prefix}  \
                             ${QC_destination_folder} ${QC_full_length_pacbio_abundance_tsv}   \
                             ${QC_short_reads_fofn} ${QC_cpus}  ${QC_report_file} \
                             ${QC_chunks} ${QC_min_ref_length} ${QC_force_id_ignore} \
                             ${QC_aligner_choice} ${QC_phylobed} \
                             ${QC_is_fusion} ${QC_orf_input} ${QC_is_fastq} \
                             ${QC_expression_matrix} ${QC_gmap_index} \
                             ${QC_coverage} ${QC_sites} ${QC_window} ${QC_genename} \
                             ${QC_saturation} 
                              
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
                
                
                ${sqanti3_filter} rules ${filter_rules_ouput_folder} \
                    ${filter_corrected_gtf} ${prefix_filter_rules} ${filter_isoforms} ${filter_isoannotgff3} \
                    ${filter_sam} ${filter_faa} ${filter_monoexonic} ${filter_skip_report} ${filter_rules_json_file} \
                    ${filter_input_classification}
            fi
            
            if [ -z ${filter_mode} ] || [ ${filter_mode} == "ml" ] || [ $filter_mode == "both" ];
            then 
                ${sqanti3_filter} ml ${filter_ml_ouput_folder} \
                ${filter_ml_prefix} ${filter_corrected_gtf} ${filter_isoforms} ${filter_isoannotgff3} \
                ${filter_sam} ${filter_faa} ${monoexonic} ${filter_monoexonic} ${filter_ml_percent_training} ${filter_ml_TP} \
                ${filter_ml_TN} ${filter_ml_threshold} \
                ${filter_ml_max_class_size} ${filter_ml_intermediate_files} ${filter_ml_intrapriming} \
                ${filter_input_classification} 
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
            ${sqanti3_rescue} ml ${rescue_ml_reference_genome} ${rescue_ml_isoforms} ${rescue_ml_gtf} ${rescue_ml_reference_classification}  \
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
