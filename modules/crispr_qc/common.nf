///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   Raw count quality control (QC)                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process sequencing_qc {
    tag "QC: sequencing QC"
    errorStrategy 'ignore'

    publishDir "${params.resultDir}/QC/sequencing/", mode: 'copy'

    input:
      tuple val( analysis_stage ), file( count_matrix ), file( library )
      file sample_mapping
      val count_indices

    output:
      path "*.tsv"
      path "*.png"

    when:
      !params.no_qc

    script:
        script_path = "${baseDir}/submodules/rcrispr/exec/raw_qc_statistics.R"
        cmd = "${params.rscript_exec} ${script_path}"

        cmd = "${cmd} -c ${count_matrix}"
        cmd = "${cmd} -i ${sample_mapping}"

        cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
        cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
        cmd = "${cmd} --count_count_column_index ${count_indices}"
        cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
        cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

        cmd = ( params.info_header ) ? cmd : "${cmd} --no_info_header"
        cmd = "${cmd} --info_delim \"${params.info_delim}\""
        cmd = "${cmd} --info_filename_column_index ${params.info_filename_column_index}"
        cmd = "${cmd} --info_label_column_index ${params.info_label_column_index}"
        cmd = "${cmd} --info_plasmid_column_index ${params.info_plasmid_column_index}"
        cmd = "${cmd} --info_control_column_index ${params.info_control_column_index}"
        cmd = "${cmd} --info_treatment_column_index ${params.info_treatment_column_index}"
        cmd = ( params.info_group_column_index ) ? "${cmd} --info_group_column_index ${params.info_group_column_index}" : cmd
        cmd = ( params.info_reads_column_index ) ? "${cmd} --info_reads_column_index ${params.info_reads_column_index}" : cmd
        
        cmd = "${cmd} --outdir \"${params.raw_qc_outdir}\""
        cmd = "${cmd} --rdata \"${params.raw_qc_rdata}\""

        """
        $cmd
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                  Intermediate quality control (QC)                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process intermediate_qc {
    tag "QC: intermediate QC"
    errorStrategy 'ignore'

    publishDir "${params.resultDir}/QC/intermediate/${input_analysis_stage}", mode: 'copy', pattern: "*.png"
    publishDir "${params.resultDir}/QC/intermediate/${input_analysis_stage}", mode: 'copy', pattern: "intermediate_summary*.tsv"

    input:
      tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), file( sgrna_fold_change_matrix ), file( gene_fold_change_matrix ), file( sample_mapping )
      val( analysis_indices )  
      val( counts_or_fc )  

    output:
      path "*.png"
      path "intermediate_summary*.tsv"

    when:
      !params.no_qc

    script:
        analysis_name = "${input_analysis_stage}.${contrast}"
        infile_id_column_index = 0
        no_check_names = false
        // Whether to use gene/sgrna or counts/fold change input
        if (counts_or_fc == "counts") {
            input_data = "${count_matrix}"
            count_indices = analysis_indices["${contrast}"]["count_lfc"]["base1_increment2"]
            infile_id_column_index = 1
            infile_gene_column_index = 2
        } else {
            input_data = "${gene_fold_change_matrix}"
            if (input_analysis_stage == "corrected") {
                count_indices = 2
                no_check_names = true
            } else {
                count_indices = analysis_indices["${contrast}"]["lfc"]["base1_increment1"]
            }
            infile_gene_column_index = 1
        }
 
        script_path = "${baseDir}/submodules/rcrispr/exec/intermediate_qc.R"

        cmd = "${params.rscript_exec} ${script_path}"
         cmd = ( counts_or_fc == "fc" ) ? "${cmd} --infile ${gene_fold_change_matrix}" : "${cmd} --infile ${count_matrix}"
        cmd = "${cmd} --info ${sample_mapping}"

        cmd = ( counts_or_fc == "fc" ) ? "${cmd} --is_gene" : "${cmd}"
        cmd = ( counts_or_fc == "fc" ) ? "${cmd} --is_fc" : "${cmd}"
        cmd = ( no_check_names ) ? "${cmd} --no_check_names" : "${cmd}"

        cmd = "${cmd} --outdir \"${params.intermediate_qc_outdir}\""
        cmd = "${cmd} --outfile \"${params.intermediate_qc_outfile}\""
        cmd = "${cmd} --rdata \"${params.intermediate_qc_rdata}\""
        cmd = "${cmd} --suffix \"${analysis_name}\""

        cmd = ( params.intermediate_qc_infile_header ) ? cmd : "${cmd} --no_infile_header"
        cmd = "${cmd} --infile_delim \"${params.intermediate_qc_infile_delim}\""
        cmd = "${cmd} --infile_id_column_index ${infile_id_column_index}"
        cmd = "${cmd} --infile_gene_column_index ${infile_gene_column_index}"
        cmd = "${cmd} --infile_data_column_index \"${count_indices}\""

        cmd = ( params.info_header ) ? cmd : "${cmd} --no_info_header"
        cmd = "${cmd} --info_delim \"${params.info_delim}\""
        cmd = "${cmd} --info_filename_column_index ${params.info_filename_column_index}"
        cmd = "${cmd} --info_label_column_index ${params.info_label_column_index}"
        cmd = "${cmd} --info_plasmid_column_index ${params.info_plasmid_column_index}"
        cmd = "${cmd} --info_control_column_index ${params.info_control_column_index}"
        cmd = "${cmd} --info_treatment_column_index ${params.info_treatment_column_index}"
        cmd = ( params.info_group_column_index ) ? "${cmd} --info_group_column_index ${params.info_group_column_index}" : cmd

        """
        $cmd
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --      Intermediate BAGEL classification quality control (QC)         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process intermediate_bagel_classification_qc {
    tag "QC: intermediate BAGEL classification QC"
    errorStrategy 'ignore'

    publishDir "${params.resultDir}/QC/intermediate/${input_analysis_stage}", mode: 'copy', pattern: "*.png"
    publishDir "${params.resultDir}/QC/intermediate/${input_analysis_stage}", mode: 'copy', pattern: "bagel_classification_summary*.tsv"

    input:
      tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), file( sgrna_fold_change_matrix ), file( gene_fold_change_matrix ), file( sample_mapping )
      val( analysis_indices )  
      val( counts_or_fc )  

    output:
      path "*.png"
      path "bagel_classification_summary*.tsv"

    when:
      !params.no_qc

    script:
        analysis_name = "${input_analysis_stage}.${contrast}"
        infile_id_column_index = 0
        no_check_names = false
        // Whether to use gene/sgrna or counts/fold change input
        if (counts_or_fc == "counts") {
            input_data = "${count_matrix}"
            count_indices = analysis_indices["${contrast}"]["count_lfc"]["base1_increment2"]
            infile_id_column_index = 1
            infile_gene_column_index = 2
        } else {
            input_data = "${gene_fold_change_matrix}"
            if (input_analysis_stage == "corrected") {
                count_indices = 2
                no_check_names = true
            } else {
                count_indices = analysis_indices["${contrast}"]["lfc"]["base1_increment1"]
            }
            infile_gene_column_index = 1
        }
 
        script_path = "${baseDir}/submodules/rcrispr/exec/bagel_classification_qc.R"

        cmd = "${params.rscript_exec} ${script_path}"
         cmd = ( counts_or_fc == "fc" ) ? "${cmd} --infile ${gene_fold_change_matrix}" : "${cmd} --infile ${count_matrix}"
        cmd = "${cmd} --info ${sample_mapping}"

        cmd = ( counts_or_fc == "fc" ) ? "${cmd} --is_gene" : "${cmd}"
        cmd = ( counts_or_fc == "fc" ) ? "${cmd} --is_fc" : "${cmd}"
        cmd = ( no_check_names ) ? "${cmd} --no_check_names" : "${cmd}"

        cmd = "${cmd} --outdir \"${params.bagel_classification_qc_outdir}\""
        cmd = "${cmd} --outfile \"${params.bagel_classification_qc_outfile}\""
        cmd = "${cmd} --rdata \"${params.bagel_classification_qc_rdata}\""
        cmd = "${cmd} --suffix \"${analysis_name}\""

        cmd = ( params.bagel_classification_qc_infile_header ) ? cmd : "${cmd} --no_infile_header"
        cmd = "${cmd} --infile_delim \"${params.bagel_classification_infile_delim}\""
        cmd = "${cmd} --infile_id_column_index ${infile_id_column_index}"
        cmd = "${cmd} --infile_gene_column_index ${infile_gene_column_index}"
        cmd = "${cmd} --infile_data_column_index \"${count_indices}\""

        cmd = ( params.info_header ) ? cmd : "${cmd} --no_info_header"
        cmd = "${cmd} --info_delim \"${params.info_delim}\""
        cmd = "${cmd} --info_filename_column_index ${params.info_filename_column_index}"
        cmd = "${cmd} --info_label_column_index ${params.info_label_column_index}"
        cmd = "${cmd} --info_plasmid_column_index ${params.info_plasmid_column_index}"
        cmd = "${cmd} --info_control_column_index ${params.info_control_column_index}"
        cmd = "${cmd} --info_treatment_column_index ${params.info_treatment_column_index}"
        cmd = ( params.info_group_column_index ) ? "${cmd} --info_group_column_index ${params.info_group_column_index}" : cmd

        cmd = "${cmd} --ess \"${params.essential_genes}\""
        cmd = "${cmd} --ess_gene_column_index ${params.ess_gene_column_index}"
        cmd = "${cmd} --no_ess_header"
        cmd = "${cmd} --ess_delim \"${params.ess_delim}\""

        cmd = "${cmd} --noness \"${params.nonessential_genes}\""
        cmd = "${cmd} --noness_gene_column_index \"${params.noness_gene_column_index}\""
        cmd = "${cmd} --no_noness_header"
        cmd = "${cmd} --noness_delim \"${params.noness_delim}\""

        """
        $cmd
        """
}