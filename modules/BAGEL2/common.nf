///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         BAGEL normalise counts                      -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
process bagel_normalise_counts {
  tag "BAGEL: normalisation"

  publishDir "${params.resultDir}/normalised", mode: 'copy', pattern: "${params.bagel_normalisation_outfile}", overwrite: true
  publishDir "${params.resultDir}/normalised/", mode: 'copy', pattern: "${params.bagel_normalisation_rdata}", overwrite: true

  input:
    tuple val(count_type), path(count_matrix), path( library )
    val count_indices

  output:
    tuple val('normalised'), path( "${params.bagel_normalisation_outfile}" ), path( library ), emit: normalised_count_matrix
    path "${params.bagel_normalisation_rdata}", emit: rdata

  when:
    params.normalisation_method == 'bagel' && !params.no_normalisation

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/BAGEL_normalisation.R"

    cmd = "${params.rscript_exec} ${script_path} -c ${count_matrix}"

    cmd = "${cmd} --outdir \"${params.bagel_normalisation_outdir}\""
    cmd = "${cmd} --outfile \"${params.bagel_normalisation_outfile}\""
    cmd = "${cmd} --rdata \"${params.bagel_normalisation_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"

    cmd = "${cmd} --pseudocount ${params.bagel_normalisation_pseudocount}"
    cmd = "${cmd} --scaling_factor ${params.bagel_normalisation_scaling_factor}"

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         BAGEL Bayes factor                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process BAGEL_bf {
  tag "BAGEL: BAGEL bf"

  publishDir "${params.resultDir}/BAGEL2/${contrast}", mode: 'copy', pattern: '*.bf', overwrite: true

  input:
    tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), path( sgrna_fold_change_matrix ), path( gene_fold_change_matrix )
    val( analysis_indices )
    val( analysis_type )

  output:
    tuple val( contrast ), path( "BAGEL.${contrast}.${analysis_type}.bf" ), emit: bagel_results

  when:
    !params.no_bagel && !params.no_analysis

  script:
    analysis_name = "BAGEL.${contrast}.${analysis_type}"
    exec_name = "BAGEL.py bf"

    // Get contrast indices
    if ( input_analysis_stage == 'corrected' ) {
      treatment_index_values = 7
    } else {
      treatment_index_values = analysis_indices["${contrast}"]["lfc"]["base1"]
    }

    output_filename = "${analysis_name}.bf"

    cmd = "${exec_name} -i ${sgrna_fold_change_matrix}"
    cmd = "${cmd} -o ${output_filename}"
    cmd = "${cmd} -e ${params.essential_genes}"
    cmd = "${cmd} -n ${params.nonessential_genes}"
    cmd = "${cmd} -c ${treatment_index_values}"
    if ( "${analysis_type}" == "sgrna" ) { cmd = "${cmd} -r" }
    cmd = (params.bagel_bf_extra_options) ? "${cmd} ${params.bagel_bf_extra_options}" : cmd

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         Calculate scaled BFs                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process scale_gene_BFs {
  tag "UTILS: calculating scaled BFs"

  publishDir "${params.resultDir}/BAGEL2/${contrast}/scaled", mode: 'copy', pattern: '*ROC*', overwrite: true
  publishDir "${params.resultDir}/BAGEL2/${contrast}/scaled", mode: 'copy', pattern: 'BF*', overwrite: true
  publishDir "${params.resultDir}/BAGEL2/${contrast}/scaled", mode: 'copy', pattern: '*scaled*', overwrite: true

  input:
    tuple val( contrast ), file( bagel_bf_gene )
    val( analysis_indices )

  output:
    tuple val( contrast ), path( "BF.gene.ROC.${analysis_suffix}.png" )
    tuple val( contrast ), path( "ROC_summary.${analysis_suffix}.tsv" )
    tuple val( contrast ), path( "BF.scaled.gene.${analysis_suffix}.tsv" )
    tuple val( contrast ), path( "BF.scaled_depletions_matrix.gene.${analysis_suffix}.tsv" )
    path "scaled_BF.${analysis_suffix}.Rdata", emit: rdata

  when:
    params.scale_bayes_factors

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/scale_lfcs_and_bfs.R"
    analysis_suffix = "BF.${contrast}"

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} --is_bf"
    cmd = "${cmd} --threshold ${params.scaled_bf_threshold}"

    cmd = "${cmd} --infile \"${bagel_bf_gene}\""
    cmd = ( params.scaled_bf_infile_header ) ? cmd : "${cmd} --no_infile_header"
    cmd = "${cmd} --infile_delim \"${params.scaled_bf_infile_delim}\""
    cmd = "${cmd} --infile_gene_column_index ${params.scaled_bf_infile_gene_column_index}"
    cmd = "${cmd} --infile_data_column_index \"${params.scaled_bf_infile_data_column_index}\""

    cmd = "${cmd} --ess \"${params.essential_genes}\""
    cmd = "${cmd} --ess_gene_column_index ${params.ess_gene_column_index}"
    cmd = "${cmd} --no_ess_header"
    cmd = "${cmd} --ess_delim \"${params.ess_delim}\""

    cmd = "${cmd} --noness \"${params.nonessential_genes}\""
    cmd = "${cmd} --noness_gene_column_index \"${params.noness_gene_column_index}\""
    cmd = "${cmd} --no_noness_header"
    cmd = "${cmd} --noness_delim \"${params.noness_delim}\""

    cmd = "${cmd} --outdir \"${params.scaled_bf_outdir}\""
    cmd = "${cmd} --suffix \"${analysis_suffix}\""
    cmd = "${cmd} --rdata \"${params.scaled_bf_rdata}\""

    """
    $cmd
    """
}
