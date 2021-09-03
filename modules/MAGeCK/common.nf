///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             MAGeCK                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process MAGeCK_normalisation {
  tag "MAGeCK: run MAGeCK test to get normalised counts"

  publishDir "${params.resultDir}/normalised", mode: 'copy', pattern: "*normalised*"

  input:
    tuple val( count_type ), path( count_matrix ), path( library )
    val( norm_method )
    val( treatment_indices )

  output:
    tuple val( 'normalised' ), path( "${analysis_name}.normalised.tsv" ), path( library ), emit: normalised_count_matrix

  when:
    !no_normalisation

  script:
    analysis_name = "count_matrix.MAGeCK_${norm_method}"
    exec_name = "mageck test"

    cmd = "${exec_name} -k ${count_matrix}"

    cmd = "${cmd} -t ${treatment_indices}"
    cmd = "${cmd} -n ${analysis_name}"

    cmd = "${cmd} --norm-method ${norm_method}"
    cmd = "${cmd} --normcounts-to-file"

    cmd = "${cmd} --remove-zero none"
    cmd = "${cmd} --remove-zero-threshold 0"

    cmd = (params.mageck_extra_options) ? "${cmd} ${params.mageck_extra_options}" : cmd

    """
    $cmd

    mv "${analysis_name}.normalized.txt" "${analysis_name}.normalised.tsv"
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             MAGeCK                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process MAGeCK_test {
  tag "MAGeCK: MAGeCK test"

  publishDir "${params.resultDir}/MAGeCK/${contrast}", mode: 'copy', pattern: "*.txt"

  input:
    tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), path( sgrna_fold_change_matrix ), path( gene_fold_change_matrix )
    val( norm_method )
    val( analysis_indices )

  output:
    tuple val(contrast), path( "${analysis_name }.sgrna_summary.txt" ), emit: sgrna_summary
    tuple val(contrast), path( "${analysis_name }.gene_summary.txt" ) , emit: gene_summary

  when:
    !params.no_mageck && !params.no_analysis

  script:
    // Get contrast indices
    contrast_to_split = contrast
    (contrast_treatment, contrast_control) = contrast_to_split.split("_vs_", 2 )
    control_index_values = analysis_indices["${contrast}"]["${contrast_control}"]["base0"]
    treatment_index_values = analysis_indices["${contrast}"]["${contrast_treatment}"]["base0"]
    analysis_name = "MAGeCK.${contrast}"

    exec_name = "mageck test"

    cmd = "${exec_name} -k ${count_matrix} -c ${control_index_values} -t ${treatment_index_values} -n ${analysis_name} --norm-method ${norm_method}"

    cmd = (params.mageck_normcounts_to_file) ? "${cmd} --normcounts-to-file" : cmd
    cmd = "${cmd} --remove-zero ${params.mageck_remove_zero}"
    cmd = "${cmd} --remove-zero-threshold ${params.mageck_remove_zero_threshold}"
    cmd = (params.mageck_extra_options) ? "${cmd} ${params.mageck_extra_options}" : cmd

    """
    $cmd
    """
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             RESULTS                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process MAGeCK_process_results {
  tag "MAGeCK: MAGeCK test"

  publishDir "${params.resultDir}/MAGeCK/${gene_contrast}", mode: 'copy', pattern: "*.png"
  publishDir "${params.resultDir}/MAGeCK/${gene_contrast}", mode: 'copy', pattern: "mageck_rra*.tsv"


  input:
    tuple val( gene_contrast ), path( mageck_gene )
    tuple val( sgrna_contrast ), path( mageck_sgrna )

  output:
    path "MAGeCK*.png"
    path "mageck_rra*.tsv"

  when:
    !params.no_mageck && !params.no_analysis

  script:
    analysis_name = "${gene_contrast}"
    script_path = "${baseDir}/submodules/rcrispr/exec/process_mageck_rra_gene_summary.R"

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} --gene_summary ${mageck_gene}"
    cmd = "${cmd} --sgrna_summary ${mageck_sgrna}"
    cmd = "${cmd} --fdr ${params.mageck_fdr_threshold}"
    cmd = "${cmd} --n_genes ${params.mageck_n_genes}"

    cmd = "${cmd} --outdir \"${params.mageck_processed_results_outdir}\""
    cmd = "${cmd} --suffix \"${analysis_name}\""
    cmd = "${cmd} --rdata \"${params.mageck_processed_results_rdata}\""

    """
    $cmd
    """
}
