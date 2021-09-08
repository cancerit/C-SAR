/*
 * Copyright (c) 2021 Genome Research Ltd
 *
 * Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
 *
 * This file is part of C-SAR.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * 1. The usage of a range of years within a copyright statement contained within
 * this distribution should be interpreted as being equivalent to a list of years
 * including the first and last year specified and all consecutive years between
 * them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
 * 2009, 2011-2012’ should be interpreted as being identical to a statement that
 * reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
 * statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
 * identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
 * 2009, 2010, 2011, 2012’.
 *
 */
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --             Convert count files to matrix                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process counts2matrix {
    tag "UTILS: convert count files to matrix"

    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "${params.counts_to_matrix_count_matrix_outfile}", overwrite: true
    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "${params.counts_to_matrix_library_outfile}", overwrite: true
    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "*.Rdata", overwrite: true

    input:
      path counts_directory
      tuple val( library_label ), file( library )
      file sample_mapping

    output:
      tuple val( 'raw' ), path ( "${params.counts_to_matrix_count_matrix_outfile}" ), path( "${params.counts_to_matrix_library_outfile}" ), emit: data
      path "${params.counts_to_matrix_rdata}", emit: rdata

    script:
      script_path = "${baseDir}/submodules/rcrispr/exec/sample_counts_to_matrix.R"

      cmd = "${params.rscript_exec} ${script_path} -c ${counts_directory} -l ${library} -i ${sample_mapping}"

      cmd = "${cmd} --outdir \"${params.counts_to_matrix_outdir}\""
      cmd = "${cmd} --count_matrix_outfile \"${params.counts_to_matrix_count_matrix_outfile}\""
      cmd = "${cmd} --library_outfile \"${params.counts_to_matrix_library_outfile}\""
      cmd = "${cmd} --rdata \"${params.counts_to_matrix_rdata}\""

      cmd = ( params.counts_header ) ? cmd : "${cmd} --no_counts_header"
      cmd = "${cmd} --counts_delim \"${params.counts_delim}\""
      cmd = "${cmd} --count_skip ${params.count_skip}"
      cmd = "${cmd} --count_id_column_index ${params.count_id_column_index}"
      cmd = "${cmd} --count_gene_column_index ${params.count_gene_column_index}"
      cmd = "${cmd} --count_count_column_index ${params.count_count_column_index}"

      cmd = ( params.library_header ) ? cmd : "${cmd} --no_library_header"
      cmd = "${cmd} --library_delim \"${params.library_delim}\""
      cmd = "${cmd} --library_id_column_index ${params.library_id_column_index}"
      cmd = "${cmd} --library_gene_column_index ${params.library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.library_end_column_index}" : cmd

      cmd = ( params.info_header ) ? cmd : "${cmd} --no_info_header"
      cmd = "${cmd} --info_delim \"${params.info_delim}\""
      cmd = "${cmd} --info_filename_column_index ${params.info_filename_column_index}"
      cmd = "${cmd} --info_label_column_index ${params.info_label_column_index}"
      cmd = "${cmd} --info_plasmid_column_index ${params.info_plasmid_column_index}"
      cmd = "${cmd} --info_control_column_index ${params.info_control_column_index}"
      cmd = "${cmd} --info_treatment_column_index ${params.info_treatment_column_index}"

      """
      $cmd
      """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          Read count matrix                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process read_count_matrix {
    tag "UTILS: read count matrix"

    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "${params.counts_to_matrix_count_matrix_outfile}", overwrite: true
    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "${params.counts_to_matrix_library_outfile}", overwrite: true
    publishDir "${params.resultDir}/raw/", mode: 'copy', pattern: "*.Rdata", overwrite: true

    input:
      path counts_directory
      tuple val( library_label ), file( library )
      file sample_mapping

    output:
      tuple val( 'raw' ), path ( "${params.counts_to_matrix_count_matrix_outfile}" ), path( "${params.counts_to_matrix_library_outfile}" ), emit: data
      path "${params.counts_to_matrix_rdata}", emit: rdata

    script:
      script_path = "${baseDir}/submodules/rcrispr/exec/read_sample_count_matrix.R"

      cmd = "${params.rscript_exec} ${script_path} -c ${counts_directory} -l ${library} -i ${sample_mapping}"

      cmd = "${cmd} --outdir \"${params.counts_to_matrix_outdir}\""
      cmd = "${cmd} --count_matrix_outfile \"${params.counts_to_matrix_count_matrix_outfile}\""
      cmd = "${cmd} --library_outfile \"${params.counts_to_matrix_library_outfile}\""
      cmd = "${cmd} --rdata \"${params.counts_to_matrix_rdata}\""

      cmd = ( params.counts_header ) ? cmd : "${cmd} --no_counts_header"
      cmd = "${cmd} --counts_delim \"${params.counts_delim}\""
      cmd = "${cmd} --count_id_column_index ${params.count_id_column_index}"
      cmd = "${cmd} --count_gene_column_index ${params.count_gene_column_index}"
      cmd = "${cmd} --count_count_column_index ${params.count_count_column_index}"

      cmd = ( params.library_header ) ? cmd : "${cmd} --no_library_header"
      cmd = "${cmd} --library_delim \"${params.library_delim}\""
      cmd = "${cmd} --library_id_column_index ${params.library_id_column_index}"
      cmd = "${cmd} --library_gene_column_index ${params.library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.library_end_column_index}" : cmd

      cmd = ( params.info_header ) ? cmd : "${cmd} --no_info_header"
      cmd = "${cmd} --info_delim \"${params.info_delim}\""
      cmd = "${cmd} --info_filename_column_index ${params.info_filename_column_index}"
      cmd = "${cmd} --info_label_column_index ${params.info_label_column_index}"
      cmd = "${cmd} --info_plasmid_column_index ${params.info_plasmid_column_index}"
      cmd = "${cmd} --info_control_column_index ${params.info_control_column_index}"
      cmd = "${cmd} --info_treatment_column_index ${params.info_treatment_column_index}"

      """
      $cmd
      """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --     Remove user defined guides from library and count matrix        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process remove_user_defined_guides {
  tag "REMOVE GUIDES: removing user-defined guides"

  publishDir "${params.resultDir}/user_removed", mode: 'copy', pattern: "${params.remove_user_guides_count_matrix_outfile}", overwrite: true
  publishDir "${params.resultDir}/user_removed", mode: 'copy', pattern: "${params.remove_user_guides_library_outfile}", overwrite: true
  publishDir "${params.resultDir}/user_removed/", mode: 'copy', pattern: "${params.remove_user_guides_rdata}", overwrite: true

  input:
    tuple val( analysis_stage ), file( count_matrix ), file( library )
    val count_indices

  output:
    tuple val('user_removed'), path( "${params.remove_user_guides_count_matrix_outfile}" ), path( "${params.remove_user_guides_library_outfile}" ), emit: user_removed_counts
    path "${params.remove_user_guides_rdata}", emit: rdata

  when:
    params.remove_guides

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/remove_guides.R"

    cmd = "${params.rscript_exec} ${script_path} -c ${count_matrix} -l ${library}"
    cmd = "${cmd} --guides_to_remove ${params.remove_guides}"

    cmd = "${cmd} --outdir \"${params.remove_user_guides_outdir}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.remove_user_guides_count_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.remove_user_guides_library_outfile}\""
    cmd = "${cmd} --rdata \"${params.remove_user_guides_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --library_id_column_index ${params.processed_library_id_column_index}"
    cmd = "${cmd} --library_gene_column_index ${params.processed_library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.processed_library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.processed_library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.processed_library_end_column_index}" : cmd
    cmd = ( params.processed_library_header ) ? cmd : "${cmd} --no_library_header"
    cmd = "${cmd} --library_delim \"${params.processed_library_delim}\""

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --       Remove duplicate guides from library and count matrix         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process remove_duplicate_guides {
  tag "REMOVE DUPLICATES: removing guides with duplicate ids"

  publishDir "${params.resultDir}/duplicates_removed", mode: 'copy', pattern: "${params.remove_duplicate_guides_count_matrix_outfile}", overwrite: true
  publishDir "${params.resultDir}/duplicates_removed", mode: 'copy', pattern: "${params.remove_duplicate_guides_library_outfile}", overwrite: true
  publishDir "${params.resultDir}/duplicates_removed/", mode: 'copy', pattern: "${params.remove_duplicate_guides_outfile}", overwrite: true
  publishDir "${params.resultDir}/duplicates_removed/", mode: 'copy', pattern: "${params.remove_duplicate_guides_rdata}", overwrite: true

  input:
    tuple val( analysis_stage ), file( count_matrix ), file( library )
    val count_indices

  output:
    tuple val('duplicates_removed'), path( "${params.remove_duplicate_guides_count_matrix_outfile}" ), path( "${params.remove_duplicate_guides_library_outfile}" ), emit: duplicates_removed_counts
    tuple val('duplicates_removed'), path( "${params.remove_duplicate_guides_outfile}" ), emit: duplicate_guides
    path "${params.remove_duplicate_guides_rdata}", emit: rdata

  when:
    !params.no_duplicate_removal

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/remove_duplicate_guides.R"

    cmd = "${params.rscript_exec} ${script_path} -c ${count_matrix} -l ${library}"

    cmd = "${cmd} --outdir \"${params.remove_duplicate_guides_outdir}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.remove_duplicate_guides_count_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.remove_duplicate_guides_library_outfile}\""
    cmd = "${cmd} --duplicate_guides_outfile \"${params.remove_duplicate_guides_outfile}\""
    cmd = "${cmd} --rdata \"${params.remove_duplicate_guides_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --library_id_column_index ${params.processed_library_id_column_index}"
    cmd = "${cmd} --library_gene_column_index ${params.processed_library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.processed_library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.processed_library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.processed_library_end_column_index}" : cmd
    cmd = ( params.processed_library_header ) ? cmd : "${cmd} --no_library_header"
    cmd = "${cmd} --library_delim \"${params.processed_library_delim}\""

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --              Filter counts based on column indices                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process filter_counts_by_indices {
  tag "FILTER: sample indices < min reads"

  publishDir "${params.resultDir}/${analysis_stage}", mode: 'copy', pattern: "${params.raw_filter_count_matrix_outfile}", overwrite: true
  publishDir "${params.resultDir}/${analysis_stage}", mode: 'copy', pattern: "${params.raw_filter_library_outfile}", overwrite: true
  publishDir "${params.resultDir}/${analysis_stage}", mode: 'copy', pattern: "${params.raw_filter_outfile}", overwrite: true
  publishDir "${params.resultDir}/${analysis_stage}", mode: 'copy', pattern: "${params.raw_filter_rdata}", overwrite: true

  input:
    tuple val( count_label ), file( count_matrix ), file( library )
    val( index_values )
    val( analysis_stage )
    val count_indices

  output:
    tuple val( analysis_stage ), path( "${params.raw_filter_count_matrix_outfile}" ), path( "${params.raw_filter_library_outfile}" ), emit: filtered_counts
    tuple val( analysis_stage ), path( "${params.raw_filter_outfile}" ), emit: filtered_guides
    path "${params.raw_filter_rdata}", emit: rdata

  when:
    !params.no_filter

  script:

    script_path = "${baseDir}/submodules/rcrispr/exec/filter_by_column_indices.R"

    cmd = "${params.rscript_exec} ${script_path} -c ${count_matrix} -l ${library}"

    cmd = "${cmd} --min_reads ${params.min_reads}"
    cmd = "${cmd} --filter_method \"${params.filter_method}\""
    cmd = "${cmd} --filter_indices \"${index_values}\""

    cmd = "${cmd} --outdir \"${params.raw_filter_outdir}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.raw_filter_count_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.raw_filter_library_outfile}\""
    cmd = "${cmd} --filtered_guides_outfile \"${params.raw_filter_outfile}\""
    cmd = "${cmd} --rdata \"${params.raw_filter_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --library_id_column_index ${params.processed_library_id_column_index}"
    cmd = "${cmd} --library_gene_column_index ${params.processed_library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.processed_library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.processed_library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.processed_library_end_column_index}" : cmd
    cmd = ( params.processed_library_header ) ? cmd : "${cmd} --no_library_header"
    cmd = "${cmd} --library_delim \"${params.processed_library_delim}\""

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                Remove guides with no coordinates                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process remove_guides_with_no_coordinates {
  tag "UTILS: remove guides with no coordinates"

  publishDir "${params.resultDir}/corrected/CRISPRcleanR/inputs", mode: 'copy', pattern: "${params.coord_filter_count_matrix_outfile}"
  publishDir "${params.resultDir}/corrected/CRISPRcleanR/inputs", mode: 'copy', pattern: "${params.coord_filter_library_outfile}"
  publishDir "${params.resultDir}/corrected/CRISPRcleanR/inputs", mode: 'copy', pattern: "${params.coord_filter_outfile}"
  publishDir "${params.resultDir}/corrected/CRISPRcleanR/inputs", mode: 'copy', pattern: "${params.coord_filter_rdata}"

  input:
    tuple val(count_type), path(count_matrix), path(library)
    val( analysis_stage )
    val count_indices

  output:
    tuple val( analysis_stage ), path( "${params.coord_filter_count_matrix_outfile}" ), path( "${params.coord_filter_library_outfile}" ) , emit: count_matrix
    tuple val( analysis_stage ), path( "${params.coord_filter_outfile}" ), emit: excluded_guides
    tuple val( analysis_stage ), path( "${params.coord_filter_rdata}" ), emit: rdata

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/remove_guides_without_coordinates.R"

    cmd = "${params.rscript_exec} ${script_path} -c ${count_matrix} -l ${library}"

    cmd = "${cmd} --outdir \"${params.coord_filter_outdir}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.coord_filter_count_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.coord_filter_library_outfile}\""
    cmd = "${cmd} --excluded_guides_outfile \"${params.coord_filter_outfile}\""
    cmd = "${cmd} --rdata \"${params.coord_filter_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --library_id_column_index ${params.processed_library_id_column_index}"
    cmd = "${cmd} --library_gene_column_index ${params.processed_library_gene_column_index}"
    cmd = ( params.library_chr_column_index ) ? "${cmd} --library_chr_column_index ${params.processed_library_chr_column_index}" : cmd
    cmd = ( params.library_start_column_index ) ? "${cmd} --library_start_column_index ${params.processed_library_start_column_index}" : cmd
    cmd = ( params.library_end_column_index ) ? "${cmd} --library_end_column_index ${params.processed_library_end_column_index}" : cmd
    cmd = ( params.processed_library_header ) ? cmd : "${cmd} --no_library_header"
    cmd = "${cmd} --library_delim \"${params.processed_library_delim}\""

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   Calculate log fold changes                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process calculate_log_fold_changes {
  tag "UTILS: calculating log fold changes"

  publishDir "${params.resultDir}/${analysis_stage}/${contrast}", mode: 'copy', pattern: "*_vs_*.tsv"
  //publishDir "${params.resultDir}/${analysis_stage}/${contrast}", mode: 'copy', pattern: "*.Rdata"


  input:
    tuple val( analysis_stage ), file( count_matrix ), file( library ), val( contrast )
    val( analysis_indices )

  output:
    tuple val( analysis_stage ), val( contrast ), path( library ), path( count_matrix ), path( "count_matrix.lfc.${contrast}.${analysis_stage}.tsv" ), path( "fold_change_matrix.sgrna.lfc.${contrast}.${analysis_stage}.tsv" ), path( "fold_change_matrix.gene.lfc.${contrast}.${analysis_stage}.tsv" ), emit: contrast_fold_changes
    //tuple val( analysis_stage ), path( "calculate_log_fc.lfc.${contrast}.${analysis_stage}.Rdata" ), emit: rdata

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/calculate_log_fold_changes.R"

    // Get contrast indices
    contrast_to_split = contrast
    (contrast_treatment, contrast_control) = contrast_to_split.split("_vs_", 2 )
    control_index_values = analysis_indices["${contrast_control}"]["base1_increment2"]
    treatment_index_values = analysis_indices["${contrast_treatment}"]["base1_increment2"]
    count_indices = analysis_indices["all"]["base1_increment2"]

    output_file_suffix = "lfc.${contrast}.${analysis_stage}"

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} -c ${count_matrix}"

    cmd = "${cmd} --control_indices \"$control_index_values\""
    cmd = "${cmd} --treatment_indices \"$treatment_index_values\""

    cmd = "${cmd} --outdir \"${params.calculate_log_fc_outdir}\""
    cmd = "${cmd} --suffix \"${output_file_suffix}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.calculate_log_fc_count_matrix_outfile}\""
    cmd = "${cmd} --sgrna_outfile \"${params.calculate_log_fc_sgrna_outfile}\""
    cmd = "${cmd} --gene_outfile \"${params.calculate_log_fc_gene_outfile}\""
    cmd = "${cmd} --rdata \"${params.calculate_log_fc_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index \"${count_indices}\""
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    """
    $cmd
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                 Calculate scaled log fold changes                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process scale_gene_log_fold_changes {
  tag "UTILS: calculating scaled log fold changes"

  publishDir "${params.resultDir}/${input_analysis_stage}/${contrast}/scaled", mode: 'copy', pattern: '*ROC*'
  publishDir "${params.resultDir}/${input_analysis_stage}/${contrast}/scaled", mode: 'copy', pattern: '*scaled*'
  publishDir "${params.resultDir}/${input_analysis_stage}/${contrast}/scaled", mode: 'copy', pattern: '*Rdata'

  input:
    tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), file( sgrna_fold_change_matrix ), file( gene_fold_change_matrix )
    val( analysis_indices )

  output:
    tuple val( input_analysis_stage ), val( contrast ), file( gene_fold_change_matrix), path( "fold_change.gene.ROC.${analysis_suffix}.png" )
    tuple val( input_analysis_stage ), val( contrast ), file( gene_fold_change_matrix), path( "ROC_summary.${analysis_suffix}.tsv" )
    tuple val( input_analysis_stage ), val( contrast ), file( gene_fold_change_matrix), path( "fold_change.scaled.gene.${analysis_suffix}.tsv" )
    tuple val( input_analysis_stage ), val( contrast ), file( gene_fold_change_matrix), path( "fold_change.scaled_depletions_matrix.gene.${analysis_suffix}.tsv" )
    path "scaled_lfc.${analysis_suffix}.Rdata", emit: rdata

  when:
    params.scale_log_fold_changes

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/scale_lfcs_and_bfs.R"
    analysis_suffix = "LFC.${input_analysis_stage}.${contrast}"

    // Get contrast indices
    // Need a specific clause for CRISPRcleanR
    // TODO: put this check into main instead of within the module (to allow for other correction software)
    if ( input_analysis_stage == 'corrected' ) {
      treatment_index_values = 2
    } else {
      contrast_to_split = contrast
      (contrast_treatment, contrast_control) = contrast_to_split.split("_vs_", 2 )
      treatment_index_values = analysis_indices["${contrast}"]["lfc"]["base1_increment1"]
    }

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} --is_fc"
    cmd = "${cmd} --threshold ${params.scaled_lfc_threshold}"

    cmd = "${cmd} --infile \"${gene_fold_change_matrix}\""
    cmd = ( params.scaled_lfc_infile_header ) ? cmd : "${cmd} --no_infile_header"
    cmd = "${cmd} --infile_delim \"${params.scaled_lfc_infile_delim}\""
    cmd = "${cmd} --infile_gene_column_index ${params.scaled_lfc_infile_gene_column_index}"
    cmd = "${cmd} --infile_data_column_index \"$treatment_index_values\""

    cmd = "${cmd} --ess \"${params.essential_genes}\""
    cmd = "${cmd} --ess_gene_column_index ${params.ess_gene_column_index}"
    cmd = "${cmd} --no_ess_header"
    cmd = "${cmd} --ess_delim \"${params.ess_delim}\""

    cmd = "${cmd} --noness \"${params.nonessential_genes}\""
    cmd = "${cmd} --noness_gene_column_index ${params.noness_gene_column_index}"
    cmd = "${cmd} --no_noness_header"
    cmd = "${cmd} --noness_delim \"${params.noness_delim}\""

    cmd = "${cmd} --outdir \"${params.scaled_lfc_outdir}\""
    cmd = "${cmd} --suffix \"${analysis_suffix}\""
    cmd = "${cmd} --rdata \"${params.scaled_lfc_rdata}\""

    """
    $cmd
    """
}
