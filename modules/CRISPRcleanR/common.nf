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
/* --                Normalise counts with CRISPRcleanR                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process crisprcleanr_normalise_counts {
  tag "CRISPRcleanR: normalisation"

  publishDir "${params.resultDir}/normalised", mode: 'copy', pattern: "${params.crisprcleanr_normalisation_count_matrix_outfile}", overwrite: true
  publishDir "${params.resultDir}/normalised/", mode: 'copy', pattern: "${params.crisprcleanr_normalisation_fc_matrix_outfile}", overwrite: true
  publishDir "${params.resultDir}/normalised/", mode: 'copy', pattern: "${params.crisprcleanr_normalisation_library_outfile}", overwrite: true
  publishDir "${params.resultDir}/normalised/", mode: 'copy', pattern: "${params.crisprcleanr_normalisation_rdata}", overwrite: true

  input:
    tuple val( count_type ), path( count_matrix ), path( library )
    val ncontrols
    val count_indices

  output:
    tuple val('normalised'), path( "${params.crisprcleanr_normalisation_count_matrix_outfile}" ), path( "${params.crisprcleanr_normalisation_library_outfile}" ), emit: normalised_count_matrix
    path "${params.crisprcleanr_normalisation_fc_matrix_outfile}", emit: fold_changes
    path "${params.crisprcleanr_normalisation_rdata}", emit: rdata

  when:
    !params.no_normalisation && params.normalisation_method == 'crisprcleanr'

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/CRISPRcleanR_normalisation.R"

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} -c ${count_matrix}"
    cmd = "${cmd} -l ${library}"
    cmd = "${cmd} --min_reads ${params.min_reads}"
    cmd = "${cmd} --n_controls ${ncontrols}"

    cmd = "${cmd} --outdir \"${params.crisprcleanr_normalisation_outdir}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.crisprcleanr_normalisation_count_matrix_outfile}\""
    cmd = "${cmd} --lfc_matrix_outfile \"${params.crisprcleanr_normalisation_fc_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.crisprcleanr_normalisation_library_outfile}\""
    cmd = "${cmd} --rdata \"${params.crisprcleanr_normalisation_rdata}\""

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
/* --            Format library and counts for CRISPRcleanR               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process format_library_and_matrices_for_crisprcleanr {
  tag "CRISPRcleanR: format library and counts for CRISPRcleanR"

  publishDir "${params.resultDir}/corrected/CRISPRcleanR/formatted_inputs", mode: 'copy', pattern: "*_input*.tsv", overwrite: true
  //publishDir "${params.resultDir}/corrected/CRISPRcleanR/formatted_inputs", mode: 'copy', pattern: "${params.crisprcleanr_formatting_rdata}", overwrite: true

  input:
    tuple val( analysis_stage ), val( contrast ), file( library ), path( fc_input_count_matrix ), path( count_matrix ), path( sgrna_fold_change_matrix ), path( gene_fold_change_matrix )
    val( analysis_indices )

  output:
    tuple val( analysis_stage ), val( contrast ), file( "library.${file_suffix}.tsv" ), path( fc_input_count_matrix ), path( "count_matrix.${file_suffix}.tsv" ), path( "fold_change_matrix.sgrna.${file_suffix}.tsv" ), path( gene_fold_change_matrix ), emit: inputdata
    //path "crisprcleanr_formatting.${file_suffix}.Rdata", emit: rdata

  when:
    !params.no_correction && !params.no_crisprcleanr

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/format_library_and_matrices_for_CRISPRcleanR.R"
    file_suffix = "CRISPRcleanR_input.${contrast}.${analysis_stage}"
    count_indices = analysis_indices["${contrast}"]["count_lfc"]["base1_increment2"]
    lfc_indices = analysis_indices["${contrast}"]["lfc"]["base1_increment2"]

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} --counts ${count_matrix}"
    cmd = "${cmd} --lfc_matrix ${sgrna_fold_change_matrix}"
    cmd = "${cmd} --library ${library}"

    cmd = "${cmd} --outdir \"${params.crisprcleanr_formatting_outdir}\""
    cmd = "${cmd} --suffix \"${file_suffix}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.crisprcleanr_formatting_count_matrix_outfile}\""
    cmd = "${cmd} --lfc_matrix_outfile \"${params.crisprcleanr_formatting_fc_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.crisprcleanr_formatting_library_outfile}\""
    cmd = "${cmd} --rdata \"${params.crisprcleanr_formatting_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --lfc_id_column_index ${params.processed_lfc_id_column_index}"
    cmd = "${cmd} --lfc_gene_column_index ${params.processed_lfc_gene_column_index}"
    cmd = "${cmd} --lfc_lfc_column_index ${lfc_indices}"
    cmd = ( params.processed_lfc_header ) ? cmd : "${cmd} --no_lfc_header"
    cmd = "${cmd} --lfc_delim \"${params.processed_lfc_delim}\""

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
/* --                       CRISPRcleanR correction                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process crisprcleanr_correction {

  tag "CRISPRcleanR: CRISPRcleanR correction"

  publishDir "${params.resultDir}/corrected/CRISPRcleanR/output", mode: 'copy', pattern: "*CRISPRcleanR_corrected*", overwrite: true
  //publishDir "${params.resultDir}/corrected/CRISPRcleanR/output", mode: 'copy', pattern: "*.Rdata", overwrite: true

  input:
    tuple val( input_analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( count_matrix ), path( sgrna_fold_change_matrix ), path( gene_fold_change_matrix )
    val( analysis_stage )
    val(analysis_indices)

  output:
    tuple val( analysis_stage ), val( contrast ), file( library ), file( fc_input_count_matrix ), file( "count_matrix.${file_suffix}.tsv" ), file( "fold_change_matrix.sgRNA.${file_suffix}.tsv" ), file( "fold_change_matrix.gene.${file_suffix}.tsv" ), emit: fold_change_matrix

  when:
    !params.no_correction && !params.no_crisprcleanr

  script:
    script_path = "${baseDir}/submodules/rcrispr/exec/CRISPRcleanR_correction.R"
    file_suffix = "CRISPRcleanR_corrected.${contrast}.${analysis_stage}"
    count_indices = analysis_indices["${contrast}"]["count_lfc"]["base1_increment2"]
    lfc_indices = analysis_indices["${contrast}"]["lfc"]["base1_increment2"]

    cmd = "${params.rscript_exec} ${script_path}"
    cmd = "${cmd} --counts ${count_matrix}"
    cmd = "${cmd} --lfc_matrix ${sgrna_fold_change_matrix}"
    cmd = "${cmd} --library ${library}"

    cmd = "${cmd} --outdir \"${params.crisprcleanr_correction_outdir}\""
    cmd = "${cmd} --suffix \"${file_suffix}\""
    cmd = "${cmd} --count_matrix_outfile \"${params.crisprcleanr_correction_count_matrix_outfile}\""
    cmd = "${cmd} --lfc_matrix_outfile \"${params.crisprcleanr_correction_lfc_matrix_outfile}\""
    cmd = "${cmd} --lfc_gene_matrix_outfile \"${params.crisprcleanr_correction_lfc_gene_matrix_outfile}\""
    cmd = "${cmd} --library_outfile \"${params.crisprcleanr_correction_library_outfile}\""
    cmd = "${cmd} --rdata \"${params.crisprcleanr_correction_rdata}\""

    cmd = "${cmd} --count_id_column_index ${params.processed_count_id_column_index}"
    cmd = "${cmd} --count_gene_column_index ${params.processed_count_gene_column_index}"
    cmd = "${cmd} --count_count_column_index ${count_indices}"
    cmd = ( params.processed_counts_header ) ? cmd : "${cmd} --no_counts_header"
    cmd = "${cmd} --counts_delim \"${params.processed_counts_delim}\""

    cmd = "${cmd} --lfc_id_column_index ${params.processed_lfc_id_column_index}"
    cmd = "${cmd} --lfc_gene_column_index ${params.processed_lfc_gene_column_index}"
    cmd = "${cmd} --lfc_lfc_column_index ${lfc_indices}"
    cmd = ( params.processed_lfc_header ) ? cmd : "${cmd} --no_lfc_header"
    cmd = "${cmd} --lfc_delim \"${params.processed_lfc_delim}\""

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
