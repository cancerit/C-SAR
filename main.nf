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
/* --                                 DSL 2                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

nextflow.enable.dsl=2

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                                 USAGE                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info """
    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run c-sar --counts <count files> --library <library file> --info <sample mapping file>


    Mandatory arguments:
      --counts                   Path to count directory or count matrix file.
      --count_type               Must be one of: single or matrix [Default: single]
      --library                  Path to library file.
      --info                     Path to sample mapping file.

    Quality control:
      --no_qc                    Don't run QC.

    Filtering:
      --no_filter                Don't run count filtering.
      --filter_type              Filter type must be one of: all, plasmid, control or treatment [Default: all].
      --filter_method            Filter method must be one of: all, or mean [Default: mean].
      --min_reads                Minimum read count [Default: 30].

    Normalisation
      --no_normalisation         Don't run normalisation.
      --normalisation_method     Normalisation method must be one of: none, crisprcleanr, mageck_median, mageck_total, mageck_control, bagel [Default: mageck_median].

    Log fold changes
      --lfc_method               Package used to generate log fold changes must be one of: none, crisprcleanr, bagel [Default: bagel].

    Correction:
      --no_correction            Don't run any correction packages.
      --no_crisprcleanr          Don't run CRISPRcleanR.

    Analysis:
      --no_analysis              Don't run any analysis packages.
      --no_mageck                Don't run MAGeCK.
      --no_bagel                 Don't BAGEL.

    Optional arguments:
      --name                     Name for the pipeline run.
      --outdir                   The output directory where the results will be saved.
      --help                     Show pipeline usage.

    For software-specific parameters, please see the relevant configuration files in ./config.
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            VALIDATE INPUTS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Show help message
if ( params.help ) {
	helpMessage()
	exit 0
}

// Has the run name been specified by the user?
// Catches both -name and --name
custom_runName = workflow.runName
if (params.name) {
	custom_runName = params.name
}
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)){
	custom_runName = workflow.runName
}

def printErr = System.err.&println

// Have counts been provided?
if ( !params.counts  ) {
	printErr("Counts must be provided.")
	exit 1
}

// Has library been provided?
if ( !params.library  ) {
	printErr("Library must be provided.")
	exit 1
}

// Has sample mapping been provided?
if ( !params.info  ) {
	printErr("Sample mapping must be provided.")
	exit 1
}

// Count type must be in array
def count_types = ['single', 'matrix']
if ( count_types.contains( params.count_type ) == false ) {
	printErr("Count type can only be one of: " + count_types.join(',') + ".")
	exit 1
}

// Filter type must be in array
def filter_types = ['all', 'plasmid','control','treatment']
if ( filter_types.contains( params.filter_type ) == false ) {
	printErr("Filter type can only be one of: " + filter_types.join(',') + ".")
	exit 1
}

// Filter method must be in array
def filter_methods = ['all', 'any', 'mean', 'median']
if ( filter_methods.contains( params.filter_method ) == false ) {
	printErr("Filter methods can only be one of: " + filter_methods.join(',') + ".")
	exit 1
}

// Normalisation method must be in array
def normalisation_method = ['none', 'crisprcleanr', 'mageck_median', 'mageck_total', 'bagel']
if ( normalisation_method.contains( params.normalisation_method ) == false ) {
	printErr("Normalisation method can only be one of: " + normalisation_method.join(',') + ".")
	exit 1
}

// LFC method must be in array
def lfc_method = ['crisprcleanr']
if ( lfc_method.contains( params.lfc_method ) == false ) {
	printErr("Log fold change method can only be one of: " + lfc_method.join(',') + ".")
	exit 1
}

// Cannot perform scaling if essential and non-essential gene lists are not provided
// @TODO: check this is working as expected
if ( params.essential_genes == '' || !params.essential_genes || params.nonessential_genes == '' || !params.nonessential_genes  ) {
  if ( params.scale_log_fold_changes || params.scale_bayes_factors ) {
    printErr("Cannot perform scaling when no essential or non-essential gene files are supplied.")
	  exit 1
  }
  if ( !params.no_BAGEL & !params.no_analysis ) {
    printErr("Cannot run BAGEL when no essential or non-essential gene files are supplied.")
	  exit 1
  }
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      Parse input files                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Gather input data into channels (counts, library and sample mapping)
Channel
  .fromPath( params.counts, checkIfExists: true )
  .set { counts_directory }

Channel
  .fromPath( params.library, checkIfExists: true )
  .set { library }

Channel
  .fromPath( params.info, checkIfExists: true )
  .set { sample_mapping }

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               MODULES                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Load module dependencies

// crispr_utils
include { counts2matrix;
          read_count_matrix;
          remove_user_defined_guides;
          remove_duplicate_guides;
          filter_counts_by_indices;
          remove_guides_with_no_coordinates;
          calculate_log_fold_changes;
          scale_gene_log_fold_changes } from './modules/crispr_utils/common' params(params)

// crispr_qc
include { sequencing_qc;
          intermediate_qc as intermediate_qc_gene_fc;
          intermediate_qc as intermediate_qc_gene_counts;
          intermediate_qc as intermediate_qc_sgrna_counts;
          intermediate_bagel_classification_qc as bagel_classification_qc_gene_fc;
          intermediate_bagel_classification_qc as bagel_classification_qc_gene_counts;
          intermediate_bagel_classification_qc as bagel_classification_sgrna_counts; } from './modules/crispr_qc/common' params(params)

// CRISPRcleanR
include { crisprcleanr_normalise_counts;
          format_library_and_matrices_for_crisprcleanr;
          crisprcleanr_correction } from './modules/CRISPRcleanR/common' params(params)

// MAGeCK
include { MAGeCK_normalisation;
          MAGeCK_test;
          MAGeCK_process_results } from './modules/MAGeCK/common' params(params)

// BAGEL2
include { bagel_normalise_counts;
          BAGEL_bf as BAGEL_bf_sgrna;
          BAGEL_bf as BAGEL_bf_gene;
          scale_gene_BFs  } from './modules/BAGEL2/common' params(params)

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      Get sample type counts                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Start by setting number of each sample type (plasmid, control and treatment to 0)
int num_plasmid_samples = 0
int num_control_samples = 0
int num_treatment_samples = 0
int total_num_samples = 0

// Read sample mapping file line by line to get number of samples of each type
// Has to be done outside of the main workflow
String line
int line_count = 0
myReader = file("${params.info}").newReader()
while( line = myReader.readLine() ) {
    if ( params.info_header && line_count == 0 ) {
        // do nothing with first line if sample mapping has no header
    } else {
        num_plasmid_samples += Integer.parseInt( line.tokenize(params.info_delim)[ Integer.parseInt("${params.info_plasmid_column_index}") - 1 ] )
        num_control_samples += Integer.parseInt( line.tokenize(params.info_delim)[ Integer.parseInt("${params.info_control_column_index}") - 1 ] )
        num_treatment_samples += Integer.parseInt( line.tokenize(params.info_delim)[ Integer.parseInt("${params.info_treatment_column_index}") - 1 ] )
        total_num_samples += 1
    }
    line_count++
}
myReader.close()

// Check that the total number of samples matches the sum of the sample type counts
if ( total_num_samples != ( num_plasmid_samples + num_control_samples + num_treatment_samples ) ) {
  printErr("Total number of samples does not match sample type count. Check samples are each only marked as one of plasmid, control and treatments")
	exit 1
}

// Set the sample totals as a list
sample_totals = [ 'plasmid': num_plasmid_samples,
                  'control': num_control_samples,
                  'treatment': num_treatment_samples ]

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   Calculalating analysis indices                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// // Use number of samples per sample type to generate indices for analyses
def get_analysis_indices( sample_totals ) {
  // pull out number of samples per sample type
  n_plasmid = sample_totals['plasmid']
  n_control = sample_totals['control']
  n_treatment = sample_totals['treatment']

  def analysis_indices = [:]
  analysis_indices << [ 'total': [  'all': ( n_plasmid + n_control + n_treatment ),
                                    'plasmid': n_plasmid,
                                    'control': n_control,
                                    'treatment': n_treatment ] ]

  analysis_indices << [ 'all': [  'base0': ( 0..( n_plasmid + n_control + n_treatment - 1 ) ).join( ',' ),
                                  'base1': ( 1..( n_plasmid + n_control + n_treatment ) ).join( ',' ),
                                  'base1_increment2': ( 3..( 2 + n_plasmid + n_control + n_treatment ) ).join( ',' ) ] ]

  if ( n_plasmid > 0 ) {
    analysis_indices << [ 'plasmid': [  'base0': ( 0..( n_plasmid - 1 ) ).join( ',' ),
                                        'base1': ( 1..( n_plasmid ) ).join( ',' ),
                                        'base1_increment2': ( 3..( 2 + n_plasmid ) ).join( ',' ) ] ]
  }
  if ( n_control > 0 ) {
    if ( n_plasmid == 0 ) {
      analysis_indices << [ 'control': [  'base0': ( 0..( n_control - 1 ) ).join( ',' ),
                                          'base1': ( 1..n_control ).join( ',' ),
                                          'base1_increment2': ( 3..(2 + n_control) ).join( ',' ) ] ]
    } else {
      analysis_indices << [ 'control': [  'base0': ( ( n_plasmid )..( n_plasmid + n_control - 1 ) ).join( ',' ),
                                          'base1': ( ( n_plasmid + 1 )..( n_plasmid + n_control ) ).join( ',' ),
                                          'base1_increment2': ( ( 2 + n_plasmid + 1 )..( 2 + n_plasmid + n_control ) ).join( ',' ) ] ]
    }
  }
  if ( n_treatment > 0 ) {
    if ( ( n_plasmid + n_control ) == 0 ) {
      analysis_indices << [ 'treatment': [  'base0': ( 0..( ( n_plasmid + n_control ) - 1 ) ).join( ',' ),
                                            'base1': ( 1..( n_plasmid + n_control ) ).join( ',' ),
                                            'base1_increment2': ( 3..( 2 + n_plasmid + n_control ) ).join( ',' ) ] ]
    } else {
      analysis_indices << [ 'treatment': [  'base0': ( ( n_plasmid + n_control )..( n_plasmid + n_control + n_treatment - 1 ) ).join( ',' ),
                                            'base1': ( ( n_plasmid + n_control + 1 )..( n_plasmid + n_control + n_treatment ) ).join( ',' ),
                                            'base1_increment2': ( ( 2 + n_plasmid + n_control + 1 )..( 2 + n_plasmid + n_control + n_treatment ) ).join( ',' ) ] ]
    }
  }

  n_contrasts = 0
  if ( n_plasmid > 0 && n_control > 0 ) {
    n_contrasts++
    analysis_indices << [ 'control_vs_plasmid': [ 'plasmid':        [   'base0': ( 0..( n_plasmid - 1 ) ).join( ',' ),
                                                                        'base1': ( 1..( n_plasmid ) ).join( ',' ) ],
                                                  'control':        [   'base0': ( ( n_plasmid )..( ( n_plasmid + n_control ) - 1 ) ).join( ',' ),
                                                                        'base1': ( ( n_plasmid + 1 )..( n_plasmid + n_control ) ).join( ',' ) ],
                                                  'count_lfc':      [   'base1_increment2': ( 3..( 2 + n_plasmid + n_control ) ).join( ',' ) ],
                                                  'lfc':            [   'base0': ( 0..( n_control - 1 ) ).join( ',' ),
                                                                        'base1': ( 1..n_control ).join( ',' ),
                                                                        'base1_increment1': ( 2..(1 + n_control) ).join( ',' ),
                                                                        'base1_increment2': ( 3..(2 + n_control) ).join( ',' ),
                                                                        'base1_scaled_lfc': ( 2..(n_control - 1) ).join( ',' ) ]
                                                ] ]
  }
  if ( n_plasmid > 0 && n_treatment > 0 ) {
    n_contrasts++
    analysis_indices << [ 'treatment_vs_plasmid': [ 'plasmid':          [   'base0': ( 0..( n_plasmid - 1 ) ).join( ',' ),
                                                                            'base1': ( 1..( n_plasmid ) ).join( ',' ) ],
                                                    'treatment':        [   'base0': ( ( n_plasmid )..( ( n_plasmid + n_treatment ) - 1 ) ).join( ',' ),
                                                                            'base1': ( ( n_plasmid + 1 )..( n_plasmid + n_treatment ) ).join( ',' ) ],
                                                    'count_lfc':        [   'base1_increment2': ( 3..( 2 + n_plasmid + n_treatment ) ).join( ',' ) ],
                                                    'lfc':              [   'base0': ( 0..( n_treatment - 1 ) ).join( ',' ),
                                                                            'base1': ( 1..n_treatment ).join( ',' ),
                                                                            'base1_increment2': ( 3..(2 + n_treatment) ).join( ',' ),
                                                                            'base1_increment1': ( 2..(1 + n_treatment) ).join( ',' ),
                                                                            'base1_scaled_lfc': ( 2..(n_treatment - 1) ).join( ',' ) ] ] ]
  }
  if ( n_control > 0 && n_treatment > 0 ) {
    n_contrasts++
    analysis_indices << [ 'treatment_vs_control': [ 'control':          [   'base0': ( 0..( n_control - 1 ) ).join( ',' ),
                                                                            'base1': ( 1..( n_control ) ).join( ',' ) ],
                                                    'treatment':        [   'base0': ( ( n_control )..( ( n_control + n_treatment ) - 1 ) ).join( ',' ),
                                                                            'base1': ( ( n_control + 1 )..( n_control + n_treatment ) ).join( ',' ) ],
                                                    'count_lfc':        [   'base1_increment2': ( 3..( 2 + n_control + n_treatment ) ).join( ',' ) ],
                                                    'lfc':              [   'base0': ( 0..( n_treatment - 1 ) ).join( ',' ),
                                                                            'base1': ( 1..n_treatment ).join( ',' ),
                                                                            'base1_increment2': ( 3..(2 + n_treatment) ).join( ',' ),
                                                                            'base1_increment1': ( 2..(1 + n_treatment) ).join( ',' ),
                                                                            'base1_scaled_lfc': ( 2..(n_treatment - 1) ).join( ',' ) ] ] ]
  }
  analysis_indices <<  [ 'n_contrasts': n_contrasts ]
  return analysis_indices
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      Validate sample types                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Cannot filter when there are no samples of the user-defined filter type
if ( num_plasmid_samples == 0 && params.filter_type == 'plasmid' ) {
	printErr("Filter type cannot be plasmid when there are no plasmid samples")
	exit 1
}
if ( num_control_samples == 0 && params.filter_type == 'control' ) {
	printErr("Filter type cannot be control when there are no control samples")
	exit 1
}
if ( num_treatment_samples == 0 && params.filter_type == 'treatment' ) {
	printErr("Filter type cannot be treatment when there are no treatment samples")
	exit 1
}

// Cannot calculate log fold changes or run analyses if there are no plasmid or control samples
if ( num_control_samples == 0 && num_plasmid_samples == 0 ) {
  // Cannot normalise with CRISPRcleanR as normalisation and LFC calculations are in same function
  if ( params.normalisation_method == 'crisprcleanr' ) {
    printErr("CRISPRcleanR cannot be used for normalisation when there are no plasmid or control samples.")
	  exit 1
  }
  // Cannot run correction
  if ( !params.no_correction ) {
    printErr("Cannot perform correction when there are no plasmid or control samples.")
	  exit 1
  }
  //Cannot run analyses
  if ( !params.no_analysis ) {
    printErr("Cannot run analyses when there are no plasmid or control samples.")
	  exit 1
  }
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         SUMMARY INFO                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Runtime information to print at the top of the log file
log.info "====================================="
log.info "${workflow.manifest.name} pipeline - version ${workflow.manifest.version}"
log.info "WORKFLOW PARAMETERS"
log.info "====================================="

def summary = [:]

// Main run parameters
summary['Run name']           = custom_runName
summary['Counts']             = params.counts
summary['Count type']         = params.count_type
summary['Library']            = params.library
summary['Sample mapping']     = params.info
summary['Results directory']  = params.outdir

// Quality control
summary['Quality control']    = !params.no_qc

// Log fold changes
summary['LFC method']         = params.lfc_method

// Duplicate removal
summary['Duplicate removal']  = !params.no_duplicate_removal

// Filtering
summary['Filtering']          = !params.no_filter
if ( !params.no_filter ) {
  summary['Filter type']      = params.filter_type
  summary['Filter method']    = params.filter_method
  summary['Minimum reads']    = params.min_reads
}

// Normalisation
if ( params.no_normalisation ) {
  summary['Normalisation'] = false
} else {
  if ( params.normalisation_method == 'none'  ) {
    summary['Normalisation'] = false
  } else {
    summary['Normalisation']         = true
    summary['Normalisation method']  = params.normalisation_method
  }
}

// Correction
if ( params.no_correction ) {
  summary['Correction']    = false
} else {
  if ( params.no_crisprcleanr ) {
    summary['Correction']  = false
    summary['CRISPRcleanR']  = false
  } else {
    summary['Correction']  = true
    summary['CRISPRcleanR']  = true
  }
}

// Gene lists
summary['Essential genes']      = params.essential_genes
summary['Non-essential genes']  = params.nonessential_genes

// Scaling of log fold changes
summary['Scale LFCs']      = params.scale_log_fold_changes
summary['Scale BFs']      = params.scale_bayes_factors

// Analysis flows
if ( params.no_analysis ) {
  summary['Analysis']  = false
  summary['MAGeCK']    = false
  summary['BAGEL2']    = false
} else {
  summary['Analysis']  = true

  // MAGeCK
  summary['MAGeCK']        = !params.no_mageck
  if ( !params.no_mageck ) {
    summary['MAGeCK normalised counts to file']         = params.mageck_normcounts_to_file
    summary['MAGeCK remove zero/low counts']            = params.mageck_remove_zero
    summary['MAGeCK remove zero/low count threshold ']  = params.mageck_remove_zero_threshold
    summary['MAGeCK extra options']                     = params.mageck_extra_options
  }

  // BAGEL2
  summary['BAGEL2']        = !params.no_bagel
  if ( !params.no_bagel ) {
    summary['BAGEL2 normalisation pseudocount']         = params.bagel_normalisation_pseudocount
    summary['BAGEL2 (fc) extra options']                = params.bagel_fc_extra_options
    summary['BAGEL2 (bf) extra options']                = params.bagel_bf_extra_options
  }
}

log.info summary.collect { k,v -> "${k.padRight(40)}: $v" }.join("\n")
log.info "====================================="

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        REUSABLE WORKFLOWS                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process write_pipeline_version {
    publishDir "${params.resultDir}", mode: 'copy', pattern: "c-sar.version"
    output:
      path "c-sar.version"
    script:
      """
      echo "${workflow.manifest.version}" > c-sar.version
      """
}

workflow prepare_contrasts {
  take:
    count_matrices
    analysis_indices

  main:
    control_vs_plasmid =Channel.empty()
    treatment_vs_plasmid = Channel.empty()
    treatment_vs_control = Channel.empty()
    analysis_list = Channel.empty()

    if ( analysis_indices['total']['plasmid'] > 0 && analysis_indices['total']['control'] > 0 ) {
      control_vs_plasmid = count_matrices.combine( ['control_vs_plasmid'] )
      analysis_list = analysis_list.mix( control_vs_plasmid )
    }
    if ( analysis_indices['total']['plasmid'] > 0 && analysis_indices['total']['treatment'] > 0 ) {
      treatment_vs_plasmid = count_matrices.combine( ['treatment_vs_plasmid'] )
      analysis_list = analysis_list.mix( treatment_vs_plasmid )
    }
    if ( analysis_indices['total']['control'] > 0 && analysis_indices['total']['treatment'] > 0 ) {
      treatment_vs_control = count_matrices.combine( ['treatment_vs_control'] )
      analysis_list = analysis_list.mix( treatment_vs_control )
    }

  emit:
    analysis_list
}

workflow  intermediate_qc_per_stage {
    take:
      matrices
      sample_mapping
      analysis_indices

    main:
      intermediate_qc_gene_fc( matrices.combine(sample_mapping), analysis_indices, 'fc' )
      intermediate_qc_sgrna_counts( matrices.combine(sample_mapping), analysis_indices, 'counts' )
      bagel_classification_qc_gene_fc( matrices.combine(sample_mapping), analysis_indices, 'fc' )
      bagel_classification_sgrna_counts( matrices.combine(sample_mapping), analysis_indices, 'counts' )
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    DUPLICATE REMOVAL WORKFLOW                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Workflow for filtering raw counts
workflow remove_duplicates {
  take:
    count_matrix
    library

  main:
    duplicate_removed_counts = Channel.empty()
    duplicate_removed_library = Channel.empty()
    duplicate_guides = Channel.empty()

    if ( !params.no_duplicate_removal ) {
      ( duplicate_removed_counts, duplicate_removed_library, duplicate_guides ) = remove_duplicate_guides( count_matrix, library )
    } else {
      duplicate_removed_counts = count_matrix
      duplicate_removed_library = library
    }

  emit:
    duplicate_removed_counts
    duplicate_removed_library
    duplicate_guides
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       FILTERING WORKFLOW                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Workflow for filtering raw counts
workflow raw_filter {
  take:
    count_matrix
    analysis_indices
    filter_type
    total_samples

  main:
    filtered_counts = Channel.empty()
    excluded_guides = Channel.empty()

    if ( params.filter_counts ) {
      ( filtered_counts, filtered_library, excluded_guides) = filter_counts_by_indices( count_matrix, analysis_indices["all"]["base1_increment2"], 'raw_filter', total_samples )
    } else {
      filtered_counts = count_matrix
    }

  emit:
    filtered_counts
    excluded_guides
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    NORMALISATION WORKFLOWS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Workflow for normalisation
workflow normalisation {
  take:
    count_matrix
    normalisation_method
    analysis_indices

  main:
    normalised_counts = Channel.empty()

    if ( params.normalisation_method == 'crisprcleanr' ) {
      if ( num_plasmid_samples == 0 && num_control_samples == 0 ) {
        // Cannot run CRISPRcleanR with only one sample type as it generates fold changes
        printErr("CRISPRcleanR cannot be used for normalisation when there are no plasmid or control samples.")
	      exit 1
      } else if ( num_plasmid_samples >= 1 ) {
        // R wrapper for CRISPRcleanR methods using plasmid samples count
        crisprcleanr_normalise_counts( count_matrix, analysis_indices['total']['plasmid'], analysis_indices["all"]['base1_increment2'] )
        normalised_counts = crisprcleanr_normalise_counts.out.normalised_count_matrix
      } else {
        // R wrapper for CRISPRcleanR methods using control sample count
        crisprcleanr_normalise_counts( count_matrix, analysis_indices['total']['control'], analysis_indices["all"]['base1_increment2'] )
        normalised_counts = crisprcleanr_normalise_counts.out.normalised_count_matrix
      }
    } else if ( params.normalisation_method =~ /^mageck/ ) {
      // Get MAGeCK normalisation method from parameter string
      mageck_norm_method = params.normalisation_method.split('_')[1]
      normalised_counts = MAGeCK_normalisation( count_matrix, mageck_norm_method, analysis_indices['treatment']['base0'] )
    } else if ( params.normalisation_method == 'bagel' ) {
      // Rscript which has been derived from BAGEL2 normalisation methods
      bagel_normalise_counts( count_matrix, analysis_indices["all"]["base1_increment2"] )
      normalised_counts = bagel_normalise_counts.out.normalised_count_matrix
    }
  emit:
    normalised_counts = normalised_counts
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       CORRECTION WORKFLOWS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

workflow CRISPRcleanR {
  take:
    all_matrices
    analysis_indices

  main:
    if ( analysis_indices["total"]["plasmid"] == 0 && analysis_indices["total"]["control"] == 0 ) {
      // Cannot run CRISPRcleanR with only one sample type
      printErr("CRISPRcleanR cannot be used for normalisation when there are no plasmid or control samples.")
      exit 1
    }
    // We only want to run this on the output files from the last stage
    matrices_for_correction = all_matrices.buffer( size: analysis_indices["n_contrasts"] ).last().flatMap{n -> n}

    // Format the library, count matrix and fold change matrix for CRISPRcleanR
    formatted_input = format_library_and_matrices_for_crisprcleanr( matrices_for_correction, analysis_indices )

    // Run CRISPRcleanR
    crisprcleanr_sgrna_fold_changes = crisprcleanr_correction( formatted_input.inputdata, 'corrected', analysis_indices )

  emit:
    crisprcleanr_sgrna_fold_changes
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      ANALYSIS WORKFLOWS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Workflow for running MAGeCK
workflow MAGeCK {
  take:
    all_matrices
    analysis_indices

  main:
    summaries_sgrna = Channel.empty()
    summaries_gene = Channel.empty()

    // We only want to run this on the output files from the last stage
    matrices_for_correction = all_matrices.buffer( size: analysis_indices['n_contrasts'] ).last().flatMap{n -> n}

    // Run MAGeCK test
    (sgrna_summary, gene_summary ) = MAGeCK_test( matrices_for_correction, 'none', analysis_indices )

    // Process MAGeCK results
    MAGeCK_process_results( gene_summary, sgrna_summary )

  emit:
    sgrna_summary = summaries_sgrna
    gene_summary = summaries_gene
}

// Workflow for running BAGEL
workflow BAGEL_bf {

  take:
    all_matrices
    analysis_indices

  main:
    bagel_gene = Channel.empty()
    bagel_sgrna = Channel.empty()

    // We only want to run this on the output files from the last stage
    matrices_for_correction = all_matrices.buffer( size: analysis_indices['n_contrasts'] ).last().flatMap{n -> n}

    // Run BAGEL bf
    bagel_gene = BAGEL_bf_gene( matrices_for_correction, analysis_indices, 'gene' )
    bagel_sgrna = BAGEL_bf_sgrna( matrices_for_correction, analysis_indices, 'sgrna' )

    // Scale BAGEL gene BFs
    scale_gene_BFs( bagel_sgrna, analysis_indices)

  emit:
    bagel_gene
    bagel_sgrna
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           MAIN WORKFLOW                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Main runner for pipeline
workflow {
  main:
    // Write pipeline version
    write_pipeline_version()

    // Create empty channels
    library = library.map { file -> tuple( 'raw', file ) }
    count_matrices = Channel.empty()
    guides_removed_in_pipeline = Channel.empty()

    // Get column indices (0-based and 1-based) of sample types for count matrices
    analysis_indices = get_analysis_indices( sample_totals )

    if (params.count_type == "single") {
        // Prepare count matrix from individual sample counts
        count_matrices = count_matrices.mix( counts2matrix( counts_directory, library, sample_mapping ).data )
    } else {
        // Prepare count matrix from individual sample counts
        count_matrices = count_matrices.mix( read_count_matrix( counts_directory, library, sample_mapping ).data )
    }

    // Sequencing QC (only run on raw counts)
    sequencing_qc( count_matrices.filter( ~/.*raw.*/ ), sample_mapping, analysis_indices["all"]["base1_increment2"])

    // Remove unwanted guides
    ( user_removed_counts, user_removed_guides ) = remove_user_defined_guides( count_matrices.last(), analysis_indices["all"]["base1_increment2"] )
    count_matrices = count_matrices.mix( user_removed_counts )
    guides_removed_in_pipeline = guides_removed_in_pipeline.mix( user_removed_guides )

    // Remove duplicate guide IDs from library and counts
    ( duplicate_removed_counts, duplicate_removed_guides ) = remove_duplicate_guides( count_matrices.last(), analysis_indices["all"]["base1_increment2"] )
    count_matrices = count_matrices.mix( duplicate_removed_counts )
    guides_removed_in_pipeline = guides_removed_in_pipeline.mix( duplicate_removed_guides )

    // Filter counts
    ( filtered_counts, raw_filter_removed_guides ) = raw_filter( count_matrices.last(), analysis_indices, params.filter_type, analysis_indices["all"]["base1_increment2"] )
    count_matrices = count_matrices.mix( filtered_counts )
    guides_removed_in_pipeline = guides_removed_in_pipeline.mix( raw_filter_removed_guides )

    // Normalise counts
    normalised_counts = normalisation( count_matrices.last(), params.normalisation_method, analysis_indices )
    count_matrices = count_matrices.mix( normalised_counts )

    // Remove guides without coordinates (only necessary if correction is being used )
    if ( !params.no_correction && !params.no_crisprcleanr ) {
      ( coord_filtered_count_matrix, coord_filter_removed_guides ) = remove_guides_with_no_coordinates( count_matrices.last(), 'coord_filter', analysis_indices["all"]["base1_increment2"] )
      count_matrices = count_matrices.mix( coord_filtered_count_matrix )
      guides_removed_in_pipeline = guides_removed_in_pipeline.mix( coord_filter_removed_guides )
    }

    // Split counts out by contrast for correction and analyses
    contrast_count_matrices = prepare_contrasts( count_matrices, analysis_indices )

    // Calculate fold changes by stage and by contrast
    contrast_combined_count_and_fold_change_matrices = calculate_log_fold_changes( contrast_count_matrices, analysis_indices )

    // CN correction
    corrected_contrasts = CRISPRcleanR( contrast_combined_count_and_fold_change_matrices, analysis_indices )
    contrast_combined_count_and_fold_change_matrices = contrast_combined_count_and_fold_change_matrices.mix( corrected_contrasts )

    // Intermediate QC
    intermediate_qc_per_stage( contrast_combined_count_and_fold_change_matrices, sample_mapping, analysis_indices )

    // Scale LFCs
    scaled_gene_lfc = scale_gene_log_fold_changes( contrast_combined_count_and_fold_change_matrices, analysis_indices )

    // Run MAGeCK
    MAGeCK( contrast_combined_count_and_fold_change_matrices, analysis_indices )

    // Run BAGEL
    BAGEL_bf( contrast_combined_count_and_fold_change_matrices, analysis_indices )
}
