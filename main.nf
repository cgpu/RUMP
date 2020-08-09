#!/usr/bin/env nextflow

/**
    RUMP: A Reproducible Untargeted Metabolomics Data Processing Pipeline
    Description : A Nextflow-based reproducible pipeline for untargeted metabolomics data analysis
    Copyright   : (C) LemasLab
    Author      : Xinsong Du
    License     : GNU GPL-v3.0 License
          
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this script. If not, see <http://www.gnu.org/licenses/>.
    
    For any bugs or problems found, please contact us at
    - xinsongdu@ufl.edu, manfiol@ufl.edu, djlemas@ufl.edu; 
    - https://github.com/lemaslab/RUMP
*/

// Those variable names which are all uppercase are channel names

version='0.0.0'
timestamp='20200226'

MZMINE = Channel.fromPath(params.mzmine_dir, type: 'dir') // The location of folder of MzMine
MZMINE.into{POS_MZMINE; NEG_MZMINE} // Duplicate the MZMINE chennel into two channels, one of which deals with positive sample while the other deals with negative sample.

POS_DATA_DIR = Channel.fromPath(params.input_dir_pos, type: 'dir') // Location of folder storing positive data
POS_DATA_DIR.into{POS_DATA_DIR_UNIT_TESTS; POS_DATA_DIR_INFO; POS_DATA_DIR_BS}
NEG_DATA_DIR = Channel.fromPath(params.input_dir_neg, type: 'dir') // Location of folder storing negative data
NEG_DATA_DIR.into{NEG_DATA_DIR_UNIT_TESTS; NEG_DATA_DIR_INFO; NEG_DATA_DIR_BS}

// Design files for positive data and negative data.
POS_DESIGN = Channel.fromPath(params.POS_design_path)
POS_DESIGN.into{POS_DESIGN_FOR_UNIT_TESTS; POS_DESIGN_FOR_AS; POS_DESIGN_FOR_BS; POS_DESIGN_FOR_PCA_NOBG; POS_DESIGN_FOR_PCA_WITHBG; POS_DESIGN_FOR_HCLUSTERING_NOBG; POS_DESIGN_FOR_HCLUSTERING_WITHBG; POS_DESIGN_FOR_VD_NOBG; POS_DESIGN_FOR_VD_WITHBG; POS_DESIGN_FOR_BARPLOT_NOBG; POS_DESIGN_FOR_BARPLOT_WITHBG}
NEG_DESIGN = Channel.fromPath(params.NEG_design_path)
NEG_DESIGN.into{NEG_DESIGN_FOR_UNIT_TESTS; NEG_DESIGN_FOR_AS; NEG_DESIGN_FOR_BS; NEG_DESIGN_FOR_PCA_NOBG; NEG_DESIGN_FOR_PCA_WITHBG; NEG_DESIGN_FOR_HCLUSTERING_NOBG; NEG_DESIGN_FOR_HCLUSTERING_WITHBG; NEG_DESIGN_FOR_VD_NOBG; NEG_DESIGN_FOR_VD_WITHBG; NEG_DESIGN_FOR_BARPLOT_NOBG; NEG_DESIGN_FOR_BARPLOT_WITHBG}

// Library
POS_LIBRARY = Channel.fromPath(params.pos_library)
NEG_LIBRARY = Channel.fromPath(params.neg_library)
POS_LIBRARY.into{POS_LIBRARY_MZMINE; POS_LIBRARY_STAT}
NEG_LIBRARY.into{NEG_LIBRARY_MZMINE; NEG_LIBRARY_STAT}

// Pre-build MultiQC report information
EXPERIMENTS_INFO = Channel.fromPath(params.experiments_info)
MQC_CONFIG = Channel.fromPath(params.mqc_config)

// Result files used by MultiQC to generate report.
// MQC_DIR = Channel.fromPath(params.mqc_dir, type: 'dir')

/**
    Prints version when asked for
*/

if (params.version) {
    System.out.println("")
    System.out.println("RUMP: A Reproducible Untargeted Metabolomics Data Processing Pipeline - Version: $version ($timestamp)")
    exit 1
}

/**
    Basic running information
*/

println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Manifest's pipeline version: $workflow.manifest.version"

/**
    Prints help when asked for
*/

if (params.help) {
    System.out.println("")
    System.out.println("RUMP: A Reproducible Untargeted Metabolomics Data Processing Pipeline - Version: $version ($timestamp)")
    System.out.println("This pipeline is distributed in the hope that it will be useful")
    System.out.println("but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.")
    System.out.println("")
    System.out.println("Please report comments and bugs to xinsongdu@ufl.edu")
    System.out.println("or at https://github.com/lemaslab/RUMP/issues.")
    System.out.println("Check https://github.com/lemaslab/RUMP for updates, and refer to")
    System.out.println("https://github.com/lemaslab/RUMP/wiki")
    System.out.println("")
    System.out.println("Usage:  ")
    System.out.println("   nextflow run_all.nf [options] -with-docker xinsongdu/lemaslab_rump:v0.0.0")
    System.out.println("")
    System.out.println("Arguments (it is mandatory to change `input_file` and `mzmine_dir` before running:")
    System.out.println("----------------------------- common parameters ----------------------------------")
    System.out.println("    --input_dir_pos                         folder location for positive data, default is 'data/POS'")
    System.out.println("    --input_dir_neg                         folder location for positive data, default is 'data/NEG'")
    System.out.println("    --POS_design_path                       location for positive design file, default is 'data/pos_design.csv'")
    System.out.println("    --NEG_design_path                       location for negative design file, default is 'data/neg_design.csv'")
    System.out.println("    --cutoff                                cutoff p-value for mummichog pathway analysis, default is 0.05")
    System.out.println("    --unknown_search                        whether do unknown search for unidentified metabolites or not, default is '1', please set it to '0' when you want to disable it")
    System.out.println("    --version                               whether to show version information or not, default is null")
    System.out.println("    --help                                  whether to show help information or not, default is null")
    System.out.println("Please refer to nextflow.config for more options.")
    System.out.println("")
    System.out.println("Container:")
    System.out.println("    Docker image to use with -with-docker|-with-singularity options is")
    System.out.println("    'docker://xinsongdu/lemaslab_rump:v0.0.0'")
    System.out.println("")
    System.out.println("RUMP supports .mzXML format files.")
    System.out.println("")
    exit 1
}

// Check appropriateness of input
process input_check {

    echo true

    input:
    file pos_data_dir from POS_DATA_DIR_UNIT_TESTS // Location of positive data
    file neg_data_dir from NEG_DATA_DIR_UNIT_TESTS // Location of negative data
    file pos_design from POS_DESIGN_FOR_UNIT_TESTS // Location of positive design
    file neg_design from NEG_DESIGN_FOR_UNIT_TESTS // Location of negative design

    script:
    """
    echo "checking input eligibility" &&
    input_check.py --pos_data ${pos_data_dir} --neg_data ${neg_data_dir} --pos_design ${pos_design} --neg_design ${neg_design}
    """
}

// Process for generating MultiQC report regarding data information
process mqc_data_info {

    publishDir "${params.outdir}/MultiQC/", mode: 'copy' //copy the output files to the folder "${params.outdir}/mqc"

    input:
    file pos_data_dir from POS_DATA_DIR_INFO // Location of positive data
    file neg_data_dir from NEG_DATA_DIR_INFO // Location of negative data

    // POS_DATA and NEG_DATA are channels containing filtered POS and NEG data, which are ready to be input to R codes.
    output:
    file params.pos_data_info_mqc into POS_DATA_INFO_MQC // file regarding positive data information that can be parsed by MultiQC
    file params.neg_data_info_mqc into NEG_DATA_INFO_MQC // file regarding negative data information that can be parsed by MultiQC

    script:
    """
    sleep 5 &&
    data_info.py -i ${pos_data_dir} -o $params.pos_data_info_mqc -n p &&
    data_info.py -i ${neg_data_dir} -o $params.neg_data_info_mqc -n n
    """
}

process batchfile_generation_mzmine {

    echo true

    input:
    file pos_data_dir from POS_DATA_DIR_BS // Location of positive data
    file neg_data_dir from NEG_DATA_DIR_BS // Location of negative data

    output:
    file params.pos_config into POS_BATCHFILE // Generated batchfile for processing positive data
    file params.neg_config into NEG_BATCHFILE // Generated batchfile for processing negative data

    script:
    """ 
    sleep 15 && 
    echo "setting parameters for MZmine" &&
    batchfile_generator_pos_253.py -x ${params.pos_config} -i ${pos_data_dir} -l $params.pos_library -o $params.pos_mzmine_peak_output &&
    batchfile_generator_neg_253.py -x ${params.neg_config} -i ${neg_data_dir} -l $params.neg_library -o $params.neg_mzmine_peak_output
    """
}

process pos_peakDetection_mzmine {

    publishDir "${params.outdir}/peakDetection_mzmine/pos/", mode: 'copy'
    echo true

    input:
    file p_b from POS_BATCHFILE // Batchfile for MzMine to process positive data.
    file pos_library from POS_LIBRARY_MZMINE // Location of library file for positive samples

    output:
    file "${params.pos_mzmine_peak_output}" optional true // MzMine processing result for positive data.

// Change "startMZmine_Linux.sh" to "startMZmine_MacOSX.command" in the following code if running locally with Mac

    script:
    """
    wget https://github.com/mzmine/mzmine2/releases/download/v2.53/MZmine-2.53-Linux.zip && \
    unzip MZmine-2.53-Linux.zip && \
    rm MZmine-2.53-Linux.zip && \
    mv MZmine-2.53-Linux/* .
    ls -l
    echo "peak detection and library matching for positive data" &&
    ./startMZmine-Linux  ${p_b}
    ls -l
    """
}

// process neg_peakDetection_mzmine {

//     publishDir "${params.outdir}/peakDetection_mzmine/neg/", mode: 'copy'

//     input:
//     file n_b from NEG_BATCHFILE // Batchfile for MzMine to process negative data.
//     file neg_library from NEG_LIBRARY_MZMINE // Location of library file for negative samples (currently still use the positive library)
//     file n_m from NEG_MZMINE // Folder of MzMine tool

//     output:
//     file "${n_m}/${params.neg_mzmine_peak_output}" into NEG_MZMINE_RESULT // MzMine processing result for negative data.
//     stdout result

//     // Change "startMZmine_Linux.sh" to "startMZmine_MacOSX.command" in the following code if running locally with Mac

//     script:
//     """ 
//     sleep 5 &&  
//     echo "peak detection and library matching for negative data" &&
//     mv ${n_b} ${n_m}  && mv ${neg_library} ${n_m}  && cd ${n_m}  && ./startMZmine-Linux ${n_b}
//     """
// }

// process add_stats {

//     publishDir "${params.outdir}/peak_table/", mode: 'copy'
//     echo true

//     input:
//     file data_pos from POS_MZMINE_RESULT
//     file pos_design from POS_DESIGN_FOR_AS
//     file data_neg from NEG_MZMINE_RESULT
//     file neg_design from NEG_DESIGN_FOR_AS
//     file pos_library from POS_LIBRARY_STAT
//     file neg_library from NEG_LIBRARY_STAT

//     output:
//     file params.pos_data_nobg into POS_DATA_NOBG
//     file params.neg_data_nobg into NEG_DATA_NOBG

//     script:
//     """   
//     add_stats.py -i ${data_pos} -d ${pos_design} -o ${params.pos_data_nobg} -l ${pos_library} &&
//     add_stats.py -i ${data_neg} -d ${neg_design} -o ${params.neg_data_nobg} -l ${neg_library}
//     """
// }

// POS_DATA_NOBG.into{POS_NOBG_FOR_BS; POS_NOBG_FOR_MQC; POS_NOBG_FOR_PCA; POS_NOBG_FOR_HCLUSTERING; POS_NOBG_FOR_VD; POS_NOBG_FOR_BARPLOT; POS_NOBG_FOR_MUMMICHOG; POS_NOBG_FOR_UNKNOWN_SEARCH}
// NEG_DATA_NOBG.into{NEG_NOBG_FOR_BS; NEG_NOBG_FOR_MQC; NEG_NOBG_FOR_PCA; NEG_NOBG_FOR_HCLUSTERING; NEG_NOBG_FOR_VD; NEG_NOBG_FOR_BARPLOT; NEG_NOBG_FOR_MUMMICHOG; NEG_NOBG_FOR_UNKNOWN_SEARCH}

// // Background subtraction
// process blank_subtraction {

//     publishDir "${params.outdir}/peak_table/", mode: 'copy'
//     echo true

//     input:
//     file data_pos from POS_NOBG_FOR_BS
//     file pos_design from POS_DESIGN_FOR_BS
//     file data_neg from NEG_NOBG_FOR_BS
//     file neg_design from NEG_DESIGN_FOR_BS

//     output:
//     file params.pos_data_withbg into POS_DATA_WITHBG
//     file params.neg_data_withbg into NEG_DATA_WITHBG

//     when:
//     params.bs == "1"

//     script:
//     """   
//     blank_subtraction.py -i ${data_pos} -d ${pos_design} -o ${params.pos_data_withbg} &&
//     blank_subtraction.py -i ${data_neg} -d ${neg_design} -o ${params.neg_data_withbg} 
//     """
// }


// // split channel content for multiple-time use
// POS_DATA_WITHBG.into{POS_WITHBG_FOR_MQC; POS_WITHBG_FOR_PCA; POS_WITHBG_FOR_HCLUSTERING; POS_WITHBG_FOR_VD; POS_WITHBG_FOR_BARPLOT; POS_WITHBG_FOR_MUMMICHOG; POS_WITHBG_FOR_UNKNOWN_SEARCH}
// NEG_DATA_WITHBG.into{NEG_WITHBG_FOR_MQC; NEG_WITHBG_FOR_PCA; NEG_WITHBG_FOR_HCLUSTERING; NEG_WITHBG_FOR_VD; NEG_WITHBG_FOR_BARPLOT; NEG_WITHBG_FOR_MUMMICHOG; NEG_WITHBG_FOR_UNKNOWN_SEARCH}

// // Process for generating files that can be parsed by MultiQC regarding peak numbers of different steps.
// process mqc_peak_number_comparison {

//     publishDir "${params.outdir}/mqc/", mode: 'copy'
//     echo true

//     input:
//     file pos_nobg from POS_NOBG_FOR_MQC
//     file neg_nobg from NEG_NOBG_FOR_MQC
//     file pos_withbg from POS_WITHBG_FOR_MQC
//     file neg_withbg from NEG_WITHBG_FOR_MQC

//     output:
//     file params.peak_number_comparison_mqc into PEAK_NUMBER_COMPARISON_MQC

//     when:
//     params.bs == "1"

//     script:
//     """
//     peak_number_comparison.py -i1 ${pos_nobg} -i2 ${neg_nobg} -i3 ${pos_withbg} -i4 ${neg_withbg} -o ${params.peak_number_comparison_mqc}
//     """
// }

// // process for PCA of "no background subtraction" results
// process pca_nobg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_NOBG_FOR_PCA
//     file pos_design from POS_DESIGN_FOR_PCA_NOBG
//     file data_neg from NEG_NOBG_FOR_PCA
//     file neg_design from NEG_DESIGN_FOR_PCA_NOBG

//     output:
//     file params.pca_pos_nobg into PCA_POS_NOBG
//     file params.pca_neg_nobg into PCA_NEG_NOBG

//     script:
//     """   
//     pca.py -i ${data_pos} -d ${pos_design} -o ${params.pca_pos_nobg} &&
//     pca.py -i ${data_neg} -d ${neg_design} -o ${params.pca_neg_nobg} 
//     """
// }

// // process for PCA of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
// process pca_withbg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_WITHBG_FOR_PCA
//     file pos_design from POS_DESIGN_FOR_PCA_WITHBG
//     file data_neg from NEG_WITHBG_FOR_PCA
//     file neg_design from NEG_DESIGN_FOR_PCA_WITHBG

//     output:
//     file params.pca_pos_withbg into PCA_POS_WITHBG
//     file params.pca_neg_withbg into PCA_NEG_WITHBG

//     when:
//     params.bs == "1"

//     script:
//     """   
//     pca.py -i ${data_pos} -d ${pos_design} -o ${params.pca_pos_withbg} &&
//     pca.py -i ${data_neg} -d ${neg_design} -o ${params.pca_neg_withbg}
//     """
// }

// // process for hierarchical clustering of "no background subtraction" results
// process h_clustering_nobg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_NOBG_FOR_HCLUSTERING
//     file pos_design from POS_DESIGN_FOR_HCLUSTERING_NOBG
//     file data_neg from NEG_NOBG_FOR_HCLUSTERING
//     file neg_design from NEG_DESIGN_FOR_HCLUSTERING_NOBG

//     output:
//     file params.hclustering_pos_nobg into HCLUSTERING_POS_NOBG
//     file params.hclustering_neg_nobg into HCLUSTERING_NEG_NOBG

//     script:
//     """   
//     h_clustering.py -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_nobg} -m 0 &&
//     h_clustering.py -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_nobg} -m 0
//     """
// }

// // process for hierarchical clustering of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
// process h_clustering_withbg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_WITHBG_FOR_HCLUSTERING
//     file pos_design from POS_DESIGN_FOR_HCLUSTERING_WITHBG
//     file data_neg from NEG_WITHBG_FOR_HCLUSTERING
//     file neg_design from NEG_DESIGN_FOR_HCLUSTERING_WITHBG

//     when:
//     params.bs == "1"

//     output:
//     file params.hclustering_pos_withbg into HCLUSTERING_POS_WITHBG
//     file params.hclustering_neg_withbg into HCLUSTERING_NEG_WITHBG

//     script:
//     """   
//     h_clustering.py -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_withbg} -m 0 &&
//     h_clustering.py -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_withbg} -m 0
//     """
// }

// // process for venn diagram of "no background subtraction" results
// process venn_diagram_nobg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_NOBG_FOR_VD
//     file pos_design from POS_DESIGN_FOR_VD_NOBG
//     file data_neg from NEG_NOBG_FOR_VD
//     file neg_design from NEG_DESIGN_FOR_VD_NOBG

//     output:
//     file params.vd_pos_nobg into VD_POS_NOBG
//     file params.vd_neg_nobg into VD_NEG_NOBG
//     file params.pos_vd_group1_nobg into POS_VD_GROUP1_NOBG
//     file params.pos_vd_group2_nobg into POS_VD_GROUP2_NOBG
//     file params.pos_vd_both_nobg into POS_VD_BOTH_NOBG
//     file params.neg_vd_group1_nobg into NEG_VD_GROUP1_NOBG
//     file params.neg_vd_group2_nobg into NEG_VD_GROUP2_NOBG
//     file params.neg_vd_both_nobg into NEG_VD_BOTH_NOBG
//     file "pos*.txt" into POS_NOBG_CUTOFFS
//     file "neg*.txt" into NEG_NOBG_CUTOFFS

//     script:
//     """   
//     venn.py -i ${data_pos} -d ${pos_design} -o ${params.vd_pos_nobg} -bs 0 -g1 ${params.pos_vd_group1_nobg} -g2 ${params.pos_vd_group2_nobg} -bt ${params.pos_vd_both_nobg} &&
//     venn.py -i ${data_neg} -d ${neg_design} -o ${params.vd_neg_nobg} -bs 0 -g1 ${params.neg_vd_group1_nobg} -g2 ${params.neg_vd_group2_nobg} -bt ${params.neg_vd_both_nobg}
//     """
// }

// // process for venn diagram of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
// process venn_diagram_withbg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_WITHBG_FOR_VD
//     file pos_design from POS_DESIGN_FOR_VD_WITHBG
//     file data_neg from NEG_WITHBG_FOR_VD
//     file neg_design from NEG_DESIGN_FOR_VD_WITHBG

//     output:
//     file params.vd_pos_withbg into VD_POS_WITHBG
//     file params.vd_neg_withbg into VD_NEG_WITHBG
//     file params.pos_vd_group1_withbg into POS_VD_GROUP1_WITHBG
//     file params.pos_vd_group2_withbg into POS_VD_GROUP2_WITHBG
//     file params.pos_vd_both_withbg into POS_VD_BOTH_WITHBG
//     file params.neg_vd_group1_withbg into NEG_VD_GROUP1_WITHBG
//     file params.neg_vd_group2_withbg into NEG_VD_GROUP2_WITHBG
//     file params.neg_vd_both_withbg into NEG_VD_BOTH_WITHBG
//     file "pos*.txt" into POS_WITHBG_CUTOFFS
//     file "neg*.txt" into NEG_WITHBG_CUTOFFS

//     when:
//     params.bs == "1"

//     script:
//     """   
//     venn.py -i ${data_pos} -d ${pos_design} -o ${params.vd_pos_withbg} -bs 1 -g1 ${params.pos_vd_group1_withbg} -g2 ${params.pos_vd_group2_withbg} -bt ${params.pos_vd_both_withbg} &&
//     venn.py -i ${data_neg} -d ${neg_design} -o ${params.vd_neg_withbg} -bs 1 -g1 ${params.neg_vd_group1_withbg} -g2 ${params.neg_vd_group2_withbg} -bt ${params.neg_vd_both_withbg}
//     """
// }

// // process for bar plot of "no background subtraction" results
// process bar_plot_nobg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_NOBG_FOR_BARPLOT
//     file pos_design from POS_DESIGN_FOR_BARPLOT_NOBG
//     file data_neg from NEG_NOBG_FOR_BARPLOT
//     file neg_design from NEG_DESIGN_FOR_BARPLOT_NOBG

//     output:
//     file params.barplot_pos_nobg into BARPLOT_POS_NOBG
//     file params.barplot_neg_nobg into BARPLOT_NEG_NOBG

//     script:
//     """   
//     bar_plot.py -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_nobg} -m 0 -bs 0 &&
//     bar_plot.py -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_nobg} -m 0 -bs 0
//     """
// }

// // process for bar plot of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
// process bar_plot_withbg {
    
//     publishDir "${params.outdir}/figs", mode: 'copy'

//     input:
//     file data_pos from POS_WITHBG_FOR_BARPLOT
//     file pos_design from POS_DESIGN_FOR_BARPLOT_WITHBG
//     file data_neg from NEG_WITHBG_FOR_BARPLOT
//     file neg_design from NEG_DESIGN_FOR_BARPLOT_WITHBG

//     output:
//     file params.barplot_pos_withbg into BARPLOT_POS_WITHBG
//     file params.barplot_neg_withbg into BARPLOT_NEG_WITHBG

//     when:
//     params.bs == "1"

//     script:
//     """   
//     bar_plot.py -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_withbg} -m 0 -bs 1 &&
//     bar_plot.py -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_withbg} -m 0 -bs 1
//     """
// }

// // unknown search for metabolites identified before blank subtraction
// process unknown_search_nobg {
    
//     publishDir "${params.outdir}/peak_table/", mode: 'copy'

//     input:
//     file data_pos from POS_NOBG_FOR_UNKNOWN_SEARCH
//     file data_neg from NEG_NOBG_FOR_UNKNOWN_SEARCH

//     output:
//     file params.unknown_search_pos_nobg into UNKNOWN_SEARCH_POS_NOBG
//     file params.unknown_search_neg_nobg into UNKNOWN_SEARCH_NEG_NOBG

//     when:
//     params.unknown_search == "1"

//     script:
//     """   
//     unknown_search.R -i ${data_pos} -n positive -c ${params.mz_col_pos_nobg} -o ${params.unknown_search_pos_nobg} &&
//     unknown_search.R -i ${data_neg} -n negative -c ${params.mz_col_neg_nobg} -o ${params.unknown_search_neg_nobg}
//     """
// }

// // unknown search for metabolites identified after blank subtraction
// process unknown_search_withbg {
    
//     publishDir "${params.outdir}/peak_table/", mode: 'copy'

//     input:
//     file data_pos from POS_WITHBG_FOR_UNKNOWN_SEARCH
//     file data_neg from NEG_WITHBG_FOR_UNKNOWN_SEARCH

//     output:
//     file params.unknown_search_pos_withbg into UNKNOWN_SEARCH_POS_WITHBG
//     file params.unknown_search_neg_withbg into UNKNOWN_SEARCH_NEG_WITHBG

//     when:
//     params.bs == "1" && params.unknown_search == "1"

//     script:
//     """   
//     unknown_search.R -i ${data_pos} -n positive -c ${params.mz_col_pos_withbg} -o ${params.unknown_search_pos_withbg} &&
//     unknown_search.R -i ${data_neg} -n negative -c ${params.mz_col_neg_withbg} -o ${params.unknown_search_neg_withbg}
//     """
// }

// process mqc_figs {

//     publishDir "${params.outdir}/MultiQC/assets", mode: 'copy'

//     input:
//     file pca_pos_nobg from PCA_POS_NOBG
//     file pca_neg_nobg from PCA_NEG_NOBG
//     file pca_pos_withbg from PCA_POS_WITHBG
//     file pca_neg_withbg from PCA_NEG_WITHBG
//     file hclustering_pos_nobg from HCLUSTERING_POS_NOBG
//     file hclustering_neg_nobg from HCLUSTERING_NEG_NOBG
//     file hclustering_pos_withbg from HCLUSTERING_POS_WITHBG
//     file hclustering_neg_withbg from HCLUSTERING_NEG_WITHBG
//     file vd_pos_nobg from VD_POS_NOBG
//     file vd_neg_nobg from VD_NEG_NOBG
//     file vd_pos_withbg from VD_POS_WITHBG
//     file vd_neg_withbg from VD_NEG_WITHBG
//     file barplot_pos_nobg from BARPLOT_POS_NOBG
//     file barplot_neg_nobg from BARPLOT_NEG_NOBG
//     file barplot_pos_withbg from BARPLOT_POS_WITHBG
//     file barplot_neg_withbg from BARPLOT_NEG_WITHBG

//     output:
//     file "*positive_with_background_subtraction_mqc.png" into MQC_FIGS

//     script:
//     """
//     mv $pca_pos_nobg "PCA_for_positive_no_background_subtraction_mqc.png" &&
//     mv $pca_neg_nobg "PCA_for_negative_no_background_subtraction_mqc.png" &&
//     mv $pca_pos_withbg "PCA_for_positive_with_background_subtraction_mqc.png" &&
//     mv $pca_neg_withbg "PCA_for_negative_with_background_subtraction_mqc.png" &&
//     mv $hclustering_pos_nobg "Hirerchical_clustering_for_positive_no_background_subtraction_mqc.png" &&
//     mv $hclustering_neg_nobg "Hirerchical_clustering_for_negative_no_background_subtraction_mqc.png" &&
//     mv $hclustering_pos_withbg "Hirerchical_clustering_for_positive_with_background_subtraction_mqc.png" &&
//     mv $hclustering_neg_withbg "Hirerchical_clustering_for_negative_with_background_subtraction_mqc.png" &&
//     mv $vd_pos_nobg "Venn_diagram_for_positive_no_background_subtraction_mqc.png" &&
//     mv $vd_neg_nobg "Venn_diagram_for_negative_no_background_subtraction_mqc.png" &&
//     mv $vd_pos_withbg "Venn_diagram_for_positive_with_background_subtraction_mqc.png" &&
//     mv $vd_neg_withbg "Venn_diagram_for_negative_with_background_subtraction_mqc.png" &&
//     mv $barplot_pos_nobg "Bar_plot_clustering_for_positive_no_background_subtraction_mqc.png" &&
//     mv $barplot_neg_nobg "Bar_plot_clustering_for_negative_no_background_subtraction_mqc.png" &&
//     mv $barplot_pos_withbg "Bar_plot_clustering_for_positive_with_background_subtraction_mqc.png" &&
//     mv $barplot_neg_withbg "Bar_plot_clustering_for_negative_with_background_subtraction_mqc.png"
//     """
// }

// // Process for running MultiQC and generating the report.
// process report_generator {

//     publishDir "${params.outdir}/MultiQC/", mode: 'copy'

//     input:
// //    file mqc_dir from MQC_DIR
//     file pos_data_info_mqc from POS_DATA_INFO_MQC
//     file neg_data_info_mqc from NEG_DATA_INFO_MQC
//     file experiments_info from EXPERIMENTS_INFO
//     file mqc_config from MQC_CONFIG
//     file peak_number_comparison_mqc from PEAK_NUMBER_COMPARISON_MQC
//     file '*' from MQC_FIGS

//     output:
//     file "multiqc_report.html" into MULTIQC_REPORT

//     script:
//     """
//     multiqc .
//     """

// }

// if (params.container != "Docker") {
//     MAT_CONFIG_DIR = Channel.from('~/.config/matplotlib/')
//     MAT_CONFIG_FILE = Channel.from('~/.config/matplotlib/matplotlibrc')
// }
// else {
//     MAT_CONFIG_DIR = Channel.from('/root/.config/matplotlib/')
//     MAT_CONFIG_FILE = Channel.from('/root/.config/matplotlib/matplotlibrc')
// }

// MAT_CONFIG_DIR.into{MAT_CONFIG_DIR_NOBG; MAT_CONFIG_DIR_WITHBG}
// MAT_CONFIG_FILE.into{MAT_CONFIG_FILE_NOBG; MAT_CONFIG_FILE_WITHBG}

// process mummichog_report_nobg {

//     publishDir "${params.outdir}/mummichog/before_blank_subtraction", mode: 'copy'

//     input:
//     file pos_vd_group1_nobg from POS_VD_GROUP1_NOBG
//     file pos_vd_group2_nobg from POS_VD_GROUP2_NOBG
//     file pos_vd_both_nobg from POS_VD_BOTH_NOBG
//     val mat_config_dir_nobg from MAT_CONFIG_DIR_NOBG
//     val mat_config_file_nobg from MAT_CONFIG_FILE_NOBG
//     file "*" from POS_NOBG_CUTOFFS

//     output:
//     file "*" into MUMMICHOG_REPORT_NOBG

//     script:
//     """
//     echo "generating mommichog report for peaks before blank subtraction" &&
//     mkdir -p ${mat_config_dir_nobg} &&
//     echo "backend: Agg" > ${mat_config_file_nobg} &&
//     mummichog_input_prepare.py -i ${pos_vd_group1_nobg} -o ${params.data_pos_nobg_group1_mummichog} &&
//     mummichog -f ${params.data_pos_nobg_group1_mummichog} -o ${params.data_pos_nobg_group1_mummichog_out} -c ${params.cutoff} &&
//     mummichog_input_prepare.py -i ${pos_vd_group2_nobg} -o ${params.data_pos_nobg_group2_mummichog} &&
//     mummichog -f ${params.data_pos_nobg_group2_mummichog} -o ${params.data_pos_nobg_group2_mummichog_out} -c ${params.cutoff} &&
//     mummichog_input_prepare.py -i ${pos_vd_both_nobg} -o ${params.data_pos_nobg_both_mummichog} &&
//     mummichog -f ${params.data_pos_nobg_both_mummichog} -o ${params.data_pos_nobg_both_mummichog_out} -c ${params.cutoff}
//     """

// }

// process mummichog_report_withbg {

//     publishDir "${params.outdir}/mummichog/after_blank_subtractionode", mode:'copy'

//     input:
//     file pos_vd_group1_withbg from POS_VD_GROUP1_WITHBG
//     file pos_vd_group2_withbg from POS_VD_GROUP2_WITHBG
//     file pos_vd_both_withbg from POS_VD_BOTH_WITHBG
//     val mat_config_dir_withbg from MAT_CONFIG_DIR_WITHBG
//     val mat_config_file_withbg from MAT_CONFIG_FILE_WITHBG
//     file "*" from POS_WITHBG_CUTOFFS

//     output:
//     file "*" into MUMMICHOG_REPORT_WITHBG

//     when:
//     params.bs == "1"

//     script:
//     """
//     echo "generating mommichog report for peaks after blank subtraction" &&
//     mkdir -p ${mat_config_dir_withbg} &&
//     echo "backend: Agg" > ${mat_config_file_withbg} &&
//     mummichog_input_prepare.py -i ${pos_vd_group1_withbg} -o ${params.data_pos_withbg_group1_mummichog} &&
//     mummichog -f ${params.data_pos_withbg_group1_mummichog} -o ${params.data_pos_withbg_group1_mummichog_out} -c ${params.cutoff} &&
//     mummichog_input_prepare.py -i ${pos_vd_group2_withbg} -o ${params.data_pos_withbg_group2_mummichog} &&
//     mummichog -f ${params.data_pos_withbg_group2_mummichog} -o ${params.data_pos_withbg_group2_mummichog_out} -c ${params.cutoff} &&
//     mummichog_input_prepare.py -i ${pos_vd_both_withbg} -o ${params.data_pos_withbg_both_mummichog} &&
//     mummichog -f ${params.data_pos_withbg_both_mummichog} -o ${params.data_pos_withbg_both_mummichog_out} -c ${params.cutoff}
//     """
// }
