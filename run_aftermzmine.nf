#!/usr/bin/env nextflow

/**
    UMPIRE: A Reproducible Untargeted Metabolomics Data Processing Pipeline
    Copyright (C) LemasLab         
          
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this script.  If not, see <http://www.gnu.org/licenses/>.
    
    For any bugs or problems found, please contact us at
    - [email placeholder]; 
    - [github placeholder]
*/

// Those variable names which are all uppercase are channel names

version='0.0.0'
timestamp='20200202'

MZMINE = Channel.fromPath(params.mzmine_dir, type: 'dir') // The location of folder of MzMine
MZMINE.into{POS_MZMINE; NEG_MZMINE} // Duplicate the MZMINE chennel into two channels, one of which deals with positive sample while the other deals with negative sample.
// BATCHFILE_GENERATOR = Channel.fromPath(params.batchfile_generator) // This channel stores Python code (~/src/batchfile_generator.py) for generating MzMine batchfile, which enables us to run MzMine in batch mode. 

POS_DATA_DIR = Channel.fromPath(params.input_dir_pos, type: 'dir') // Location of folder storing positive data
NEG_DATA_DIR = Channel.fromPath(params.input_dir_neg, type: 'dir') // Location of folder storing negative data

PYTHON_PCA = Channel.fromPath(params.python_pca) // Chennel of Python code for principle component analysis
PYTHON_PCA.into{PYTHON_PCA_NOBG; PYTHON_PCA_WITHBG} // Duplicate the above chennel to two channels, one the them processes result without background substraction, the other one processes processes result with background subtraction.

PYTHON_HCLUSTERING = Channel.fromPath(params.python_hclustering) // Chennel of Python code for hierarchical clustering
PYTHON_HCLUSTERING.into{PYTHON_HCLUSTERING_NOBG; PYTHON_HCLUSTERING_WITHBG}

PYTHON_VD = Channel.fromPath(params.python_vd) // Chennel of Python code for venn diagram
PYTHON_VD.into{PYTHON_VD_NOBG; PYTHON_VD_WITHBG}

PYTHON_BARPLOT = Channel.fromPath(params.python_barplot) // Chennel of Python code for venn diagram
PYTHON_BARPLOT.into{PYTHON_BARPLOT_NOBG; PYTHON_BARPLOT_WITHBG}

PYTHON_ADDSTATS = Channel.fromPath(params.python_addstats)

PYTHON_DATA_INFO = Channel.fromPath(params.data_info) // Python code for generating MultiQC file regarding data information including file name and file size.
PYTHON_PEAK_NUMBER_COMPARISON = Channel.fromPath(params.peak_number_comparison_path) // Python code for generating MultiQC file ragarding peak numbers for different background subtraction threshold.
PYTHON_MUMMICHOG_INPUT_PREPARE = Channel.fromPath(params.python_mummichog_input_prepare)
PYTHON_MUMMICHOG_INPUT_PREPARE.into{PYTHON_MUMMICHOG_INPUT_PREPARE_NOBG; PYTHON_MUMMICHOG_INPUT_PREPARE_WITHBG}

// Following is Python code for background subtraction.
PYTHON_BS = Channel.fromPath(params.python_bs)

POS_MZMINE_RESULT = Channel.fromPath(params.pos_mzmine_peak_output)
// POS_NOBG.into{POS_NOBG_FOR_AS; POS_NOBG_FOR_BS; POS_NOBG_FOR_PCA; POS_NOBG_FOR_HCLUSTERING; POS_NOBG_FOR_VD; POS_NOBG_FOR_BARPLOT}
NEG_MZMINE_RESULT = Channel.fromPath(params.neg_mzmine_peak_output)
// NEG_NOBG.into{NEG_NOBG_FOR_AS; NEG_NOBG_FOR_BS; NEG_NOBG_FOR_PCA; NEG_NOBG_FOR_HCLUSTERING; NEG_NOBG_FOR_VD; NEG_NOBG_FOR_BARPLOT}

// Library
POS_LIBRARY = Channel.fromPath(params.pos_library)
NEG_LIBRARY = Channel.fromPath(params.neg_library)
POS_LIBRARY.into{POS_LIBRARY_MZMINE; POS_LIBRARY_STAT}
NEG_LIBRARY.into{NEG_LIBRARY_MZMINE; NEG_LIBRARY_STAT}

// Design files for positive data and negative data.
POS_DESIGN = Channel.fromPath(params.POS_design_path)
POS_DESIGN.into{POS_DESIGN_FOR_AS; POS_DESIGN_FOR_BS; POS_DESIGN_FOR_PCA_NOBG; POS_DESIGN_FOR_PCA_WITHBG; POS_DESIGN_FOR_HCLUSTERING_NOBG; POS_DESIGN_FOR_HCLUSTERING_WITHBG; POS_DESIGN_FOR_VD_NOBG; POS_DESIGN_FOR_VD_WITHBG; POS_DESIGN_FOR_BARPLOT_NOBG; POS_DESIGN_FOR_BARPLOT_WITHBG}
NEG_DESIGN = Channel.fromPath(params.NEG_design_path)
NEG_DESIGN.into{NEG_DESIGN_FOR_AS; NEG_DESIGN_FOR_BS; NEG_DESIGN_FOR_PCA_NOBG; NEG_DESIGN_FOR_PCA_WITHBG; NEG_DESIGN_FOR_HCLUSTERING_NOBG; NEG_DESIGN_FOR_HCLUSTERING_WITHBG; NEG_DESIGN_FOR_VD_NOBG; NEG_DESIGN_FOR_VD_WITHBG; NEG_DESIGN_FOR_BARPLOT_NOBG; NEG_DESIGN_FOR_BARPLOT_WITHBG}

// Pre-build MultiQC report information
// EXPERIMENTS_INFO = Channel.fromPath(params.experiments_info)
// MQC_CONFIG = Channel.fromPath(params.mqc_config)

// Result files used by MultiQC to generate report.
// MQC_DIR = Channel.fromPath(params.mqc_dir, type: 'dir')

/**
    Prints version when asked for
*/

if (params.version) {
    System.out.println("")
    System.out.println("UMPIRE: A Reproducible Untargeted Metabolomics Data Processing Pipeline - Version: $version ($timestamp)")
    exit 1
}

/**
    Prints help when asked for
*/

if (params.help) {
    System.out.println("")
    System.out.println("UMPIRE: A Reproducible Untargeted Metabolomics Data Processing Pipeline - Version: $version ($timestamp)")
    System.out.println("This pipeline is distributed in the hope that it will be useful")
    System.out.println("but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.")
    System.out.println("")
    System.out.println("Please report comments and bugs to alessia.visconti@kcl.ac.uk")
    System.out.println("or at https://github.com/GalaxyDream/metabolomics_data_processing/issues.")
    System.out.println("Check https://github.com/GalaxyDream/metabolomics_data_processing for updates, and refer to")
    System.out.println("[wiki placeholder]")
    System.out.println("")
    System.out.println("Usage:  ")
    System.out.println("   nextflow run_all.nf [options] -with-docker galaxydream/metabolomics_pipeline")
    System.out.println("")
    System.out.println("Arguments (it is mandatory to change `input_file` and `mzmine_dir` before running:")
    System.out.println("----------------------------- common parameters ----------------------------------")
    System.out.println("    --input_dir_pos                         folder location for positive data, default is 'data/POS'")
    System.out.println("    --input_dir_neg                         folder location for positive data, default is 'data/NEG'")
    System.out.println("    --POS_design_path                       location for positive design file, default is 'data/pos_design.csv'")
    System.out.println("    --NEG_design_path                       location for negative design file, default is 'data/neg_design.csv'")
    System.out.println("    --pos_mzmine_peak_output                location for positive peak table generated by MZmine-2.53, default is 'pos_data.csv'")
    System.out.println("    --neg_mzmine_peak_output                location for negative peak table generated by MZmine-2.53, default is 'neg_data.csv'")
    System.out.println("    --version                               whether to show version information or not, default is null")
    System.out.println("    --help                                  whether to show help information or not, default is null")
    System.out.println("Please refer to nextflow.config for more options.")
    System.out.println("")
    System.out.println("Container:")
    System.out.println("    Docker image to use with -with-docker|-with-singularity options is")
    System.out.println("    'docker://galaxydream/metabolomics_pipeline'")
    System.out.println("")
    System.out.println("UMPIRE supports .mzXML format files.")
    System.out.println("")
    exit 1
}

// Prepare data information for multiQC.
process mqc_data_info {

    publishDir './results/mqc/', mode: 'copy'

    input:
    file get_data_info from PYTHON_DATA_INFO
    file pos_data_dir from POS_DATA_DIR
    file neg_data_dir from NEG_DATA_DIR


    // POS_DATA and NEG_DATA are channels containing filtered POS and NEG data, which are ready to be input to R codes.
    output:
    file params.pos_data_info_mqc into POS_DATA_INFO_MQC
    file params.neg_data_info_mqc into NEG_DATA_INFO_MQC

    shell:
    """
    python3 ${get_data_info} -i ${pos_data_dir} -o $params.pos_data_info_mqc -n p &&
    python3 ${get_data_info} -i ${neg_data_dir} -o $params.neg_data_info_mqc -n n
    """
}

process add_stats {

    publishDir './results/peak_table/', mode: 'copy'
    echo true

    input:
    file python_addstats from PYTHON_ADDSTATS
    file data_pos from POS_MZMINE_RESULT
    file pos_design from POS_DESIGN_FOR_AS
    file data_neg from NEG_MZMINE_RESULT
    file neg_design from NEG_DESIGN_FOR_AS
    file pos_library from POS_LIBRARY_STAT
    file neg_library from NEG_LIBRARY_STAT

    output:
    file params.pos_data_nobg into POS_DATA_NOBG
    file params.neg_data_nobg into NEG_DATA_NOBG

    shell:
    """   
    python3 ${python_addstats} -i ${data_pos} -d ${pos_design} -o ${params.pos_data_nobg} -l ${pos_library} &&
    python3 ${python_addstats} -i ${data_neg} -d ${neg_design} -o ${params.neg_data_nobg} -l ${neg_library}

    """
}

// split channel content for multiple-time use
POS_DATA_NOBG.into{POS_NOBG_FOR_BS; POS_NOBG_FOR_MQC; POS_NOBG_FOR_PCA; POS_NOBG_FOR_HCLUSTERING; POS_NOBG_FOR_VD; POS_NOBG_FOR_BARPLOT; POS_NOBG_FOR_MUMMICHOG}
NEG_DATA_NOBG.into{NEG_NOBG_FOR_BS; NEG_NOBG_FOR_MQC; NEG_NOBG_FOR_PCA; NEG_NOBG_FOR_HCLUSTERING; NEG_NOBG_FOR_VD; NEG_NOBG_FOR_BARPLOT; NEG_NOBG_FOR_MUMMICHOG}

// Background subtraction
process blank_subtraction {

    publishDir './results/peak_table/', mode: 'copy'

    input:
    file python_bs from PYTHON_BS
    file data_pos from POS_NOBG_FOR_BS
    file pos_design from POS_DESIGN_FOR_BS
    file data_neg from NEG_NOBG_FOR_BS
    file neg_design from NEG_DESIGN_FOR_BS

    output:
    file params.pos_data_withbg into POS_DATA_WITHBG
    file params.neg_data_withbg into NEG_DATA_WITHBG

    shell:
    """   
    python3 ${python_bs} -i ${data_pos} -d ${pos_design} -o ${params.pos_data_withbg} &&
    python3 ${python_bs} -i ${data_neg} -d ${neg_design} -o ${params.neg_data_withbg} 

    """
}

// split channel content for multiple-time use
POS_DATA_WITHBG.into{POS_WITHBG_FOR_MQC; POS_WITHBG_FOR_PCA; POS_WITHBG_FOR_HCLUSTERING; POS_WITHBG_FOR_VD; POS_WITHBG_FOR_BARPLOT; POS_WITHBG_FOR_MUMMICHOG}
NEG_DATA_WITHBG.into{NEG_WITHBG_FOR_MQC; NEG_WITHBG_FOR_PCA; NEG_WITHBG_FOR_HCLUSTERING; NEG_WITHBG_FOR_VD; NEG_WITHBG_FOR_BARPLOT; NEG_WITHBG_FOR_MUMMICHOG}

// Process for generating files that can be parsed by MultiQC regarding peak numbers of different steps.
process mqc_peak_number_comparison {

    publishDir './results/mqc/', mode: 'copy'
    echo true

    input:
    file get_peak_number_comparison from PYTHON_PEAK_NUMBER_COMPARISON
    file pos_nobg from POS_NOBG_FOR_MQC
    file neg_nobg from NEG_NOBG_FOR_MQC
    file pos_withbg from POS_WITHBG_FOR_MQC
    file neg_withbg from NEG_WITHBG_FOR_MQC

    output:
    file params.peak_number_comparison_mqc into PEAK_NUMBER_COMPARISON_MQC

    when:
    params.bs == "1"

    shell:
    """
    python3 ${get_peak_number_comparison} -i1 ${pos_nobg} -i2 ${neg_nobg} -i3 ${pos_withbg} -i4 ${neg_withbg} -o ${params.peak_number_comparison_mqc}

    """
}

// process for PCA of "no background subtraction" results
process pca_nobg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_NOBG_FOR_PCA
    file pos_design from POS_DESIGN_FOR_PCA_NOBG
    file data_neg from NEG_NOBG_FOR_PCA
    file neg_design from NEG_DESIGN_FOR_PCA_NOBG
    file python_pca from PYTHON_PCA_NOBG

    output:
    file params.pca_pos_nobg into PCA_POS_NOBG
    file params.pca_neg_nobg into PCA_NEG_NOBG

    shell:
    """   
    python3 ${python_pca} -i ${data_pos} -d ${pos_design} -o ${params.pca_pos_nobg} &&
    python3 ${python_pca} -i ${data_neg} -d ${neg_design} -o ${params.pca_neg_nobg} 

    """

}

// process for PCA of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
process pca_withbg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_WITHBG_FOR_PCA
    file pos_design from POS_DESIGN_FOR_PCA_WITHBG
    file data_neg from NEG_WITHBG_FOR_PCA
    file neg_design from NEG_DESIGN_FOR_PCA_WITHBG
    file python_pca from PYTHON_PCA_WITHBG

    output:
    file params.pca_pos_withbg into PCA_POS_WITHBG
    file params.pca_neg_withbg into PCA_NEG_WITHBG

    shell:
    """   
    python3 ${python_pca} -i ${data_pos} -d ${pos_design} -o ${params.pca_pos_withbg} &&
    python3 ${python_pca} -i ${data_neg} -d ${neg_design} -o ${params.pca_neg_withbg} 

    """

}

// process for hierarchical clustering of "no background subtraction" results
process h_clustering_nobg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_NOBG_FOR_HCLUSTERING
    file pos_design from POS_DESIGN_FOR_HCLUSTERING_NOBG
    file data_neg from NEG_NOBG_FOR_HCLUSTERING
    file neg_design from NEG_DESIGN_FOR_HCLUSTERING_NOBG
    file python_hclustering from PYTHON_HCLUSTERING_NOBG

    output:
    file params.hclustering_pos_nobg into HCLUSTERING_POS_NOBG
    file params.hclustering_neg_nobg into HCLUSTERING_NEG_NOBG

    shell:
    """   
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_nobg} -m 0 -bs 0 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_nobg} -m 0 -bs 0

    """

}

// process for hierarchical clustering of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
process h_clustering_withbg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_WITHBG_FOR_HCLUSTERING
    file pos_design from POS_DESIGN_FOR_HCLUSTERING_WITHBG
    file data_neg from NEG_WITHBG_FOR_HCLUSTERING
    file neg_design from NEG_DESIGN_FOR_HCLUSTERING_WITHBG
    file python_hclustering from PYTHON_HCLUSTERING_WITHBG

    output:
    file params.hclustering_pos_withbg into HCLUSTERING_POS_WITHBG
    file params.hclustering_neg_withbg into HCLUSTERING_NEG_WITHBG

    shell:
    """   
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_withbg} -m 0 -bs 1 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_withbg} -m 0 -bs 1

    """

}

// process for venn diagram of "no background subtraction" results
process venn_diagram_nobg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_NOBG_FOR_VD
    file pos_design from POS_DESIGN_FOR_VD_NOBG
    file data_neg from NEG_NOBG_FOR_VD
    file neg_design from NEG_DESIGN_FOR_VD_NOBG
    file python_vd from PYTHON_VD_NOBG

    output:
    file params.vd_pos_nobg into VD_POS_NOBG
    file params.vd_neg_nobg into VD_NEG_NOBG
    file params.pos_vd_group1_nobg into POS_VD_GROUP1_NOBG
    file params.pos_vd_group2_nobg into POS_VD_GROUP2_NOBG
    file params.pos_vd_both_nobg into POS_VD_BOTH_NOBG
    file params.neg_vd_group1_nobg into NEG_VD_GROUP1_NOBG
    file params.neg_vd_group2_nobg into NEG_VD_GROUP2_NOBG
    file params.neg_vd_both_nobg into NEG_VD_BOTH_NOBG
    file "pos*.txt" into POS_NOBG_CUTOFFS
    file "neg*.txt" into NEG_NOBG_CUTOFFS

    shell:
    """   
    python3 ${python_vd} -i ${data_pos} -d ${pos_design} -o ${params.vd_pos_nobg} -bs 0 -g1 ${params.pos_vd_group1_nobg} -g2 ${params.pos_vd_group2_nobg} -bt ${params.pos_vd_both_nobg} &&
    python3 ${python_vd} -i ${data_neg} -d ${neg_design} -o ${params.vd_neg_nobg} -bs 0 -g1 ${params.neg_vd_group1_nobg} -g2 ${params.neg_vd_group2_nobg} -bt ${params.neg_vd_both_nobg}

    """

}

// process for venn diagram of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
process venn_diagram_withbg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_WITHBG_FOR_VD
    file pos_design from POS_DESIGN_FOR_VD_WITHBG
    file data_neg from NEG_WITHBG_FOR_VD
    file neg_design from NEG_DESIGN_FOR_VD_WITHBG
    file python_vd from PYTHON_VD_WITHBG

    output:
    file params.vd_pos_withbg into VD_POS_WITHBG
    file params.vd_neg_withbg into VD_NEG_WITHBG
    file params.pos_vd_group1_withbg into POS_VD_GROUP1_WITHBG
    file params.pos_vd_group2_withbg into POS_VD_GROUP2_WITHBG
    file params.pos_vd_both_withbg into POS_VD_BOTH_WITHBG
    file params.neg_vd_group1_withbg into NEG_VD_GROUP1_WITHBG
    file params.neg_vd_group2_withbg into NEG_VD_GROUP2_WITHBG
    file params.neg_vd_both_withbg into NEG_VD_BOTH_WITHBG
    file "pos*.txt" into POS_WITHBG_CUTOFFS
    file "neg*.txt" into NEG_WITHBG_CUTOFFS

    shell:
    """   
    python3 ${python_vd} -i ${data_pos} -d ${pos_design} -o ${params.vd_pos_withbg} -bs 1 -g1 ${params.pos_vd_group1_withbg} -g2 ${params.pos_vd_group2_withbg} -bt ${params.pos_vd_both_withbg} &&
    python3 ${python_vd} -i ${data_neg} -d ${neg_design} -o ${params.vd_neg_withbg} -bs 1 -g1 ${params.neg_vd_group1_withbg} -g2 ${params.neg_vd_group2_withbg} -bt ${params.neg_vd_both_withbg}

    """

}

// process for bar plot of "no background subtraction" results
process bar_plot_nobg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_NOBG_FOR_BARPLOT
    file pos_design from POS_DESIGN_FOR_BARPLOT_NOBG
    file data_neg from NEG_NOBG_FOR_BARPLOT
    file neg_design from NEG_DESIGN_FOR_BARPLOT_NOBG
    file python_barplot from PYTHON_BARPLOT_NOBG

    output:
    file params.barplot_pos_nobg into BARPLOT_POS_NOBG
    file params.barplot_neg_nobg into BARPLOT_NEG_NOBG

    shell:
    """   
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_nobg} -m 0 -bs 0 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_nobg} -m 0 -bs 0

    """

}

// process for bar plot of "with background subtraction" results, here we use 100 as the threshold of background subtraction.
process bar_plot_withbg {
    
    publishDir './results/figs', mode: 'copy'

    input:
    file data_pos from POS_WITHBG_FOR_BARPLOT
    file pos_design from POS_DESIGN_FOR_BARPLOT_WITHBG
    file data_neg from NEG_WITHBG_FOR_BARPLOT
    file neg_design from NEG_DESIGN_FOR_BARPLOT_WITHBG
    file python_barplot from PYTHON_BARPLOT_WITHBG

    output:
    file params.barplot_pos_withbg into BARPLOT_POS_WITHBG
    file params.barplot_neg_withbg into BARPLOT_NEG_WITHBG

    shell:
    """   
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_withbg} -m 0 -bs 1 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_withbg} -m 0 -bs 1

    """

}

if (params.use_singularity == "1") {
    mat_config_dir = "~/.config/matplotlib/"
}
else {
    mat_config_dir = "/root/.config/matplotlib/"
}

mat_confit_file = mat_config_dir + "matplotlibrc"

process mummichog_report_nobg {

    publishDir './results/mummichog/before_blank_subtraction', mode: 'copy'

    input:

    file python_mummichog_input_prepare from PYTHON_MUMMICHOG_INPUT_PREPARE_NOBG
    file pos_vd_group1_nobg from POS_VD_GROUP1_NOBG
    file pos_vd_group2_nobg from POS_VD_GROUP2_NOBG
    file pos_vd_both_nobg from POS_VD_BOTH_NOBG
    file "*" from POS_NOBG_CUTOFFS

    output:
    file "*" into MUMMICHOG_REPORT_NOBG

    shell:
    """
    echo "generating mommichog report for peaks before blank subtraction" &&
    mkdir -p !{mat_config_dir} &&
    echo "backend: Agg" > !{mat_config_file} &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_group1_nobg} -o !{params.data_pos_nobg_group1_mummichog} &&
    mummichog -f !{params.data_pos_nobg_group1_mummichog} -o !{params.data_pos_nobg_group1_mummichog_out} -c 0.05 &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_group2_nobg} -o !{params.data_pos_nobg_group2_mummichog} &&
    mummichog -f !{params.data_pos_nobg_group2_mummichog} -o !{params.data_pos_nobg_group2_mummichog_out} -c 0.05 &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_both_nobg} -o !{params.data_pos_nobg_both_mummichog} &&
    mummichog -f !{params.data_pos_nobg_both_mummichog} -o !{params.data_pos_nobg_both_mummichog_out} -c 0.05
    """

}

process mummichog_report_withbg {

    publishDir './results/mummichog/after_blank_subtraction', mode: 'copy'

    input:

    file python_mummichog_input_prepare from PYTHON_MUMMICHOG_INPUT_PREPARE_WITHBG
    file pos_vd_group1_withbg from POS_VD_GROUP1_WITHBG
    file pos_vd_group2_withbg from POS_VD_GROUP2_WITHBG
    file pos_vd_both_withbg from POS_VD_BOTH_WITHBG
    file "*" from POS_WITHBG_CUTOFFS

    output:
    file "*" into MUMMICHOG_REPORT_WITHBG

    when:
    params.bs == "1"

    shell:
    """
    echo "generating mommichog report for peaks after blank subtraction" &&
    mkdir -p !{mat_config_dir} &&
    echo "backend: Agg" > !{mat_config_file} &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_group1_withbg} -o !{params.data_pos_withbg_group1_mummichog} &&
    mummichog -f !{params.data_pos_withbg_group1_mummichog} -o !{params.data_pos_withbg_group1_mummichog_out} -c 0.05 &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_group2_withbg} -o !{params.data_pos_withbg_group2_mummichog} &&
    mummichog -f !{params.data_pos_withbg_group2_mummichog} -o !{params.data_pos_withbg_group2_mummichog_out} -c 0.05 &&
    python3 !{python_mummichog_input_prepare} -i !{pos_vd_both_withbg} -o !{params.data_pos_withbg_both_mummichog} &&
    mummichog -f !{params.data_pos_withbg_both_mummichog} -o !{params.data_pos_withbg_both_mummichog_out} -c 0.05
    """

}