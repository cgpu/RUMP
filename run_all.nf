#!/usr/bin/env nextflow

/**
    Metabolomics Pipeline for Suptercomputers (MPS)
    Copyright (C) lemas-research-group          
          
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
timestamp='20190916'

MZMINE = Channel.fromPath(params.mzmine_dir, type: 'dir') // The location of folder of MzMine
MZMINE.into{POS_MZMINE; NEG_MZMINE} // Duplicate the MZMINE chennel into two channels, one of which deals with positive sample while the other deals with negative sample.
BATCHFILE_GENERATOR_POS = Channel.fromPath(params.batchfile_generator_pos) // This channel stores Python code (~/src/batchfile_generator_pos.py) for generating MzMine batchfile for positive samples, which enables us to run MzMine in batch mode. 
BATCHFILE_GENERATOR_NEG = Channel.fromPath(params.batchfile_generator_neg) // This channel stores Python code (~/src/batchfile_generator_neg.py) for generating MzMine batchfile for negative samples, which enables us to run MzMine in batch mode. 

// MZMINE_RESULT_EXTRACTION = Channel.fromPath(params.mzmine_result_extraction) // This channel stores Python code for formatting the result generated by MzMine enabling the intermidiate result to be feeded appropriately into background subtraction code.

POS_DATA_DIR = Channel.fromPath(params.input_dir_pos, type: 'dir') // Location of folder storing positive data
POS_DATA_DIR.into{POS_DATA_DIR_INFO; POS_DATA_DIR_BS}
NEG_DATA_DIR = Channel.fromPath(params.input_dir_neg, type: 'dir') // Location of folder storing negative data
NEG_DATA_DIR.into{NEG_DATA_DIR_INFO; NEG_DATA_DIR_BS}

PYTHON_VD = Channel.fromPath(params.python_vd) // Chennel of Python code for venn diagram
PYTHON_VD.into{PYTHON_VD_NOBG; PYTHON_VD_WITHBG}

PYTHON_BARPLOT = Channel.fromPath(params.python_barplot) // Chennel of Python code for venn diagram
PYTHON_BARPLOT.into{PYTHON_BARPLOT_NOBG; PYTHON_BARPLOT_WITHBG}

PYTHON_ADDSTATS = Channel.fromPath(params.python_addstats)

PYTHON_PCA = Channel.fromPath(params.python_pca) // Chennel of Python code for principle component analysis
PYTHON_PCA.into{PYTHON_PCA_NOBG; PYTHON_PCA_WITHBG} // Duplicate the above chennel to two channels, one the them processes result without background substraction, the other one processes processes result with background subtraction.

PYTHON_HCLUSTERING = Channel.fromPath(params.python_hclustering) // Chennel of Python code for hierarchical clustering
PYTHON_HCLUSTERING.into{PYTHON_HCLUSTERING_NOBG; PYTHON_HCLUSTERING_WITHBG}

PYTHON_DATA_INFO = Channel.fromPath(params.data_info) // Python code for generating MultiQC file regarding data information including file name and file size.
PYTHON_PEAK_NUMBER_COMPARISON = Channel.fromPath(params.peak_number_comparison_path) // Python code for generating MultiQC file ragarding peak numbers for different background subtraction threshold.
PYTHON_MUMMICHOG_INPUT_PREPARE = Channel.fromPath(params.python_mummichog_input_prepare)

// Following is Python code for background subtraction.
PYTHON_BS = Channel.fromPath(params.python_bs)

// Design files for positive data and negative data.
POS_DESIGN = Channel.fromPath(params.POS_design_path)
POS_DESIGN.into{POS_DESIGN_FOR_AS; POS_DESIGN_FOR_BS; POS_DESIGN_FOR_PCA_NOBG; POS_DESIGN_FOR_PCA_WITHBG; POS_DESIGN_FOR_HCLUSTERING_NOBG; POS_DESIGN_FOR_HCLUSTERING_WITHBG; POS_DESIGN_FOR_VD_NOBG; POS_DESIGN_FOR_VD_WITHBG; POS_DESIGN_FOR_BARPLOT_NOBG; POS_DESIGN_FOR_BARPLOT_WITHBG}
NEG_DESIGN = Channel.fromPath(params.NEG_design_path)
NEG_DESIGN.into{NEG_DESIGN_FOR_AS; NEG_DESIGN_FOR_BS; NEG_DESIGN_FOR_PCA_NOBG; NEG_DESIGN_FOR_PCA_WITHBG; NEG_DESIGN_FOR_HCLUSTERING_NOBG; NEG_DESIGN_FOR_HCLUSTERING_WITHBG; NEG_DESIGN_FOR_VD_NOBG; NEG_DESIGN_FOR_VD_WITHBG; NEG_DESIGN_FOR_BARPLOT_NOBG; NEG_DESIGN_FOR_BARPLOT_WITHBG}

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
    System.out.println("METAGENOMIC PIPELINE FOR SUPERCOMPUTERS (MPS) - Version: $version ($timestamp)")
    exit 1
}

/**
    Prints help when asked for
*/

if (params.help) {
    System.out.println("")
    System.out.println("Metabolomics Pipeline for Suptercomputers (MPS) - Version: $version ($timestamp)")
    System.out.println("This pipeline is distributed in the hope that it will be useful")
    System.out.println("but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.")
    System.out.println("")
    System.out.println("Please report comments and bugs to alessia.visconti@kcl.ac.uk")
    System.out.println("or at https://github.com/GalaxyDream/metabolomics_data_processing/issues.")
    System.out.println("Check https://github.com/GalaxyDream/metabolomics_data_processing for updates, and refer to")
    System.out.println("[wiki placeholder]")
    System.out.println("")
    System.out.println("Usage:  ")
    System.out.println("   nextflow run.nf [options] -with-docker galaxydream/mzmine_oldversion")
    System.out.println("")
    System.out.println("Arguments (it is mandatory to change `input_file` and `mzmine_dir` before running:")
    System.out.println("----------------------------- common parameters ----------------------------------")
    System.out.println("    --version                               whether to show version information or not, default is null")
    System.out.println("    --help                                  whether to show help information or not, default is null")
    System.out.println("    --input_dir_pos                         location of your input metabolomics data (positive part) folder")
    System.out.println("    --input_dir_neg                         location of your input metabolomics data (negative part) folder")
    System.out.println("    --mzmine_dir                            location of your MzMine software, here we use MzMine-2.28")
    System.out.println("    --library                               location of your customized library for matching matabolites")
    System.out.println("    --pos_config                            location of generated MzMine batchfile for positive data")
    System.out.println("    --neg_config                            location of generated MzMine batchfile for negative data")
    System.out.println("    --batchfile_generator                   location of the Python code for generating MzMine batchfile")
    System.out.println("    --mzmine_result_extraction              location of the Python code for reformating MzMine processed result so that those result files can be input appropriately to R scripts for background subtraction")
    System.out.println("    --raw_stats_merge_path                  location of the Python code for merging files of raw intension and the corresponding stats")
    System.out.println("    --venn_path                             location of the Python code for generating venn diagram")
    System.out.println("    --python_pca_path                       location of the Python code for generating PCA analysis result")
    System.out.println("    --data_info                             location of the Python code for generating file regarding data information that can be parsed by MultiQC")
    System.out.println("    --peak_number_comparison_path           location of the Python code for generating file that can be parsed by MultiQC regarding identified number of peaks corresponding to different background subtraction threshold")
    System.out.println("    --R_01_path                             Part of the R script for background subtraction")
    System.out.println("    --R_02_threshold_path                   Part of the R script for background subtraction")
    System.out.println("    --R_02_blankbased_path                  Part of the R script for background subtraction")
    System.out.println("    --R_03_anova_path                       Part of the R script for background subtraction")
    System.out.println("    --R_03_kruskalwallis_path               Part of the R script for background subtraction")
    System.out.println("    --R_04_path                             Part of the R script for background subtraction")
    System.out.println("    --mqc_dir                               path of the folder containing all files that need to be parsed by MultiQC to generate the report")
    System.out.println("    --pos_data_info_mqc                     name of the file ending with '.yaml' regarding positive data information that can be parsed by MultiQC")
    System.out.println("    --neg_data_info_mqc                     name of the file ending with '.yaml' regarding negative data information that can be parsed by MultiQC")
    System.out.println("    --POS_design_path                       location of the input file design file for positive data background subtraction")
    System.out.println("    --NEG_design_path                       location of the input file design file for negative data background subtraction")
    System.out.println("    --pos_nobg_ready                        name of the positive peak detection result (csv file) without background subtraction")
    System.out.println("    --neg_nobg_ready                        name of the negative peak detection result (csv file) without background subtraction")
    System.out.println("    --pos_005_withbg_ready                  name of the positive peak detection result (csv file) with a background subtraction threshold of 005")
    System.out.println("    --neg_005_withbg_ready                  name of the negative peak detection result (csv file) with a background subtraction threshold of 005")
    System.out.println("    --pos_100_withbg_ready                  name of the positive peak detection result (csv file) with a background subtraction threshold of 100")
    System.out.println("    --neg_100_withbg_ready                  name of the negative peak detection result (csv file) with a background subtraction threshold of 100")
    System.out.println("    --pos_200_withbg_ready                  name of the positive peak detection result (csv file) with a background subtraction threshold of 200")
    System.out.println("    --neg_200_withbg_ready                  name of the negative peak detection result (csv file) with a background subtraction threshold of 200")
    System.out.println("Please refer to nextflow.config for more options.")
    System.out.println("")
    System.out.println("Container:")
    System.out.println("    Docker image to use with -with-docker|-with-singularity options is")
    System.out.println("    'docker://galaxydream/bioconductor_metabolomics'")
    System.out.println("")
    System.out.println("MPS supports mzXML format files.")
    System.out.println("")
    exit 1
}

// Process for generating MultiQC report regarding data information
process mqc_data_info {

    publishDir './results/mqc/', mode: 'copy'

    input:
    file get_data_info from PYTHON_DATA_INFO // Python code for generating MultiQC file regarding data information including file name and file size.
    file pos_data_dir from POS_DATA_DIR_INFO // Location of positive data
    file neg_data_dir from NEG_DATA_DIR_INFO // Location of negative data

    // POS_DATA and NEG_DATA are channels containing filtered POS and NEG data, which are ready to be input to R codes.
    output:
    file params.pos_data_info_mqc into POS_DATA_INFO_MQC // file regarding positive data information that can be parsed by MultiQC
    file params.neg_data_info_mqc into NEG_DATA_INFO_MQC // file regarding negative data information that can be parsed by MultiQC

    shell:
    """
    python3 ${get_data_info} -i ${pos_data_dir} -o $params.pos_data_info_mqc -n p &&
    python3 ${get_data_info} -i ${neg_data_dir} -o $params.neg_data_info_mqc -n n
    """
}

process batchfile_generation_mzmine {

    echo true

    input:
    file batchfile_generator_pos from BATCHFILE_GENERATOR_POS 
    file batchfile_generator_neg from BATCHFILE_GENERATOR_NEG
    file pos_data_dir from POS_DATA_DIR_BS // Location of positive data
    file neg_data_dir from NEG_DATA_DIR_BS // Location of negative data

    output:
    file params.pos_config into POS_BATCHFILE // Generated batchfile for processing positive data
    file params.neg_config into NEG_BATCHFILE // Generated batchfile for processing negative data

    shell:
    """  
    echo "setting parameters for MZmine" &&
    python ${batchfile_generator_pos} -x ${params.pos_config} -i ${pos_data_dir} -l $params.library -o $params.pos_mzmine_peak_output &&
    python ${batchfile_generator_neg} -x ${params.neg_config} -i ${neg_data_dir} -l $params.library -o $params.neg_mzmine_peak_output

    """
}

process pos_peakDetection_mzmine {

    echo true

    input:
    file p_b from POS_BATCHFILE // Batchfile for MzMine to process positive data.
    file p_m from POS_MZMINE // Folder of MzMine tool

    output:
    file "MZmine-2.53-Linux/${params.pos_mzmine_peak_output}" into POS_MZMINE_RESULT // MzMine processing result for positive data.

// Change "startMZmine_Linux.sh" to "startMZmine_MacOSX.command" in the following code if running locally with Mac

    shell:
    """   
    echo "peak detection and library matching for positive data" &&
    mv ${p_b} ${p_m} && cd ${p_m} && ./startMZmine-Linux ${p_b}
    """
}

process neg_peakDetection_mzmine {

    input:
    file n_b from NEG_BATCHFILE // Batchfile for MzMine to process negative data.
    file n_m from NEG_MZMINE // Folder of MzMine tool

    output:
    file "MZmine-2.53-Linux/${params.neg_mzmine_peak_output}" into NEG_MZMINE_RESULT // MzMine processing result for negative data.
    stdout result

// Change "startMZmine_Linux.sh" to "startMZmine_MacOSX.command" in the following code if running locally with Mac

    shell:
    """   
    echo "peak detection and library matching for negative data" &&
    mv ${n_b} ${n_m} && cd ${n_m} && ./startMZmine-Linux ${n_b}
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

    output:
    file params.pos_data_nobg into POS_DATA_NOBG
    file params.neg_data_nobg into NEG_DATA_NOBG

    shell:
    """   
    python3 ${python_addstats} -i ${data_pos} -d ${pos_design} -o ${params.pos_data_nobg} -l $params.library &&
    python3 ${python_addstats} -i ${data_neg} -d ${neg_design} -o ${params.neg_data_nobg} -l $params.library

    """
}

POS_DATA_NOBG.into{POS_NOBG_FOR_BS; POS_NOBG_FOR_MQC; POS_NOBG_FOR_PCA; POS_NOBG_FOR_HCLUSTERING; POS_NOBG_FOR_VD; POS_NOBG_FOR_BARPLOT; POS_NOBG_FOR_MUMMICHOG}
NEG_DATA_NOBG.into{NEG_NOBG_FOR_BS; NEG_NOBG_FOR_MQC; NEG_NOBG_FOR_PCA; NEG_NOBG_FOR_HCLUSTERING; NEG_NOBG_FOR_VD; NEG_NOBG_FOR_BARPLOT; NEG_NOBG_FOR_MUMMICHOG}

// Background subtraction
process blank_subtraction {

    publishDir './results/peak_table/', mode: 'copy'
    echo true

    input:
    file python_bs from PYTHON_BS
    file data_pos from POS_NOBG_FOR_BS
    file pos_design from POS_DESIGN_FOR_BS
    file data_neg from NEG_NOBG_FOR_BS
    file neg_design from NEG_DESIGN_FOR_BS

    when:
    params.bs == "1"

    output:
    file params.pos_data_withbg into POS_DATA_WITHBG
    file params.neg_data_withbg into NEG_DATA_WITHBG

    when:
    params.bs == "1"

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

    when:
    params.bs == "1"

    output:
    file params.pca_pos_withbg into PCA_POS_WITHBG
    file params.pca_neg_withbg into PCA_NEG_WITHBG

    when:
    params.bs == "1"

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
    file params.hclustering_pos_nobg_om into HCLUSTERING_POS_NOBG_OM
    file params.hclustering_neg_nobg into HCLUSTERING_NEG_NOBG
    file params.hclustering_neg_nobg_om into HCLUSTERING_NEG_NOBG_OM

    shell:
    """   
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_nobg} -m 0 -bs 0 &&
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_nobg_om} -m 1 -bs 0 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_nobg} -m 0 -bs 0 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_nobg_om} -m 1 -bs 0

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

    when:
    params.bs == "1"

    output:
    file params.hclustering_pos_withbg into HCLUSTERING_POS_WITHBG
    file params.hclustering_pos_withbg_om into HCLUSTERING_POS_WITHBG_OM
    file params.hclustering_neg_withbg into HCLUSTERING_NEG_WITHBG
    file params.hclustering_neg_withbg_om into HCLUSTERING_NEG_WITHBG_OM

    shell:
    """   
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_withbg} -m 0 -bs 1 &&
    python3 ${python_hclustering} -i ${data_pos} -d ${pos_design} -o ${params.hclustering_pos_withbg_om} -m 1 -bs 1 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_withbg} -m 0 -bs 1 &&
    python3 ${python_hclustering} -i ${data_neg} -d ${neg_design} -o ${params.hclustering_neg_withbg_om} -m 1 -bs 1

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
    file params.barplot_pos_nobg_om into BARPLOT_POS_NOBG_OM
    file params.barplot_neg_nobg into BARPLOT_NEG_NOBG
    file params.barplot_neg_nobg_om into BARPLOT_NEG_NOBG_OM

    shell:
    """   
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_nobg} -m 0 -bs 0 &&
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_nobg_om} -m 1 -bs 0 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_nobg} -m 0 -bs 0 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_nobg_om} -m 1 -bs 0 &&

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
    file params.barplot_pos_withbg_om into BARPLOT_POS_WITHBG_OM
    file params.barplot_neg_withbg into BARPLOT_NEG_WITHBG
    file params.barplot_neg_withbg_om into BARPLOT_NEG_WITHBG_OM

    shell:
    """   
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_withbg} -m 0 -bs 1 &&
    python3 ${python_barplot} -i ${data_pos} -d ${pos_design} -o ${params.barplot_pos_withbg_om} -m 1 -bs 1 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_withbg} -m 0 -bs 1 &&
    python3 ${python_barplot} -i ${data_neg} -d ${neg_design} -o ${params.barplot_neg_withbg_om} -m 1 -bs 1

    """

}

process mqc_figs {

    publishDir './results/mqc/', mode: 'copy'

    input:
    file pca_pos_nobg from PCA_POS_NOBG
    file pca_neg_nobg from PCA_NEG_NOBG
    file pca_pos_withbg from PCA_POS_WITHBG
    file pca_neg_withbg from PCA_NEG_WITHBG
    file hclustering_pos_nobg from HCLUSTERING_POS_NOBG
    file hclustering_pos_nobg_om from HCLUSTERING_POS_NOBG_OM
    file hclustering_neg_nobg from HCLUSTERING_NEG_NOBG
    file hclustering_neg_nobg_om from HCLUSTERING_NEG_NOBG_OM
    file hclustering_pos_withbg from HCLUSTERING_POS_WITHBG
    file hclustering_pos_withbg_om from HCLUSTERING_POS_WITHBG_OM
    file hclustering_neg_withbg from HCLUSTERING_NEG_WITHBG
    file hclustering_neg_withbg_om from HCLUSTERING_NEG_WITHBG_OM
    file vd_pos_nobg from VD_POS_NOBG
    file vd_neg_nobg from VD_NEG_NOBG
    file vd_pos_withbg from VD_POS_WITHBG
    file vd_neg_withbg from VD_NEG_WITHBG
    file barplot_pos_nobg from BARPLOT_POS_NOBG
    file barplot_pos_nobg_om from BARPLOT_POS_NOBG_OM
    file barplot_neg_nobg from BARPLOT_NEG_NOBG
    file barplot_neg_nobg_om from BARPLOT_NEG_NOBG_OM
    file barplot_pos_withbg from BARPLOT_POS_WITHBG
    file barplot_pos_withbg_om from BARPLOT_POS_WITHBG_OM
    file barplot_neg_withbg from BARPLOT_NEG_WITHBG
    file barplot_neg_withbg_om from BARPLOT_NEG_WITHBG_OM

    output:
    file "*positive_with_background_subtraction_mqc.png" into MQC_FIGS

    when:
    params.bs == "1"

    shell:
    """
    mv $pca_pos_nobg "PCA_for_positive_no_background_subtraction_mqc.png" &&
    mv $pca_neg_nobg "PCA_for_negative_no_background_subtraction_mqc.png" &&
    mv $pca_pos_withbg "PCA_for_positive_with_background_subtraction_mqc.png" &&
    mv $pca_neg_withbg "PCA_for_negative_with_background_subtraction_mqc.png" &&
    mv $hclustering_pos_nobg "Hirerchical_clustering_for_positive_no_background_subtraction_mqc.png" &&
    mv $hclustering_pos_nobg_om "Hirerchical_clustering_for_positive_no_background_subtraction_matched_only_mqc.png" &&
    mv $hclustering_neg_nobg "Hirerchical_clustering_for_negative_no_background_subtraction_mqc.png" &&
    mv $hclustering_neg_nobg_om "Hirerchical_clustering_for_negative_no_background_subtraction_matched_only_mqc.png" &&
    mv $hclustering_pos_withbg "Hirerchical_clustering_for_positive_with_background_subtraction_mqc.png" &&
    mv $hclustering_pos_withbg_om "Hirerchical_clustering_for_positive_with_background_subtraction_matched_only_mqc.png" &&
    mv $hclustering_neg_withbg "Hirerchical_clustering_for_negative_with_background_subtraction_mqc.png" &&
    mv $hclustering_neg_withbg_om "Hirerchical_clustering_for_negative_with_background_subtraction_matched_only_mqc.png" &&
    mv $vd_pos_nobg "Venn_diagram_for_positive_no_background_subtraction_mqc.png" &&
    mv $vd_neg_nobg "Venn_diagram_for_negative_no_background_subtraction_mqc.png" &&
    mv $vd_pos_withbg "Venn_diagram_for_positive_with_background_subtraction_mqc.png" &&
    mv $vd_neg_withbg "Venn_diagram_for_negative_with_background_subtraction_mqc.png" &&
    mv $barplot_pos_nobg "Bar_plot_clustering_for_positive_no_background_subtraction_mqc.png" &&
    mv $barplot_pos_nobg_om "Bar_plot_clustering_for_positive_no_background_subtraction_matched_only_mqc.png" &&
    mv $barplot_neg_nobg "Bar_plot_clustering_for_negative_no_background_subtraction_mqc.png" &&
    mv $barplot_neg_nobg_om "Bar_plot_clustering_for_negative_no_background_subtraction_matched_only_mqc.png" &&
    mv $barplot_pos_withbg "Bar_plot_clustering_for_positive_with_background_subtraction_mqc.png" &&
    mv $barplot_pos_withbg_om "Bar_plot_clustering_for_positive_with_background_subtraction_matched_only_mqc.png" &&
    mv $barplot_neg_withbg "Bar_plot_clustering_for_negative_with_background_subtraction_mqc.png" &&
    mv $barplot_neg_withbg_om "Bar_plot_clustering_for_negative_with_background_subtraction_matched_only_mqc.png"
    """
}

// Process for running MultiQC and generating the report.
process report_generator {

    publishDir './results/mqc/', mode: 'copy'

    input:
//    file mqc_dir from MQC_DIR
    file pos_data_info_mqc from POS_DATA_INFO_MQC
    file neg_data_info_mqc from NEG_DATA_INFO_MQC
    file experiments_info from EXPERIMENTS_INFO
    file mqc_config from MQC_CONFIG
    file peak_number_comparison_mqc from PEAK_NUMBER_COMPARISON_MQC
    file '*' from MQC_FIGS

    output:
    file "multiqc_report.html" into MULTIQC_REPORT

    when:
    params.bs == "1"

    shell:
    """
    multiqc .
    """

}

process mummichog_report {

    publishDir './results/mummichog/', mode: 'copy'

    input:

    file python_mummichog_input_prepare from PYTHON_MUMMICHOG_INPUT_PREPARE
    file pos_data from POS_WITHBG_FOR_MUMMICHOG

    output:
    file "*" into MUMMICHOG_REPORT

    when:
    params.bs == "1"

    shell:
    """
    echo "generating mommichog report"
    python3 ${python_mummichog_input_prepare} -i ${pos_data} -o ${params.data_pos_withbg_mummichog} &&
    mommichog -f ${params.data_pos_withbg_mummichog} -o ${params.data_pos_withbg_mummichog_out} -c ${params.cutoff}
    """

}