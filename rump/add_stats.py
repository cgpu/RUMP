#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description : This code do basic statistical tests (i.e., student t-test, fold change,
              Benjamini-Hochberg false discovery rate adjustment) for peak table generated
              by MZmine-2.53
Copyright   : (c) LemasLab, 02/23/2020
Author      : Xinsong Du
License     : GNU GPL-v3.0 License
Maintainer  : xinsongdu@ufl.edu, manfiol@ufl.edu, djlemas@ufl.edu
Usage       : python add_stats.py -i $input_peak_table
                                  -d $design_file_location
                                  -o $output_peak_table
                                  -l $library_location
"""

import warnings
import logging
import logging.handlers
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s]: %(levelname)s: %(message)s')

warnings.filterwarnings('ignore')

def add_threshold(row, names):
    """Add threshold for blank subtraction algorithm.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        names: column names in the peak table of a certain group of samples

    # Returns:
        threshold value
    """

    value = np.mean(row[names]) + 3*np.std(row[names])
    return value if value > 0 else 5000

def blank_subtraction_flag(row, name_group, name_threshold, bar):
    """Blank subtraction function.

    Blank subtraction algorithm:
    - Calculate mean (mean_blank) and standard deviation (sd_blank)
      of peak intensions in blank samples.
    - Threshold ← mean_blank+3*sd_blank
    - If threshold <=0, then replace it with 5,000 (why 5,000?)
    - Calculate mean peak intension in fat (mean_fat), whole (mean_whole)
      and skim (mean_skim) samples.
    - ratio_fat ← (mean_fat-threshold)/threshold;
      ratio_whole ← (mean_whole-threshold)/threshold;
      ratio_skim ← (mean_skim-threshold)/threshold
    - If ratio_fat<self_defined_number (e.g. 100) and
      ratio_whole<self_defined_number and ratio_skim<self_defined_number,
      then drop the peak.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        name_group: name of the group.
        name_threshold: name of the threshold column.
        bar: bar value of blank subtraction algorithm.

    # Returns:
        If a certain peak of this group still exist after blank subtraction

    """
    return (np.mean(row[name_group]) - row[name_threshold])/row[name_threshold] > bar

# Judge whether certain peak intensity of a sample is 0 or not
def zero_intensity_flag(row, name_group):
    """Check if the mean intensity of certain group of samples is zero. If zero, then
    the metabolite is not existed in that material.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        name_group: name of the group.

    # Returns:
        True (the mean intensity is zero) or False (the mean intensity is not zero).
    """
    return np.mean(row[name_group]) <= 0

# Add p-value for student t-test between two groups of samples
def add_pvalue(row, left_names, right_names):
    """Add p value for two group comparison based on student t-test.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        left_names: column names in the peak table of the first group of samples.
        right_names: column names in the peak table of the second group of samples.

    # Returns:
        p value of student t-test
    """
    _, p = stats.ttest_ind(row[left_names], row[right_names])
    return p

# Add t-value for student t-test between two groups of samples
def add_tvalue(row, left_names, right_names):
    """Add t value for two group comparison based on student t-test.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        left_names: column names in the peak table of the first group of samples.
        right_names: column names in the peak table of the second group of samples.

    # Returns:
        t value of student t-test
    """
    t, _ = stats.ttest_ind(row[left_names], row[right_names])
    return t

# Add fold-change for the mean values of two groups of samples
def fold_change(row, left, right):
    """Add fold change value for two group comparison.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        left: column name in the peak table of the mean intensity of first group of samples.
        right: column name in the peak table of the mean intensity of second group of samples.

    # Returns:
        fold change value.
    """
    if row[right] == 0:
        return np.inf
    if row[left] == 0:
        return -np.inf
    result = row[left]/row[right]
    return result if result >= 1 else -1/result

# Absolute value of fold-change
def abs_fold_change(row, fold_change_column):
    """Add absolute fold change value for two group comparison.

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        fold_change_column: column name in the peak table of the fold change value.

    # Returns:
        absolute fold change value.
    """
    return abs(row[fold_change_column])

# Add ppm value for identified metabolites.
## The library search result produced by MZmine may exceed 5 ppm,
## so those beyond 5 ppm should be filtered out
def add_ppm(row, library_df):
    """Add part per million (ppm) value for library matching. The library matching done by
    MZmine may not follow the threshold strictly (i.e., when setting the ppm to 5, some
    metabolites with ppm of more than 5 may also appear in the peak table).

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        library_df: library dataframe.

    # Returns:
        ppm value of the matched metabolite in the row.
    """
    if pd.isnull(row['row identity (main ID)']):
        return None
    mzs = list(library_df[library_df.Name.str.strip() == row['row identity (main ID)']]['M/Z'])
    mz_observe = row["row m/z"]
    diff = []
    for mz in mzs:
        diff.append(abs(mz_observe - mz))
    if len(diff) == 0:
        return None
    mz_theoretical = mzs[diff.index(min(diff))]
    return abs((mz_observe-mz_theoretical)*10e5/mz_theoretical)

def add_label(row, group1_name, group2_name):
    """Add label for metabolite represented by the row.
    Format: "m_z/retention_time/fold_change".

    # Arguments:
        row: certain row of peak table (pandas dataframe).
        group1_name: name of the group of first group of samples.
        group2_name: name of the group of second group of samples.

    # Returns:
        label (string type).
    """
    if pd.isnull(row["row identity (main ID)"]) or \
       row["row identity (main ID)"] == "nan" or \
       row["row identity (main ID)"] == None:
        return str(round(row["row m/z"], 2)) + "/" + \
               str(round(row["row retention time"], 2)) + \
               "/" + str(round(row["fold_change" + \
               "(" + str(group1_name) + " versus " + \
               str(group2_name) + ")"], 2))

    return str(row["row identity (main ID)"]) + "/" + \
           str(round(row["fold_change" + "(" + str(group1_name) + \
           " versus " + str(group2_name) + ")"], 2))

def add_stats(data_file="data_pos_ph.csv", design_file="design", \
              output_file="pos_withstats.csv", \
              library="positive_library.csv"):
    """Add basic statistics to peak table produced by MZmine.

    # Arguments:
        data_file: peak table.
        design_file: design file corresponding to the peak table.
        output_file: the name of processed file.
        library: location and file name of the library used to identify metabolites.

    # Returns:
        list of identified metabolites.

    # Outputs:
        prosessed peak table
    """

    data = pd.read_csv(data_file)
    data["row identity (main ID)"] = data["row identity (main ID)"].apply(str)
    data = data[~(data["row identity (main ID)"].str.contains("adduct|Complex", \
                na=False, regex=True))]
    data["number of comparisons"] = len(data)

    data_library = pd.read_csv(library)
    data["ppm"] = data.apply(lambda row: add_ppm(row, data_library), axis=1)

    blank_group_name = "zero-blank"
    design = pd.read_csv(design_file)
    ratio_bar = 100

    group_names = list(set(design['group']))
    group_names.sort()

    group1_name = group_names[0]
    group2_name = group_names[1]

    group1_columns = design[design.group == group1_name].sampleID.tolist()
    group2_columns = design[design.group == group2_name].sampleID.tolist()

#    group1_columns = list(set(group1_columns) & set(data.columns))
#    group2_columns = list(set(group2_columns) & set(data.columns))

    data[str(group1_name) + '_mean'] = data[group1_columns].mean(axis=1)
    data[str(group2_name) + '_mean'] = data[group2_columns].mean(axis=1)

    logger.info("calculating fold change")

    data['fold_change' + '(' + str(group1_name) + ' versus ' + str(group2_name) + ')'] = \
        data.apply(lambda row: \
        fold_change(row, str(group1_name) + '_mean', str(group2_name) + '_mean'), axis=1)
    data['log2_fold_change' + '(' + str(group1_name) + ' versus ' + str(group2_name) + ')'] = \
        np.log2(data[str(group1_name) + '_mean']/data[str(group2_name) + '_mean'])
    data['abs_fold_change' + "(" + str(group1_name) + " versus " + str(group2_name) + ")"] = \
        data.apply(lambda row: abs_fold_change(row, \
        'fold_change' + "(" + str(group1_name) + " versus " + str(group2_name) + ")"), axis=1)

    logger.info("calculating t-test p-value")

    data['p_value'] = data.apply(lambda row: \
        add_pvalue(row, group1_columns, group2_columns), axis=1)

    data['t_value'] = data.apply(lambda row: \
        add_tvalue(row, group1_columns, group2_columns), axis=1)

    data.dropna(subset=["p_value"], inplace=True)

    data[str(group1_name) + '_zero'] = data.apply(lambda row: \
        zero_intensity_flag(row, group1_columns), axis=1)

    data[str(group2_name) + '_zero'] = data.apply(lambda row: \
        zero_intensity_flag(row, group2_columns), axis=1)

    data['label'] = data.apply(lambda row: add_label(row, group1_name, group2_name), axis=1)

    _, adjust_p = multipletests(pvals=data.p_value.tolist(), alpha=0.05, method="fdr_bh")[:2]

    data["adjusted_p_value"] = adjust_p

    if blank_group_name in group_names:
        blank_columns = design[design.group == blank_group_name].sampleID.tolist()
        data['threshold'] = data.apply(lambda row: add_threshold(row, blank_columns), axis=1)

        data[str(group1_name) + "_selected"] = data.apply(lambda row: \
            blank_subtraction_flag(row, group1_columns, "threshold", ratio_bar), axis=1)

        data[str(group2_name) + "_selected"] = data.apply(lambda row: \
            blank_subtraction_flag(row, group2_columns, "threshold", ratio_bar), axis=1)

    data.to_csv(output_file, index=False)

    # the following return is used for unit tests
    return list(data["row identity (main ID)"])

if __name__ == '__main__':

    logger.info('generating venn diagram...')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', help="define the location of input csv file;", \
        default="data_pos_ph.csv", required=False)
    parser.add_argument(
        '-d', '--design', help="define the location of input design csv file;", \
        default="pos_design.csv", dest="design", required=False)
    parser.add_argument(
        '-o', '--output', help="define the location of output figure;", \
        default="pos_withstats.csv", required=False)
    parser.add_argument(
        '-l', '--library', help="define the location of library file;", \
        default="library.csv", dest="library", required=False)

    args = parser.parse_args()
    add_stats(args.input, args.design, args.output, args.library)
