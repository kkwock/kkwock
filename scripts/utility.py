# Import Modules
from datetime import datetime
from multiprocessing import Pool
import subprocess
from subprocess import Popen
from functools import reduce
import json
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings
import traceback
import string
import itertools
from itertools import chain
from pathlib import Path
import configparser

warnings.filterwarnings("ignore")


# Utility Functions

def check_vital(boolean_expression, message, input_header_list):
    '''
    Behaves like assert(), but it calls sys.exit() if the expression is false
    Prints a message to stdout in this case
    '''
    if not boolean_expression:
        traceback.print_stack()
        sys.exit("Error: {}".format(message))


def check_ts_vitals(input_filename, input_header_list):
    check_vital(str('compactPeakTable.csv') in str(input_filename),
                'Input file name must contain \'compactPeakTable.csv\' and be a compact peak table csv file.',
                input_header_list)
    check_vital('Well' in input_header_list, 'Input File Format Incorrect - Well Column Missing.',
                input_header_list)
    check_vital('Size [bp]' in input_header_list, 'Input File Format Incorrect - Size [bp] Column Missing.',
                input_header_list)
    check_vital('Peak Molarity [nmol/l]' in input_header_list,
                'Input File Format Incorrect - Peak Molarity [nmol/l] Column Missing.', input_header_list)


def filter_peaks(ts_data, config, mode):
    # -- Filter Peaks -- #
    # Will filter ts_data for peaks of interest (POI's) based on configs
    # PCR quants will filter for peaks based on bp range.
    # NTC quants will filter for peaks based on bp range & concentration threshold

    peaks_found = 0  # sets peaks found to false
    peak_dict = pd.DataFrame()  # will contain only peaks of interest

    # Config Input
    if mode == 'pcr':
        pcr = config['pcr']
        bp_min = pcr['bp_min']
        bp_max = pcr['bp_max']
    elif mode == 'ntc':
        ntc = config['ntc']
        bp_min = ntc['bp_min']
        bp_max = ntc['bp_max']
        conc = ntc['conc']
        conc_operator = ntc['conc_operator']

    # Parse for Peaks based for PCR or NTC
    for peak in range(0, len(ts_data.Well)):
        if mode == 'pcr':
            if ts_data.sizebp[peak] >= int(bp_min) and ts_data.sizebp[peak] <= int(bp_max):
                peak_dict = peak_dict.append(ts_data.loc[peak])
                peaks_found = 1
                ts_data.peak_presence[peak] = ts_data.peak_molarity[peak]
            else:
                ts_data.peak_presence[peak] = 0.00
        elif mode == 'ntc':
            for peak in range(0, len(ts_data.Well)):
                if int(bp_min) <= ts_data.sizebp[peak] <= int(bp_max) and \
                        pd.eval(f"{ts_data.concentration[peak]}{conc_operator}{conc}"):
                    peak_dict = peak_dict.append(ts_data.loc[peak])
                    peaks_found = 1
                    ts_data.peak_presence[peak] = ts_data.peak_molarity[peak]
                else:
                    ts_data.peak_presence[peak] = 0.00
    ts_data.reset_index(inplace=True, drop=True)  # resets indexes to use in loops later

    if peaks_found == 0:
        peak_wells = list('NA')
    else:
        peak_wells = list(peak_dict.Well)

    # Remove duplicate wells
    filtered_ts = ts_data[['Well', 'peak_presence']].drop_duplicates()
    filtered_ts.reset_index(inplace=True, drop=True)

    # Remove Electronic Ladder (EL) and Duplicate Peaks
    for i in range(0, len(filtered_ts.Well)):
        if filtered_ts.Well[i].startswith('EL'):
            filtered_ts = filtered_ts.drop(i, axis=0)
    filtered_ts.reset_index(inplace=True, drop=True)

    # Removes duplicate POI data
    for i in range(0, len(filtered_ts.Well)):
        if filtered_ts.Well[i] in peak_wells and filtered_ts.peak_presence[i] == 0.00:
            filtered_ts = filtered_ts.drop(i, axis=0)
    filtered_ts.reset_index(inplace=True, drop=True)

    # create matrix
    hdf = pd.DataFrame()
    ts_matrix = filtered_ts.groupby(['Well']).agg({'peak_presence': np.sum}).reset_index()

    # creates a 2D array from ts_data df
    for index in range(0, len(ts_matrix.Well)):
        row_char = str(ts_matrix.Well[index][0])
        col_num = int(ts_matrix.Well[index][1:])
        hdf.loc[row_char, col_num] = ts_matrix.peak_presence[index]

    # Fill NA wells with 0
    hdf = hdf.sort_index(axis=1, ascending=True)
    hdf = hdf.fillna(0)

    # asserts matrix ranges
    assert len(hdf) != 0

    return (filtered_ts, hdf)


def calc_sum(input_values):
    return np.sum(input_values)


def read_config(file):
    # embedded modules
    import os
    import configparser
    from pathlib import Path

    # path variables
    filename = Path(os.path.realpath(file))
    path = Path(filename).parent
    config = configparser.ConfigParser()

    if 'config' in file:
        config_dir = os.path.join(path, 'configs')

        if ("/configfile.ini" in file):
            config.read(file)
        else:
            config.read(f"{config_dir}/configfile.ini")

    elif 'results' in file:
        config_dir = os.path.join(path, 'src')
        if ("results.ini" in file):
            config.read(file)
        else:
            config.read(f"{config_dir}/results.ini")

    return config


def ts_parse(tapestation_csv_input, config, mode):
    csv = pd.read_csv(tapestation_csv_input, encoding='latin1')
    header_list = list(csv.columns)

    # input vitals check
    check_ts_vitals(tapestation_csv_input, header_list)

    # prepare TS data
    ts_data = csv.rename(columns={"Size [bp]": "sizebp",
                                  "Calibrated Conc. [ng/µl]": "concentration",
                                  "Peak Molarity [nmol/l]": "peak_molarity"})

    ts_data[["peak_presence"]] = ''

    # filter peaks
    filtered_ts, hdf = filter_peaks(ts_data, config, mode)

    if mode == 'pcr':
        pass_fail, test_iqr = passing_criteria_check(hdf, config, mode)
        return filtered_ts, hdf, pass_fail, test_iqr

    elif mode == 'ntc':
        pass_fail, ref_grade = passing_criteria_check(hdf, config, mode)

        return filtered_ts, hdf, pass_fail, ref_grade


def passing_criteria_check(hdf, config, mode):
    # Passing Criteria Check
    # input hdf - well matrix of tapestation data
    # config - configParsed config.ini
    # mode - either pcr or ntc for TS data
    # output: pass_fail

    if mode == 'pcr':
        ref = list(range(6, 8))
        test = list(range(6))
        pcr = config['pcr']
        iqr = pcr['iqr_min']
        iqr_operator = pcr['iqr_operator']
    elif mode == 'ntc':
        ref = [0]
        test = list(range(1, 7))

    # Ref and Test Data
    ref_data = (hdf[hdf.columns[ref]]).values.tolist()
    test_data = (hdf[hdf.columns[test]]).values.tolist()

    ref_list = []
    test_list = []

    # Flatten Datasets
    for peak in ref_data:
        # appending elements to the flat_list
        ref_list += peak

    for peak in test_data:
        # appending elements to the flat_list
        test_list += peak

    # NTC and PCR Passing Criteria
    pass_fail = ''
    ntc_peaks = 0
    test_grade = ''
    ref_grade = ''
    test_iqr = ''

    if mode == 'pcr':
        test_iqr = np.quantile(test_list, 0.25).round(4)
        if pd.eval(f"{test_iqr} {iqr_operator} {iqr}"):
            pass_fail = 'Pass'
        else:
            pass_fail = 'Fail'
        return pass_fail, test_iqr

    elif mode == 'ntc':
        if 0 in ref_list:
            ref_grade = 'Fail'
        for i in test_list:
            if i > 0:
                ntc_peaks += 1
        if ntc_peaks > 0:
            test_grade = 'Fail'
        if ref_grade or test_grade == 'Fail':
            pass_fail = 'Fail'
        else:
            pass_fail = 'Pass'
        return pass_fail, ref_grade


def ts_output(filtered_ts, config, mode, job_id):
    # matrix sums up all the peak values
    ts_matrix = filtered_ts.groupby(['Well']).agg({'peak_presence': np.sum}).reset_index()

    # Variables
    if mode == 'pcr':
        total_columns = list(range(1, 9))
        ref_columns = list(range(7, 9))
        output_df = pd.DataFrame(index=range(0, len(filtered_ts)),
                                 columns=['Job',
                                          'Well',
                                          'Sample Description',
                                          'Peak Molarity [nM]'])
    elif mode == 'ntc':
        total_columns = list(range(1, 8))
        ref_columns = [1]
        output_df = pd.DataFrame(index=range(0, len(filtered_ts)),
                                 columns=['Job',
                                          'Well',
                                          'Sample Description',
                                          'Peak Molarity [nM]',
                                          'Pass/Fail'])

    row_labels = list(string.ascii_uppercase[:8])
    ref_iter = list(itertools.product(row_labels, ref_columns))
    total_iter = list(itertools.product(row_labels, total_columns))
    ref_well = []
    total_well = []

    for a, b in ref_iter:
        ref_well.append(f"{a}{b}")

    for a, b in total_iter:
        total_well.append(f"{a}{b}")

    # Checker for incorrect well numbers
    missing_wells = list(set(list(total_well)) - set(list(ts_matrix['Well'])))

    ntc_ref_list = []
    ntc_test_list = []
    for peak in range(0, len(ts_matrix)):
        output_df['Well'][peak] = ts_matrix.Well[peak]
        output_df['Peak Molarity [nM]'][peak] = f"{ts_matrix.peak_presence[peak]}"
        output_df['Job'][peak] = job_id

        if output_df['Well'][peak] in ref_well:
            output_df['Sample Description'][peak] = 'ref'
        else:
            output_df['Sample Description'][peak] = 'test'

        # Failure for NTC
        if mode == 'ntc':
            if output_df['Sample Description'][peak] == 'ref':
                if ts_matrix.peak_presence[peak] > 0.00:
                    output_df['Pass/Fail'][peak] = 'Pass'
                else:
                    output_df['Pass/Fail'][peak] = 'Fail'
                    ntc_ref_list.append(output_df['Well'][peak])
            elif output_df['Sample Description'][peak] == 'test':
                if ts_matrix.peak_presence[peak] != 0.00:
                    output_df['Pass/Fail'][peak] = 'Fail'
                    ntc_test_list.append(output_df['Well'][peak])
                else:
                    output_df['Pass/Fail'][peak] = 'Pass'

    output_df['new'] = output_df['Well'].str.extract('(\d+)').astype(int)
    output_df = output_df.sort_values(by=['new', 'Well']).drop('new', axis=1)

    if mode == 'pcr':
        return (output_df)
    if mode == 'ntc':
        return (output_df, ntc_ref_list, ntc_test_list)


def pcr_heatmap(hdf, test_iqr, config, pass_fail, job_id, output_filepath):
    # Creates heatmap visual for TS data
    # Variables
    pcr = config['pcr']
    bp_min = pcr['bp_min']
    bp_max = pcr['bp_max']
    iqr_min = pcr['iqr_min']
    iqr_operator = pcr['iqr_operator']
    row_labels = list(string.ascii_uppercase[0:8])

    ax = plt.subplots(figsize=(12, 8))
    ax = sns.heatmap(hdf, annot=True, fmt=".2f", annot_kws={"size": 10}, cmap='Reds', vmin=0, vmax=150,
                     yticklabels=row_labels)
    ax.set_ylim(top=8, bottom=.01)
    ax.invert_yaxis()
    ax.set_yticklabels(labels=row_labels, rotation=0)
    plt.title(f"G360 EIO iQC - PCR Productivity - Heatmap \n"
              f"Peak Molarity nM ({bp_min}-{bp_max}bp)\n"
              f"{job_id}")
    plt.text(s=f'IQR Threshold: {iqr_operator}{iqr_min}\n'
               f'Lower IQR: {test_iqr}\n'
               f'Pass/Fail: {pass_fail}',
             ha='left', va='bottom', fontdict={'fontsize': 9}, y=9, x=0)

    plt.savefig(f'{output_filepath}{job_id}_pcr_heatmap.png')


def pcrprod_generate_figures(heatmap_df, test_qr_value, OUTPUT_FILE_PATH, JOBID):
    f, ax = plt.subplots(figsize=[15, 8])
    ax = sns.scatterplot(y='peak_presence', x='Well', hue='type', data=heatmap_df)
    plt.plot([0, 48], [100, 100], 'k-', lw=2)
    plt.xticks(rotation=45, size=9.5)
    plt.ylim(0, 250)
    plt.ylabel('Average Peak Molarity nM (130-165bp)')
    plt.title(
        'G360 EIO IQC - PCR Productivity [Standard]: TapeStation - Scatterplot \nAverage Peak Molarity nM (130-165bp)\n{0}'.format(
            JOBID))
    ROW_LABELS = ["A", "B", "C", "D", "E", "F", "G", "H"]
    if test_qr_value < 100:
        plt.text(s='FAIL,{0} - {1},Lower IQR > 100nM, Lower IQR = {2}\n'.format(130, 165, test_qr_value), ha='left',
                 va='bottom', fontdict={'fontsize': 9}, y=1, x=0)
    else:
        plt.text(s='PASS,{0} - {1},Lower IQR > 100nM, Lower IQR = {2}\n'.format(130, 165, test_qr_value), ha='left',
                 va='bottom', fontdict={'fontsize': 9}, y=1, x=0)

    plt.savefig(OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_scatterplot.png'.format(JOBID))

    hdf = pd.DataFrame()

    for index in range(0, len(heatmap_df.Well)):
        individual_row_char = str(heatmap_df.Well[index][0])
        individual_col_int = int(heatmap_df.Well[index][1:])
        hdf.loc[individual_row_char, individual_col_int] = float(heatmap_df.peak_presence[index])

    # hdf.reset_index(inplace=True, drop=True)
    # hdf.sort(columns=list(range(0,12)))
    hdf = hdf.sort_index(axis=1, ascending=True)

    # asserts matrix ranges
    assert len(hdf) != 0
    assert len(hdf.count(axis='rows')) * len(hdf.count(axis='columns')) <= 96

    ax1 = plt.subplots(figsize=(12, 8))
    ax1 = sns.heatmap(hdf, annot=True, fmt=".2f", annot_kws={"size": 10}, cmap='Reds', vmin=0, vmax=150,
                      yticklabels=ROW_LABELS)  # heatmap uses hdf 2D dataframe to populate
    ax1.set_ylim(top=8, bottom=.01)
    ax1.invert_yaxis()
    ax1.set_yticklabels(labels=ROW_LABELS, rotation=0)
    plt.title(
        'G360 EIO IQC - PCR Productivity [Standard] : TapeStation - Heatmap \nAverage Peak Molarity nM (130-165bp)\n{0}'.format(
            JOBID))

    if test_qr_value < 100:
        plt.text(s='FAIL,{0} - {1},Lower IQR > 100nM, Lower IQR = {2}\n'.format(130, 165, test_qr_value), ha='left',
                 va='bottom', fontdict={'fontsize': 9}, y=9, x=0)
    else:
        plt.text(s='PASS,{0} - {1},Lower IQR > 100nM, Lower IQR = {2}\n'.format(130, 165, test_qr_value), ha='left',
                 va='bottom', fontdict={'fontsize': 9}, y=9, x=0)
    plt.savefig(OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_heatmap.png'.format(JOBID))

    # print('URL to the PCR Productivity Output Scatterplot: https://bifs.ghdna.io' + OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_scatterplot.png\n'.format(JOBID))

    # print('URL to the PCR Productivity Output Heatmap: https://bifs.ghdna.io'+ OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_heatmap.png\n'.format(JOBID))


# Sequencing Step Functions (index_oligo_qc.py)

def read_oligo(oligo_file, directory, fcid):
    sample_sheet_file = directory + fcid + '_SampleSheet.csv'
    config = pd.read_csv(oligo_file)

    if 'G' in config['Sample_ID'][0]:
        config.loc[:, 'Sample_Project'] = fcid + '_G360_EIO_QC'
        assay = "G"
    else:
        if 'O' in config['Sample_ID'][0]:
            config.loc[:, 'Sample_Project'] = fcid + '_OMNI_LIO_QC'
            assay = "O"
        else:
            if 'L01' in config['Sample_ID'][0]:
                config.loc[:, 'Sample_Project'] = fcid + '_LUNAR1_EIO_QC'
                assay = "L01"
            else:
                config.loc[:, 'Sample_Project'] = fcid + '_LUNAR2_EIO_QC'
                assay = "L02"
    config.to_csv(sample_sheet_file, index=False)

    return sample_sheet_file, assay


def run_configuration(fcid_path, assay):
    cmd = "grep -i CompletedRead1Cycles " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    read1_cycles = o.decode("utf-8")
    read1_cycles = read1_cycles.split('>')[1]
    read1_cycles = read1_cycles.split('<')[0]

    cmd = "grep -i CompletedRead2Cycles " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    read2_cycles = o.decode("utf-8")
    read2_cycles = read2_cycles.split('>')[1]
    read2_cycles = read2_cycles.split('<')[0]

    cmd = "grep -i CompletedIndex1ReadCycles " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    index1_cycles = o.decode("utf-8")
    index1_cycles = index1_cycles.split('>')[1]
    index1_cycles = index1_cycles.split('<')[0]

    cmd = "grep -i CompletedIndex2ReadCycles " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    index2_cycles = o.decode("utf-8")
    index2_cycles = index2_cycles.split('>')[1]
    index2_cycles = index2_cycles.split('<')[0]

    if assay == 'G':
        if (read1_cycles == '26') & (read2_cycles == '0') & (index1_cycles == '7'):  # & (index2_cycles == '0'):
            return "Correct"
        else:
            return "Incorrect"
    else:
        if assay == '0':
            if (read1_cycles == '26') & (read2_cycles == '0') & (index1_cycles == '9') & (index2_cycles == '9'):
                return "Correct"
            else:
                return "Incorrect"
        else:
            if (assay == 'L01'):
                if (read1_cycles == '26') & (read2_cycles == '0') & (index1_cycles == '7') & (index2_cycles == '7'):
                    return "Correct"
                else:
                    return "Incorrect"
            else:
                if (read1_cycles == '26') & (read2_cycles == '0') & (index1_cycles == '10') & (index2_cycles == '10'):
                    return "Correct"
                else:
                    return "Incorrect"


def seq_metrics(fcid_path):
    cmd = "grep -i ClusterDensity " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    cluster_density = o.decode("utf-8")
    cluster_density = cluster_density.split('>')[1]
    cluster_density = cluster_density.split('<')[0]

    cmd = "grep -i ClustersPassingFilter " + fcid_path + "/RunCompletionStatus.xml"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    pass_filter = o.decode("utf-8")
    pass_filter = pass_filter.split('>')[1]
    pass_filter = pass_filter.split('<')[0]

    return int(round(float(cluster_density), 0)), int(round(float(pass_filter), 0))


def find_flowcell(fcid, directory, assay, paths):
    # paths = ['/ghds/ivd/raw', '/ghds/raw']
    for path in paths:
        cmd = "ls -d " + path + "/*" + fcid
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        o, e = ps.communicate()
        output = o.decode("utf-8")
        fcid_path = output.split('\n')[0]
        if "No such file" in fcid_path:
            next
        else:
            if run_configuration(fcid_path, assay) == 'Incorrect':
                return ('Incorrect')
            else:
                return fcid_path
    return (0)


def check_fcid_path(fcid_path, fcid):
    if fcid_path == 0:
        print("FCID " + fcid + " could not be found")
        return None
    if fcid_path == "Incorrect":
        print("Incorrect flowcell run configuration")
        return None


def seq_pass_fail(directory, oligo_file, fcid_path, config, url, job_id):
    # Adapted from initial_report()
    working_directory = os.getcwd()
    if directory.startswith(working_directory):
        output_directory = directory
    else:
        output_directory = os.path.join(working_directory, directory)

    # Variables
    seq_outpath = f"{directory}{job_id}_g360-eioqc_seq.csv"
    output_directory = url + output_directory

    # -- Config Thresholds--#
    seq = config['seq']

    # Cluster Density
    cd_min_cutoff_val = str(seq['cd_min'])
    cd_min_cutoff_op = str(seq['cd_min_operator'])
    cd_max_cutoff_val = str(seq['cd_max'])
    cd_max_cutoff_op = str(seq['cd_max_operator'])

    # Cluster Passing Filter
    pf_cutoff_val = str(seq['pf_cutoff'])
    pf_cutoff_op = str(seq['pf_operator'])

    oligo = pd.read_csv(oligo_file)
    with open(seq_outpath, 'w') as fout:

        ### Sequencing Quality Metrics Section
        cluster_density, cluster_pf = seq_metrics(fcid_path)

        # This evaluates if cluster density is within thresholds
        if pd.eval(f'{cluster_density}{cd_min_cutoff_op}{cd_min_cutoff_val}'):
            cluster_density_status = 'Pass'
        else:
            cluster_density_status = 'Fail'

        # This evaluates if cluster pass filter is >= 80
        if pd.eval(f'{cluster_pf}{pf_cutoff_op}{pf_cutoff_val}'):
            cluster_pf_status = "Pass"
        else:
            cluster_pf_status = "Fail"

        if (cluster_density_status == "Pass") & (cluster_pf_status == "Pass"):
            seq_qual_status = "Pass"
        else:
            seq_qual_status = "Fail"

    return (seq_qual_status, cluster_density, cluster_pf, directory + job_id + '_g360-eioqc_seq.csv')


def results_config(repo_path, mode, seq_input, overall_seq_passfail, seq_passfail,
                   seq_cd, seq_cf, xc_passfail, xc_list,
                   pcr_input, pcr_passfail, pcr_iqr, ntc_input,
                   ntc_passfail, ntc_ctrl_passfail, ntc_ref_list, ntc_test_list):
    results = configparser.ConfigParser()
    results_path = f'{repo_path}/results.ini'
    results.read(results_path)

    if mode not in results:
        results.add_section(f"{mode}")

    if mode == 'seq':
        results.set('seq', 'seq_input', f'{seq_input}')
        results.set('seq', 'overall_seq_passfail', f'{overall_seq_passfail}')
        results.set('seq', 'xc_passfail', f'{xc_passfail}')
        results.set('seq', 'xc_list', f'{xc_list}')
        results.set('seq', 'seq_passfail', f'{seq_passfail}')
        results.set('seq', 'seq_cd', f'{seq_cd}')
        results.set('seq', 'seq_cf', f'{seq_cf}')
    elif mode == 'pcr':
        results.set('pcr', 'pcr_input', f'{pcr_input}')
        results.set('pcr', 'pcr_passfail', f'{pcr_passfail}')
        results.set('pcr', 'pcr_iqr', f'{pcr_iqr}')
    elif mode == 'ntc':
        results.set('ntc', 'ntc_input', f'{ntc_input}')
        results.set('ntc', 'ntc_passfail', f'{ntc_passfail}')
        results.set('ntc', 'ntc_ctrl_passfail', f'{ntc_ctrl_passfail}')
        results.set('ntc', 'ntc_ref_list', f'{ntc_ref_list}')
        results.set('ntc', 'ntc_test_list', f'{ntc_test_list}')

    with open(f"{repo_path}/results.ini", 'w') as result_file:
        results.write(result_file)


def initial_report(directory, oligo_file, fcid, fcid_path, config, url, job_id):
    working_directory = os.getcwd()
    if directory.startswith(working_directory):
        output_directory = directory
    else:
        output_directory = os.path.join(working_directory, directory)

    output_directory = url + output_directory
    seq = config['seq']

    ########## Get acceptance criteria ############
    cd_min_cutoff_val = str(seq['cd_min'])
    cd_min_cutoff_op = str(seq['cd_min_operator'])
    cd_max_cutoff_val = str(seq['cd_max'])
    cd_max_cutoff_op = str(seq['cd_max_operator'])
    pf_cutoff_val = str(seq['pf_cutoff'])
    pf_cutoff_op = str(seq['pf_operator'])
    ###### Finish getting acceptance criteria #####

    oligo = pd.read_csv(oligo_file)
    with open(directory + job_id + '_g360-eioqc_seq.csv', 'w') as fout:
        fout.write('Purpose:,' + ' '.join(oligo['Sample_Project'].iloc[0].split('_')[1:]) + ' Validation\n')
        fout.write('URL Output Files:,' + output_directory + '\n')
        fout.write('Flowcell ID:,' + fcid + '\n')
        fout.write('Script Name:,g360_index_oligo_qc.py' + script_name + '\n')
        fout.write('Script Version:,1.0' + version + '\n')
        fout.write('Location of input files:,' + fcid_path + '\n')
        fout.write(
            'Location of report:,' + output_directory + oligo['Sample_Project'].iloc[0] + '_g360-eioqc_seq.csv' + '\n')
        fout.write('Date of analysis:,' + str(datetime.today().strftime("%m/%d/%Y")) + "\n")

        ### Sequencing Quality Metrics Section
        cluster_density, cluster_pf = seq_metrics(fcid_path)
        fout.write('\nSequencing Quality Metrics:\n')
        fout.write('Metrics,Passing Criteria,Value      ,Pass/Fail\n')

        # This evaluates if cluster density is within thresholds
        if eval(str(cluster_density) + " " + cd_min_cutoff_op + " " + cd_min_cutoff_val):
            cluster_density_status = 'Pass'
        else:
            cluster_density_status = 'Fail'

        fout.write('Cluster Density,' + cd_min_cutoff_op + cd_min_cutoff_val + " K/mm^2," + str(
            cluster_density) + " K/mm^2," + cluster_density_status + '\n')

        # This evaluates if cluster pass filter is >= 80
        if eval(str(cluster_pf) + " " + pf_cutoff_op + " " + pf_cutoff_val):
            cluster_pf_status = "Pass"
        else:
            cluster_pf_status = "Fail"

        fout.write('Cluster Passing Filter,' + pf_cutoff_op + pf_cutoff_val + "%," + str(
            cluster_pf) + "%," + cluster_pf_status + '\n')

        if (cluster_density_status == "Pass") & (cluster_pf_status == "Pass"):
            seq_qual_status = "Pass"
        else:
            seq_qual_status = "Fail"

    return (seq_qual_status, directory + job_id + 'seq.csv')


def generate_sample_sheet(directory, input_dir, oligo_csv, name):
    '''Used to generate sample sheet from Xiaobo's spreadsheet.

    Parameters
    ----------
    oligo_tsv : str, location of tab separated spreadsheet of Xiabo's oligo sequences
    name : str, prefix to name output file

    Output
    ------
    Sample sheet in csv format

    Returns
    -------
    str, name of sample sheet
    '''

    # outputs a sample sheet that is unused - >

    # input columns: Sample_ID, Sequence, i7_Index, i5_Index, Index_ID, Sample_Project
    oligos = pd.read_csv(oligo_csv)
    oligos = oligos.drop_duplicates(subset=['Index_ID'])

    samp_sheet = pd.DataFrame(index=oligos.index)

    samp_sheet.loc[:, 'Lane'] = 1
    samp_sheet.loc[:, 'Sample_ID'] = oligos['Index_ID']
    samp_sheet.loc[:, 'Sample_Name'] = oligos['Index_ID']

    samp_sheet.loc[:, 'Sample_Plate'] = ''
    samp_sheet.loc[:, 'Sample_Well'] = ''
    samp_sheet.loc[:, 'I7_Index_ID'] = ''

    samp_sheet.loc[:, 'index'] = oligos['i7_Index']

    mode = 'single'

    if 'i5_Index' in oligos.columns:
        mode = 'dual'
        samp_sheet.loc[:, 'index2'] = oligos['i5_Index']
        samp_sheet.loc[:, 'I5_Index_ID'] = ''

    samp_sheet.loc[:, 'Sample_Project'] = oligos['Sample_Project']
    samp_sheet.loc[:, 'Description'] = ''

    samp_sheet = samp_sheet[~samp_sheet['index'].isna()]
    samp_sheet.to_csv(directory + name + '.tmp', index=False)

    # print(directory+name+'.csv')

    with open(directory + name + '.tmp') as f:
        with open(directory + name + '.csv', 'w') as fout:
            fout.write('[Data]' + ',' * (len(samp_sheet.columns) - 1) + '\n')
            for line in f:
                fout.write(line)

    os.remove(directory + name + '.tmp')

    return directory + name + '.csv', oligos['Sample_Project'].iloc[0], mode


def call_bcl2fastq(sample_csv, in_dir, out_dir, threads=1, mismatches=0):
    '''
    Calls bcl2fastq with no mismatches and enforcing a read length of 20 (default)
    '''

    file_out = out_dir + "bcl2fastq.log"

    args = f"bcl2fastq -R {in_dir} " \
           f"-o {out_dir} " \
           f"--sample-sheet {sample_csv} " \
           f"-p {threads} " \
           f"--barcode-mismatches {mismatches} " \
           f"--ignore-missing-bcls"

    # print('Demuxing...\n')

    with open(file_out, "w+") as f:
        p = Popen(args.split(' '), stdout=f, stderr=f)
        code = p.wait()

    return code


def count_matches_sub(params):
    '''Count exact matches to each oligo in fastq.gz files

    Parameters
    ----------
    fastq : str, location of gzipped fastq file
    df : pd.DataFrame, zero filled with index values as index and fastq file names as columns
        This parameter is only used to make sure the index is consistent between all series
    seq_df : pd.DataFrame, dataframe read from oligo tsv file

    Returns
    -------
    fq_ID : str, beggining of fastq file name
    s : pd.Series, counts for that fastq file for each oligo

    '''

    fastq, seq_df = params
    fq_ID = '_'.join(fastq.split('/')[-1].split('_')[0:2])

    # open fastq file, get only seq line + up to 25th char, sort, get count of uniques, sort decreasing reads and write out to tsv file
    tsv_file = fastq.split(".")[0] + ".tsv"
    cmd = "zcat " + fastq + " | awk 'NR%4==2' | cut -c 1-25 | sort | uniq -c | sort -nr > " + tsv_file
    ps = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # removes leading whitespaces
    cmd = "sed -i 's/^ *//' " + tsv_file
    ps = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # create dataframe from the tsv file called reads_df
    reads_df = pd.read_csv(tsv_file, sep=" ", names=['Count', 'Sequence'])

    seq_df_select = seq_df[['Sample_ID', 'Sequence']]

    # merge based on perfect match to oligo
    df_merge = pd.merge(left=reads_df, right=seq_df_select, left_on='Sequence', right_on='Sequence', how='right')
    df_merge_select = df_merge[['Sample_ID', 'Count']]
    df_merge_select = df_merge_select.sort_values(by='Sample_ID')
    df_merge_select.rename(columns={'Count': str(fq_ID)}, inplace=True)
    s = df_merge_select.fillna(0)
    s = s[s['Sample_ID'] != 0]
    return s


def count_mixed_dual(fastq, oligo_df):
    '''Count exact matches to each oligo in Undertermined fastq.gz file.
    These are the reads that didn't demux bc not perfect macth to both i7 and i5 index.
    These reads will have a mixture of different i7 and i5 indices

    Parameters
    ----------
    fastq : str, location of Undetermined fastq file
    oligo_df : pd.DataFrame, dataframe read from oligo tsv file

    Returns
    -------
    DATAFRAME THAT TABULATES ALL I7 AND I5 PERMUTATION

    '''

    # create tsv file with the indexes from undetermined fastq file
    index_tsv_file = fastq.split(".")[0] + "_index.tsv"
    cmd = "zcat " + fastq + " | awk 'NR%4==1' | cut -d ':' -f 10  > " + index_tsv_file
    ps = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    index_df = pd.read_csv(index_tsv_file, sep="+", names=['I7_Index_ID', 'I5_Index_ID'])

    # create tsv file with the reads from undetermined fastq file
    seq_tsv_file = fastq.split(".")[0] + "_seq.tsv"
    cmd = "zcat " + fastq + " | awk 'NR%4==2' | cut -c 1-25 > " + seq_tsv_file
    ps = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    seq_df = pd.read_csv(seq_tsv_file, names=['Sequence'])

    # merge the two tsv files
    mixed_df = pd.merge(index_df, seq_df, left_index=True, right_index=True)
    mixed_df = mixed_df.groupby(['I7_Index_ID', 'I5_Index_ID', 'Sequence']).size().reset_index(name='Count')

    # merge the mixed_df with the known oligos
    seq_df_select = oligo_df[['Sample_ID', 'Sequence']]
    seq_df_select = seq_df_select.drop_duplicates()
    mixed_df = pd.merge(left=mixed_df, right=seq_df_select, left_on='Sequence', right_on='Sequence', how='right')

    # merge the mixed_df with the known i7
    seq_df_select = oligo_df[['Index_ID', 'i7_Index']]
    seq_df_select = seq_df_select.drop_duplicates()
    mixed_df = pd.merge(left=mixed_df, right=seq_df_select, left_on='I7_Index_ID', right_on='i7_Index', how='right')

    # merge the mixed_df with the known i5
    seq_df_select = oligo_df[['Index_ID', 'i5_Index']]
    seq_df_select = seq_df_select.drop_duplicates()
    mixed_df = pd.merge(left=mixed_df, right=seq_df_select, left_on='I5_Index_ID', right_on='i5_Index', how='right')

    mixed_df = mixed_df[['Index_ID_x', 'Index_ID_y', 'Sample_ID', 'Count']]
    mixed_df.rename(columns={'Index_ID_x': 'i7_Index', 'Index_ID_y': 'i5_Index', 'Sample_ID': 'Oligo'}, inplace=True)

    # write out the list of oligos that matched and which combo of i7 and i5 indexes it was amplified with
    mixed_df.to_csv(fastq.split(".")[0] + "_mixed_index.tsv", index=False)

    mixed_df['i7_i5'] = mixed_df['i7_Index'] + ':' + mixed_df['i5_Index']
    mixed_select = mixed_df[['i7_i5', 'Oligo', 'Count']]
    mixed_select = mixed_select.pivot_table(index=['Oligo'], columns='i7_i5', values='Count', fill_value=0)

    return (mixed_select)


def count_matches(directory, config_file, mode, threads=6):
    '''Multiprocessing setup for count_matches_sub, requires oligo tsv in Xiaobo's format

    Parameters
    ----------
    directory : str, location of fastq files
    config_file : str, location of tsv file containing oligos and expected indexes
    mode : single or dual indices
    threads : int, default 6, number of processors to Pool

    Returns
    -------
    df : pd.DataFrame containing read counts for each oligo for each fastq file (index)

    Output
    ------
    *read_counts.csv dataframe
    '''
    seq_df = pd.read_csv(config_file)
    project_name = seq_df['Sample_Project'].iloc[0]

    fq_list = [x for x in os.listdir(directory + project_name) if x.endswith('fastq.gz')]

    df = pd.DataFrame()
    params = [(directory + project_name + '/' + x, seq_df) for x in fq_list]
    p = Pool(threads)
    s_list = p.map(count_matches_sub, params)

    tsv_list = [x for x in os.listdir(directory + project_name) if x.endswith('.tsv')]

    if len(fq_list) != len(tsv_list):
        print("Analysis files not created. Please check permissions")
        exit

    df = reduce(lambda x, y: pd.merge(x, y, on='Sample_ID'), s_list)

    # add in eio columns if there was no fastq file for that expected eio
    eios = seq_df['Index_ID'].drop_duplicates()
    for i in eios:
        if i in df.columns:
            continue
        else:
            df[i] = 0

    # Merge the fastq read count df + undetermined fastq df
    if mode != 'single':
        mixed_dual = count_mixed_dual(directory + '/Undetermined_S0_L001_R1_001.fastq.gz', seq_df)
        mixed_dual = mixed_dual.reset_index(drop=False)
        df = pd.merge(left=df, right=mixed_dual, left_on='Sample_ID', right_on='Oligo', how='right')
        df = df.drop('Oligo', axis=1)
        # print('Finished processing undetermined fastq file')

    df.to_csv(directory + project_name + '_read_counts.csv', index=False)

    return df


def calculate_contamination(directory, read_counts_df, job_id, config):
    # Config Variables
    seq = config['seq']
    assay = seq['assay']
    slope = seq['slope']
    intercept = seq['intercept']
    contam_cutoff = seq['contam_cutoff']
    contam_operator = seq['contam_operator']

    #
    m = float(slope)
    b = float(intercept)
    df = read_counts_df.copy()
    # Get the sum of reads for each oligo = OLIGO_SUM column
    df['OLIGO_SUM'] = df.sum(axis=1, skipna=True)

    # Split the Sample_ID column
    df['TYPE'] = df['Sample_ID'].str.split('_').str[0]
    df['SET'] = df['Sample_ID'].str.split('_').str[2]
    df['EIO'] = df['Sample_ID'].str.split('_').str[3]

    # Melt dataframe and get percentage of reads
    df = df.drop(columns=["Sample_ID"])
    df = pd.melt(df, id_vars=["TYPE", "SET", "EIO", "OLIGO_SUM"])
    df['i7'] = df['variable'].str.split(':').str[0].str.split('_').str[1]
    df['perc_reads'] = 100 * df['value'] / df['OLIGO_SUM']

    # Only keep testing oligos (TYPE == T)
    df = df[df['TYPE'] == "T"]

    # Remove if both i7 and i5 are equal to the oligo
    if assay != 'G':
        df['i5'] = df['variable'].str.split(':').str[1].str.split('_').str[1]
        df.loc[df['i5'].isnull(), 'i5'] = df['i7']
        df = df[(df['EIO'] != df['i7']) | (df['EIO'] != df['i5'])]
        df = df[['TYPE', 'SET', 'EIO', 'i7', 'i5', 'value', 'OLIGO_SUM', 'perc_reads']]
    else:
        df = df[(df['EIO'] != df['i7'])]
        df = df[['TYPE', 'SET', 'EIO', 'i7', 'value', 'OLIGO_SUM', 'perc_reads']]

    df.perc_reads = df.perc_reads.round(2)

    # do not apply correction factor if the perc_reads is 0
    df['predicted_perc_contam'] = np.where(df['perc_reads'] > 0,
                                           df['perc_reads'].subtract(b).divide(m),
                                           df['perc_reads'])

    # round the percentage of predicted reads to two decimal places
    df.predicted_perc_contam = df.predicted_perc_contam.round(2)

    # write out contamination csv file
    df.rename(columns={'value': 'oligo_reads', 'OLIGO_SUM': 'total_oligo_reads'}, inplace=True)
    df.to_csv(f'{directory}{job_id}_test.csv', index=False)

    # Use only the testing oligos
    # df = df.drop(columns=["oligo_reads", "total_oligo_reads", "perc_reads", "TYPE"])
    df = df.drop(columns=["TYPE"])

    if assay == 'G':
        # flag known index pairs that will have eio bleeding for G360:
        df = df.pivot_table(index=['EIO', 'i7'],
                            columns='SET',
                            values='predicted_perc_contam',
                            fill_value=0)
        df['Median'] = df.median(axis=1)
        df['Job ID'] = job_id
        df['Contamination'] = ''
        df = df.reset_index(drop=False)
        df['Bleed_Pair'] = np.where((((df['i7'] == '26') & (df['EIO'] == '14')) |
                                     ((df['i7'] == '09') & (df['EIO'] == "41")) |
                                     ((df['i7'] == '16') & (df['EIO'] == "23"))),
                                    'TRUE', 'FALSE')
        eio_bleeding_df = df[df['Bleed_Pair'] == "TRUE"]
        df = df[df['Bleed_Pair'] != "TRUE"]
        df = df.append(eio_bleeding_df)

        # Insert Contamination Column
        for val in range(len(df)):
            if pd.eval(f"{df.Median[val]}{contam_operator}{contam_cutoff}"):
                if df.Bleed_Pair[val] != 'TRUE':
                    df['Contamination'][val] = f"EIO{df['i7'][val]} in EIO{df['EIO'][val]}"

        df = df.set_index(['EIO', 'i7'])
    else:
        if assay != 'O':
            df = df.pivot_table(index=['EIO', 'i7', 'i5'],
                                columns='SET',
                                values='predicted_perc_contam',
                                fill_value=0)
            # use mad?
            df['Median'] = df.median(axis=1)
        else:
            # REP will be either A, B, C, D
            df['REP'] = df['SET'].apply(lambda x: x[-1])
            # SET will be either Set1, Set2, Set3
            df['SET'] = df['SET'].apply(lambda x: x[:-1])
            # Pivot table
            df = df.pivot_table(index=['i7', 'i5', 'EIO', 'REP'], columns='SET', values='predicted_perc_contam',
                                fill_value=0)
            df['Median'] = df.median(axis=1)

    df = df.reset_index(drop=False)
    df = df.rename_axis(None, axis=1)
    df = df[['Job ID', 'Bleed_Pair', 'EIO', 'i7', 'Set1', 'Set2', 'Set3', 'Median', 'Pass/Fail']]
    df = df.rename(columns={"Set1": "Rep1", "Set2": "Rep2", "Set3": "Rep3", "Median": "Median_PC"})
    df.index.names = ['ix']
    df = df.sort_values('Bleed_Pair', ascending=False)

    df.to_csv(directory + str(job_id) + '_contamination.csv', index=True, index_label='ix')
    return df


def count_oligos(directory, job_id, read_counts_df):
    # Returns dataframe of counted oligos.
    df = read_counts_df.copy()

    # Remove Set4
    df = df.loc[df["Sample_ID"].str.contains("Set4") == False]

    # Sum each set
    split_cols = df["Sample_ID"].str.split("_", n=-1, expand=True)
    df['Sample_ID'] = "EIO_" + split_cols[3]
    df = df.groupby(df["Sample_ID"], as_index=False).sum()

    # Merge EIOs
    df['Count'] = [df.iloc[i][d] for i, d in enumerate(df['Sample_ID'])]
    df = df[['Sample_ID', 'Count']]
    df['Total_Reads'] = sum(df['Count'])
    df['Read_Frac'] = pd.eval(df['Count'] / df['Total_Reads'])

    df.to_csv(directory + str(job_id) + '_read-fractions.csv', index=True, index_label='ix')

    return(df)


def contam_pass_fail(directory, read_counts_df, job_id, config):
    # Config Variables
    seq = config['seq']
    assay = seq['assay']
    slope = seq['slope']
    intercept = seq['intercept']
    contam_cutoff = seq['contam_cutoff']
    contam_operator = seq['contam_operator']

    #
    m = float(slope)
    b = float(intercept)
    df = read_counts_df.copy()
    # Get the sum of reads for each oligo = OLIGO_SUM column
    df['OLIGO_SUM'] = df.sum(axis=1, skipna=True)

    # Split the Sample_ID column
    df['TYPE'] = df['Sample_ID'].str.split('_').str[0]
    df['SET'] = df['Sample_ID'].str.split('_').str[2]
    df['EIO'] = df['Sample_ID'].str.split('_').str[3]

    # Melt dataframe and get percentage of reads
    df = df.drop(columns=["Sample_ID"])
    df = pd.melt(df, id_vars=["TYPE", "SET", "EIO", "OLIGO_SUM"])
    df['i7'] = df['variable'].str.split(':').str[0].str.split('_').str[1]
    df['perc_reads'] = 100 * df['value'] / df['OLIGO_SUM']

    # Only keep testing oligos (TYPE == T)
    df = df[df['TYPE'] == "T"]
    df = df[(df['EIO'] != df['i7'])]
    df = df[['TYPE', 'SET', 'EIO', 'i7', 'value', 'OLIGO_SUM', 'perc_reads']]
    df.perc_reads = df.perc_reads.round(2)

    # do not apply correction factor if the perc_reads is 0
    df['predicted_perc_contam'] = np.where(df['perc_reads'] > 0,
                                           df['perc_reads'].subtract(b).divide(m),
                                           df['perc_reads'])

    df.predicted_perc_contam = df.predicted_perc_contam.round(2)

    # write out per sample pair oligo reads
    df.rename(columns={'value': 'oligo_reads', 'OLIGO_SUM': 'total_oligo_reads'}, inplace=True)
    # df.to_csv(f'{directory}{job_id}_oligo-reads.csv', index=False)


    # Use only the testing oligos
    df = df.drop(columns=["oligo_reads", "total_oligo_reads", "perc_reads", "TYPE"])

    # Build DF and Flag Bleeding Pairs
    df = df.pivot_table(index=['EIO', 'i7'],
                        columns='SET',
                        values='predicted_perc_contam',
                        fill_value=0)
    df['Median'] = df.median(axis=1)
    df['Job ID'] = job_id
    df['Contamination'] = ''
    df = df.reset_index(drop=False)
    df['Bleed_Pair'] = np.where((((df['i7'] == '26') & (df['EIO'] == '14')) |
                                 ((df['i7'] == '09') & (df['EIO'] == "41")) |
                                 ((df['i7'] == '16') & (df['EIO'] == "23"))),
                                'TRUE', 'FALSE')
    eio_bleeding_df = df[df['Bleed_Pair'] == "TRUE"]
    df = df[df['Bleed_Pair'] != "TRUE"]
    df = df.append(eio_bleeding_df)

    # Insert Contamination Column
    for val in range(len(df)):
        if pd.eval(f"{df.Median[val]}{contam_operator}{contam_cutoff}"):
            if df.Bleed_Pair[val] != 'TRUE':
                df['Contamination'][val] = f"EIO{df['i7'][val]} in EIO{df['EIO'][val]}"

    contam_count = 0
    contam_list = []
    for a in df['Contamination']:
        if a != '':
            contam_list.append(a)
            contam_count = + 1

    if contam_count > 0:
        contam_pass_fail = 'Fail'
    else:
        contam_pass_fail = 'Pass'

    df = df.set_index(['EIO', 'i7'])
    df = df.reset_index(drop=False)
    df = df.rename_axis(None, axis=1)
    df = df[['Job ID', 'Bleed_Pair', 'EIO', 'i7', 'Set1', 'Set2', 'Set3', 'Median', 'Contamination']]
    df = df.rename(columns={"Set1": "Rep1", "Set2": "Rep2", "Set3": "Rep3", "Median": "Median_PC"})
    df.index.names = ['ix']
    df = df.sort_values('Bleed_Pair', ascending=False)
    df.to_csv(directory + str(job_id) + '_contamination.csv', index=True, index_label='ix')

    return df, contam_pass_fail, contam_list


def calculate_performance(directory, project_name, read_counts_df):
    df = read_counts_df.copy()
    df['SET'] = df['Sample_ID'].str.split('_').str[2]
    df = df.set_index('Sample_ID')
    # remove columns that have a ":" in them; this is for dropping columns w/ mixed i7 and i5
    df = df[df.columns.drop(list(df.filter(regex=":")))]
    total = df.drop('SET', axis=1).sum().sum()

    # compile the proportion of reads for each index in each set
    performance_df = pd.DataFrame(columns=['SET', 'Index', 'perc_reads'])
    for s in df['SET'].drop_duplicates().sort_values():
        for n in df.index.str.split('_').str[-1].drop_duplicates():
            n_df = df[(df['SET'] == s) & (df.index.str.endswith(n))]
            # perc reads per set per index out of total reads with assigned oligo
            pr = n_df.drop('SET', axis=1).sum().sum() / float(total) * 100
            if pr > 0:
                row = pd.Series([s, n, pr], index=performance_df.columns)
                performance_df = performance_df.append(row, ignore_index=True)

    # calculate mean, std dev and z scores
    control_indexes = df[df.index.str.startswith('R')].index
    # need to see all 16 controls
    if len(control_indexes) != 16:
        performance_df = []
        return (performance_df)
    avg = performance_df.loc[performance_df.Index.isin(control_indexes.str[-2:]), 'perc_reads'].mean()
    std = performance_df.loc[performance_df.Index.isin(control_indexes.str[-2:]), 'perc_reads'].std(ddof=0)
    performance_df.loc[:, 'z_score'] = performance_df['perc_reads'].subtract(avg).divide(std)
    performance_df.loc[:, 'perc_of_mean'] = performance_df['perc_reads'].divide(avg).multiply(100)
    performance_df.to_csv(directory + project_name + '_performance.csv')

    performance_df = performance_df.drop(columns=["z_score", "perc_reads"])
    performance_df = performance_df.pivot_table(index=['Index'], columns='SET', values='perc_of_mean', fill_value=0)
    performance_df = performance_df.drop(columns=["Set4"])
    performance_df['Median'] = performance_df.median(axis=1)
    performance_df = performance_df.reset_index(drop=False)

    performance_df = performance_df.rename_axis(None, axis=1)
    performance_df.index.names = ['ix']
    performance_df.to_csv(directory + project_name + '_performance_set.csv', index=True, index_label='ix')

    return performance_df


def process_barcodes(directory):
    with open(directory + 'Stats/Stats.json') as f:
        run_stats = json.loads(f.read())
    RPF = run_stats['ConversionResults'][0]['TotalClustersPF']

    reads_per_exp = []
    for s in run_stats['ConversionResults'][0]['DemuxResults']:
        reads_per_exp.append(s['NumberReads'] / float(RPF) * 100)

    unexp_ix = []
    for bc, reads in run_stats['UnknownBarcodes'][0]['Barcodes'].items():
        pct = reads / float(RPF) * 100
        if pct > min(reads_per_exp):
            unexp_ix.append((bc, pct))
    return unexp_ix


def print_contam(cont_df):
    with open(report_path, 'a') as fout:
        fout.write('\nEIO Cross-contamination quality metrics:\n')

        ### Cross Contamination Metrics Section (MK-TAG LEVEL / SAMPLE LEVEL)
        df_i7 = cont_df.set_index('i7')
        for index, new_df in df_i7.groupby(level=0):
            max_row = new_df[new_df['Median'] == new_df['Median'].max()]
            i7 = 'Sample_' + str(index)
            oligo = 'Sample_' + str((max_row.iloc[0]['OLIGO']))
            median = max_row.iloc[0]['Median']
            if eval(str(median) + " " + contam_cutoff_op + " " + contam_cutoff_val):
                pass
            else:
                failed_sample_list.append(oligo)
                cont_sample_fail += 1

            ##################################################################################################
            ### Cross Contamination Metrics Section (INDEX LEVEL)
            # fout.write('\nCross contamination quality metrics (EIO Contamination):\n')
            fout.write('Index,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
            df_oligo = cont_df.set_index('OLIGO')
            for index, new_df in df_oligo.groupby(level=0):
                max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                oligo = 'EIO_' + str(index)
                i7 = 'EIO_' + str((max_row.iloc[0]['i7']))
                median = max_row.iloc[0]['Median']
                if eval(str(median) + " " + contam_cutoff_op + " " + contam_cutoff_val):
                    fout.write(
                        oligo + "," + contam_cutoff_op + contam_cutoff_val + "%," + "%0.2f" % median + "%,Pass\n")
                else:
                    fout.write(
                        oligo + "," + contam_cutoff_op + contam_cutoff_val + "%," + "%0.2f" % median + "%,Fail," + i7 + "\n")
                    failed_i7_list.append(oligo)
                    cont_contam_fail += 1


def final_report(report_path, cont_df, performance_df, seq_qual_status, config):
    # Read Config Variables
    seq = config['seq']
    assay = seq['assay']

    ######## TYPE OF INDEX OLIGO  #########
    if assay == 'O':
        io_type = 'LIO'
    else:
        io_type = 'EIO'

    ######### GET ACCEPTANCE CRITERIA ##########
    contam_cutoff_val = seq['contam_cutoff']
    contam_cutoff_op = seq['contam_operator']
    perf_cutoff_val = seq['perf_cutoff']
    perf_cutoff_op = seq['perf_operator']

    #

    cont_contam_fail = 0
    cont_sample_fail = 0
    failed_sample_list = list()
    failed_i7_list = list()

    if assay == 'G':
        # Separate known eio bleeding pairs from the other pairs
        eio_bleed = cont_df[cont_df['eio_bleeding_pair'] == "TRUE"]
        cont_df = cont_df[cont_df['eio_bleeding_pair'] != "TRUE"]

    # Writing out report
    with open(report_path, 'a') as fout:
        fout.write('\nEIO Cross-contamination quality metrics:\n')

        ### Cross Contamination Metrics Section (MK-TAG LEVEL / SAMPLE LEVEL)
        # fout.write('\nCross contamination quality metrics (Sample Contamination):\n')
        if assay == 'G':
            # fout.write('Sample,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
            df_i7 = cont_df.set_index('i7')
            for index, new_df in df_i7.groupby(level=0):
                max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                i7 = 'Sample_' + str(index)
                oligo = 'Sample_' + str((max_row.iloc[0]['OLIGO']))
                median = max_row.iloc[0]['Median']
                if eval(str(median) + " " + contam_cutoff_op + " " + contam_cutoff_val):
                    # fout.write(i7 + "," + contam_cutoff_op +  contam_cutoff_val + "%," + "%0.2f" %median + "%,Pass\n")
                    pass
                else:
                    # fout.write(i7 + "," + contam_cutoff_op +  contam_cutoff_val + "%," + "%0.2f" %median + "%,Fail," + oligo + "\n")
                    failed_sample_list.append(oligo)
                    cont_sample_fail += 1

        else:
            # for dual indices, only look at when i7 == i5, these will be the misassigned reads --> translates to mk-tag contamination
            df_dual = cont_df[cont_df['i7'] == cont_df['i5']]
            df_dual = df_dual.set_index('i7')

            if assay == 'O':
                rep = cont_df['REP'].unique()
            else:
                rep = ['A']

            for x in rep:
                if len(rep) == 1:
                    df_dual_rep = df_dual
                else:
                    df_dual_rep = df_dual[df_dual['REP'] == x]
                    # fout.write(io_type + ' Replicate ' + x + ' Sample Contamination\n')

                # fout.write('Sample,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
                for index, new_df in df_dual_rep.groupby(level=0):
                    if index in df_dual_rep[
                        'OLIGO'].unique():  # Only print out for results of Samples that did have an Oligo
                        max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                        i7 = 'Sample_' + str(index)
                        oligo = 'Sample_' + str((max_row.iloc[0]['OLIGO']))
                        median = max_row.iloc[0]['Median']
                        if eval(str(median) + " " + contam_cutoff_op + " " + contam_cutoff_val):
                            # fout.write(i7 + "," + contam_cutoff_op +  contam_cutoff_val + "%," + "%0.2f" %median + "%,Pass\n")
                            pass
                        else:
                            # fout.write(i7 + "," + contam_cutoff_op +  contam_cutoff_val + "%," + "%0.2f" %median + "%,Fail," + oligo + "\n")
                            cont_sample_fail += 1

        ##################################################################################################
        ### Cross Contamination Metrics Section (INDEX LEVEL)
        if assay == 'G':
            # fout.write('\nCross contamination quality metrics (EIO Contamination):\n')
            fout.write('Index,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
            df_oligo = cont_df.set_index('OLIGO')
            for index, new_df in df_oligo.groupby(level=0):
                max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                oligo = 'EIO_' + str(index)
                i7 = 'EIO_' + str((max_row.iloc[0]['i7']))
                median = max_row.iloc[0]['Median']
                if eval(str(median) + " " + contam_cutoff_op + " " + contam_cutoff_val):
                    fout.write(
                        oligo + "," + contam_cutoff_op + contam_cutoff_val + "%," + "%0.2f" % median + "%,Pass\n")
                else:
                    fout.write(
                        oligo + "," + contam_cutoff_op + contam_cutoff_val + "%," + "%0.2f" % median + "%,Fail," + i7 + "\n")
                    failed_i7_list.append(oligo)
                    cont_contam_fail += 1

            ### Known EIO Bleeding (only applies to G360 EIOs)
            fout.write('\nCross contamination quality metrics (Known EIO-bleeding. Informational Only):\n')
            fout.write('Index,Passing Criteria,Value,Source\n')
            for index, row in eio_bleed.iterrows():
                oligo = 'EIO_' + row['OLIGO']
                eio = 'EIO_' + row['i7']
                fout.write(eio + ',None,' + "%0.2f" % row['Median'] + "%," + oligo + "\n")
        else:
            # For dual indices, print out 2 sections (1) i7-level contamination and (2) i5-level contamination

            # i7-level contamination:
            fout.write('\nCross contamination quality metrics (' + io_type + ' Contamination: i7 Contamination):\n')
            df_oligo_i7 = cont_df[cont_df['i7'] != cont_df['OLIGO']]
            df_oligo_i7 = df_oligo_i7.set_index('OLIGO')
            if assay == 'O':
                rep = cont_df['REP'].unique()
            else:
                rep = ['A']
            for x in rep:
                if len(rep) == 1:
                    df_oligo_i7_rep = df_oligo_i7
                else:
                    df_oligo_i7_rep = df_oligo_i7[df_oligo_i7['REP'] == x]
                    fout.write(io_type + ' Replicate ' + x + ' i7 Contamination\n')
                    fout.write('Index,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
                for index, new_df in df_oligo_i7_rep.groupby(level=0):
                    max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                    oligo = io_type + '_' + str(index)
                    i7 = io_type + '_' + str((max_row.iloc[0]['i7']))
                    median = max_row.iloc[0]['Median']
                    if eval(str(median) + " " + index_cutoff_op + " " + index_cutoff_val):
                        fout.write(
                            oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Pass\n")
                    else:
                        fout.write(
                            oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Fail," + i7 + "\n")
                        cont_contam_fail += 1

        #     # i5-level contamination:
        #     fout.write('\nCross contamination quality metrics (' + io_type + ' Contamination: i5 Contamination):\n')
        #     df_oligo_i5 = cont_df[cont_df['i5'] != cont_df['OLIGO']]
        #     df_oligo_i5 = df_oligo_i5.set_index('OLIGO')
        #     if assay == 'O':
        #         rep = cont_df['REP'].unique()
        #     else:
        #         rep = ['A']
        #     for x in rep:
        #         if len(rep) == 1:
        #             df_oligo_i5_rep = df_oligo_i5
        #         else:
        #             df_oligo_i5_rep = df_oligo_i5[df_oligo_i5['REP'] == x]
        #             fout.write(io_type + ' Replicate ' + x + ' i5 Contamination\n')
        #         fout.write('Index,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
        #         for index, new_df in df_oligo_i5_rep.groupby(level=0):
        #             max_row = new_df[new_df['Median'] == new_df['Median'].max()]
        #             oligo = io_type + '_' + str(index)
        #             i5 = io_type + '_' + str((max_row.iloc[0]['i5']))
        #             median = max_row.iloc[0]['Median']
        #             if eval(str(median) + " " + index_cutoff_op + " " + index_cutoff_val):
        #                 fout.write(
        #                     oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Pass\n")
        #             else:
        #                 fout.write(
        #                     oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Fail," + i5 + "\n")
        #                 cont_contam_fail += 1
        # '''
        ################################################################################################################################
        ### Concentration Metrics Section
        fout.write('\nPerformance quality metrics:\n')
        fout.write('Performance quality metrics (Concentration):\n')
        if len(performance_df) != 0:
            fout.write('Index,Passing Criteria,Value,Pass/Fail\n')
            performance_fail = 0
            for index, row in performance_df.iterrows():
                i7 = io_type + '_' + row['Index']
                median = row['Median']
                if eval(str(median) + " " + perf_cutoff_op + " " + perf_cutoff_val):
                    fout.write(
                        i7 + "," + perf_cutoff_op + perf_cutoff_val + '%,' + "%0.2f" % row['Median'] + "%,Pass\n")
                else:
                    fout.write(
                        i7 + "," + perf_cutoff_op + perf_cutoff_val + '%,' + "%0.2f" % row['Median'] + "%,Fail\n")
                    performance_fail += 1
        else:
            performance_fail = 'NA'
            fout.write('Reference oligos not detected in run\n')
        '''
        ################################################################################################################################################
        ### Summary Section
        # fout.write('\nCross-Contamination Qualification Summary:')
        # fout.write('\nQC Metric,Passing Criteria,Comments,Status\n')

        # ## Overall Status: Pass/Fail/Repeat
        # if seq_qual_status != "Pass":
        #     overall_status = "FAIL"
        # elif cont_contam_fail != 0:
        #     overall_status = "FAIL"
        # else:
        #     overall_status = "PASS"
        # 
        # if seq_qual_status != "Pass":
        #     seq_qual_status = "FAIL"
        # else:
        #     seq_qual_status = 'PASS'
        # 
        # if overall_status != "PASS" and seq_qual_status == "FAIL":
        #     fout.write('Overall Cross-Contamination QC,All sections pass,Sequence Qualification Fail,' + str(
        #         overall_status) + '\n')
        # elif overall_status != "PASS" and cont_contam_fail != 0:
        #     fout.write('Overall Cross-Contamination QC,All sections pass,Cross-Contamination Qualification Fail,' + str(
        #         overall_status) + '\n')
        # else:
        #     fout.write('Overall Cross-Contamination QC,All sections pass,NA,' + str(overall_status) + '\n')
        # ################################################################################################################################################
        # # fout.write('Sequencing Qualification Status:,Both Cluster Density and Cluster Passing Filter pass,,' + seq_qual_status + '\n')
        # fout.write('Sequencing Qualification,Both Cluster Density and Cluster Passing Filter pass,NA,' + str(
        #     seq_qual_status) + '\n')
        # if cont_contam_fail == 0:
        #     fout.write('Cross Contamination Qualification,All ' + io_type + 's pass,NA,PASS\n')
        # else:
        #     fout.write('Cross Contamination Qualification,All ' + io_type + 's pass,' + io_type + 's ' + str(
        #         failed_i7_list).replace(',', ' ').replace(']', '').replace('[', '') + ' failed,FAIL\n')

        if performance_fail == 0:
            fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,Pass\n')
        else:
            if performance_fail == 'NA':
                fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,NA\n')
            else:
                fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,Fail\n')

        '''
        # fout.write('\nPerformed By:\n')
        # fout.write('Performed Date:\n')

    fout.close()
    # print('\n\nURL to Cross-Contamination Output CSV: https://bifs.ghdna.io{0}\n'.format(report_path))
    # f = open(report_path)
    # for line in f:
    #    print(line)
    # f.close()

    return (seq_qual_status, cont_contam_fail, overall_status)


def return_report_path(output_path):
    # this function was created primarily to make it easy for SQA to test the report csv file
    # performs ls to find report file and returns it
    cmd = "ls " + output_path + "/*SEQ_RESULT.csv"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    o, e = ps.communicate()
    output = o.decode("utf-8")
    report_path = output.split('\n')[0]
    if "No such file" in report_path:
        return ('Could not find report')
    else:
        return (report_path)

    ### Write Output


def g360_generate_report(output_file, fcid, fc_path, config, results, lot, part, url):
    # Uses a results.ini file with results input to parse data and create report.

    SCRIPT_NAME = os.path.basename(__file__)
    SCRIPT_DIR = Path(__file__).parent

    # Version Number
    if os.path.exists(f"{SCRIPT_DIR}/../VERSION.txt"):
        with open(SCRIPT_DIR / "../VERSION.txt", "r") as file:
            result = file.readline().strip().split(" ")
            GIT_VERSION = result[0]
            if len(result) == 2:
                VERSION = result[1]
            else:
                VERSION = ""
    else:
        # EXTRACT MOST RECENT GIT COMMIT HASH FOR SCRIPT
        result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE)
        GIT_VERSION = result.stdout.decode("utf-8").strip()

        # EXRACT GIT TAG (SCRIPT VERSION) FOR MOST RECENT GIT COMMIT HASH
        result = subprocess.run(["git", "tag", "--list", "--contains", GIT_VERSION], stdout=subprocess.PIPE)
        tag = result.stdout.decode("utf-8").strip().split("/")[-1]
        result = subprocess.run(["git", "rev-list", "--count", "HEAD"], stdout=subprocess.PIPE)
        n_builds = result.stdout.decode("utf-8").strip()
        VERSION = tag + "-" + n_builds + "-" + GIT_VERSION

    now = datetime.now()
    TODAY_STRING = now.strftime("%m/%d/%Y")

    # ---RESULT VARIABLES---#
    # creating dictionary for threshold values and setting to NA
    r_dict = {}
    default_value = 'NA'
    config_list = ['seq_input', 'overall_seq_passfail', 'xc_passfail', 'xc_list', 'seq_passfail', 'seq_cd', 'seq_cf',
                   'pcr_input', 'pcr_passfail', 'pcr_iqr', 'ntc_input', 'ntc_passfail', 'ntc_ctrl_passfail']
    r_dict = dict.fromkeys(config_list, default_value)

    # Set Variables
    mode_list = ['pcr', 'ntc', 'seq']
    for mode in mode_list:
        if mode in results:
            if mode == 'seq':
                r_seq = results['seq']
                seq_input = r_dict['seq_input'] = r_seq['seq_input']
                r_dict['overall_seq_passfail'] = r_seq['overall_seq_passfail']
                r_dict['xc_passfail'] = r_seq['xc_passfail']
                r_dict['xc_list'] = r_seq['xc_list']
                r_dict['seq_passfail'] = r_seq['seq_passfail']
                r_dict['seq_cd'] = r_seq['seq_cd']
                r_dict['seq_cf'] = r_seq['seq_cf']
            if mode == 'pcr':
                r_pcr = results['pcr']
                pcr_input = r_dict['pcr_input'] = r_pcr['pcr_input']
                r_dict['pcr_passfail'] = r_pcr['pcr_passfail']
                r_dict['pcr_iqr'] = r_pcr['pcr_iqr']
            if mode == 'ntc':
                r_ntc = results['ntc']
                ntc_input = r_dict['ntc_input'] = r_ntc['ntc_input']
                r_dict['ntc_passfail'] = r_ntc['ntc_passfail']
                r_dict['ntc_ctrl_passfail'] = r_ntc['ntc_ctrl_passfail']
                r_dict['ntc_ref_list'] = r_ntc['ntc_ref_list']
                r_dict['ntc_test_list'] = r_ntc['ntc_test_list']

    # ---SEQUENCING METRICS---#
    seq = config['seq']
    contam_cutoff = str(seq['contam_cutoff'])
    contam_op = str(seq['contam_operator'])
    cd_min = str(seq['cd_min'])
    cd_min_op = str(seq['cd_min_operator'])
    cd_max = str(seq['cd_max'])
    cd_max_op = str(seq['cd_max_operator'])
    pf_cutoff = str(seq['pf_cutoff'])
    pf_op = str(seq['pf_operator'])

    # ---PCR METRICS---#
    pcr = config['pcr']
    iqr_min = pcr['iqr_min']
    iqr_op = pcr['iqr_operator']
    pcr_bp_min = pcr['bp_min']
    pcr_bp_max = pcr['bp_max']

    # ---NTC METRICS---#
    ntc = config['ntc']
    conc = ntc['conc']
    conc_op = ntc['conc_operator']
    ntc_bp_min = ntc['bp_min']
    ntc_bp_max = ntc['bp_max']

    # ---Build Results Dataframes---#
    qc_report = pd.DataFrame({"Metric": ["Cross-Contamination", "PCR Productivity", "NTC", "Overall Qualification"],
                              "Outcome": "NA"
                              })

    # ---SEQUENCING METRICS---#
    seq_df = pd.DataFrame({"Lot": lot,
                           "Metric": ["Cluster Density", "Cluster Passing Filter", "Cross-Contamination"],
                           "Threshold": [f"{cd_min_op}{cd_min}K/mm^2", f"{pf_op}{pf_cutoff}%",
                                         f"{contam_op}{contam_cutoff}%"],
                           "Value": [f"{r_dict['seq_cd']}", f"{r_dict['seq_cf']}", "NA"],
                           "Outcome": "NA",
                           "Comment": ""})

    if pd.eval(f"{r_dict['seq_cd']}{cd_min_op}{cd_min}"):
        seq_df.loc[0]['Outcome'] = 'Pass'
    else:
        seq_df.loc[0]['Outcome'] = 'Fail'

    if pd.eval(f"{r_dict['seq_cf']}{pf_op} {pf_cutoff}"):
        seq_df.loc[1]['Outcome'] = 'Pass'
    else:
        seq_df.loc[1]['Outcome'] = 'Fail'

    if r_dict['xc_passfail'] == 'Fail':
        xc_list = (str(r_dict['xc_list']).strip('[]').replace('\'', ''))
        seq_df.loc[2]['Comment'] = f"{xc_list}"

    seq_pass = seq_df['Outcome'] == 'Pass'
    if seq_pass.all():
        qc_report.loc[0]['Outcome'] = 'Pass'
    else:
        qc_report.loc[0]['Outcome'] = 'Fail'

    # ---PCR METRICS---#
    pcr_df = pd.DataFrame({"Lot": [lot],
                           "Metric": ["Lower iQR"],
                           "Threshold": [f'{iqr_op}{iqr_min}nM'],
                           "Value": f"{r_dict['pcr_iqr']}",
                           "Outcome": f"{r_dict['pcr_passfail']}",
                           "Comment": ""})

    if r_dict['pcr_passfail'] == 'Pass':
        pcr_df.loc[0]['Comment'] = 'All Wells Pass'
        qc_report.loc[1]['Outcome'] = 'Pass'
    else:
        qc_report.loc[1]['Outcome'] = 'Fail'

    # ---NTC METRICS---#
    ntc_df = pd.DataFrame({"Lot": [lot, ""],
                           "Metric": ["test", "ref"],
                           "Threshold": f"{ntc_bp_min}-{ntc_bp_max}bp,{conc_op}{conc}ng/uL",
                           "Condition": ["No peaks present", "No missing peaks"],
                           "Outcome": [f"{r_dict['ntc_passfail']}", f"{r_dict['ntc_ctrl_passfail']}"],
                           "Comment": ""})

    if r_dict['ntc_passfail'] == 'Pass':
        ntc_df.loc[0]['Comment'] == 'All Wells Pass'
    else:
        test_list = (str(r_dict['ntc_test_list']).strip('[]').replace('\'', ''))
        ntc_df.loc[0]['Comment'] = f"{test_list}"
    if r_dict['ntc_ctrl_passfail'] == 'Pass':
        ntc_df.loc[1]['Comment'] == 'All Wells Pass'
    else:
        ref_list = (str(r_dict['ntc_ref_list']).strip('[]').replace('\'', ''))
        ntc_df.loc[1]['Comment'] = f"{ref_list}"

    ntc_pass = ntc_df['Outcome'] == 'Pass'
    if ntc_pass.all():
        qc_report.loc[2]['Outcome'] = 'Pass'
    else:
        qc_report.loc[2]['Outcome'] = 'Fail'

    # ---OVERALL METRICS---#
    # Overall Qualification eioqc_passfail
    overall_pass = qc_report['Outcome'] == 'Pass'
    if overall_pass.all():
        qc_report.loc[3]['Outcome'] = 'Pass'
    else:
        qc_report.loc[3]['Outcome'] = 'Fail'

    header = f"Script name:,{SCRIPT_NAME}\nScript version:,{VERSION}\n" \
             f"Script github commit version:,{str(GIT_VERSION)}\n" \
             f"GH part number:,{part}\n" \
             f"Date of analysis:,{TODAY_STRING}\n" \
             f"Flowcell ID:, {fcid}\n" \
             f"Lot number:,{lot}\n" \
             f"Flowcell filepath:,{fc_path}\n" \
             f"PCR input filepath:,{pcr_input}\n" \
             f"NTC input filepath:,{ntc_input}\n" \
             f"Output filepath:,{output_file}\n\n" \
             f"Summary of qualification status:\n"
    footer = f"\nPerformed by:\nPerformed date:\nVerified by:\nVerified on:\n"

    # ---PRINT OUTPUTS---#
    with open(output_file, 'w') as fout:
        fout.write(header)
    qc_report.to_csv(output_file, sep=",", mode="a", index=False)

    # ---SEQUENCING OUTPUT---#
    with open(output_file, "a") as fout:
        fout.write('\nCross-contamination metric:\n')
    seq_df.to_csv(output_file, sep=",", mode="a", index=False)

    # ---PCR OUTPUT---#
    with open(output_file, "a") as fout:
        fout.write('\nPCR Productivity metric:\n')
    pcr_df.to_csv(output_file, sep=",", mode="a", index=False)

    # ---NTC OUTPUT---#
    with open(output_file, "a") as fout:
        fout.write('\nNTC metric:\n')
    ntc_df.to_csv(output_file, sep=",", mode="a", index=False)

    with open(output_file, "a") as fout:
        fout.write(footer)

    print(f"Output CSV: {output_file}")
    output_file = output_file.replace("^/ghess/", "/ghds/")
    output_url = f"{url}{output_file}"
    print(f"Output URL: {output_url}")

    return (0)


# Print Results
def print_results(output_df, output_filepath, job_id, mode):
    # output_filename
    if mode == 'ntc':
        report_name = 'g360-eioqc_ntc'
    elif mode == 'pcr':
        report_name = 'g360-eioqc_pcr'

    output_filename = (f'{output_filepath}{job_id}_{report_name}.csv')

    script_dir = Path(__file__).parent

    with open(output_filename, "w") as output_file:
        output_df.to_csv(output_file, sep=",", mode="a", index=False)
    output_file.close()


