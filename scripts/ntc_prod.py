# G360 EIO NTC Tool v1.0 [02/03/2021]
from datetime import datetime
import argparse
import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings
import traceback

warnings.filterwarnings("ignore")






PERFORMED_BY = ' '
PERFORMED_DATE = ' '

TODAY_STRING = datetime.now().strftime("%m/%d/%Y")
ROW_LABELS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']




# def below courtesy of Isaac Ho
def check_vital(boolean_expression, message, input_header_list):
    '''
    Behaves like assert(), but it calls sys.exit() if the expression is false
    Prints a message to stdout in this case
    '''
    if not boolean_expression:
        traceback.print_stack()
        sys.exit("Error: {}".format(message))


# def below adopted from code written by Isaac Ho
# formats information to output csv file
def print_results(output_df, col_header_list, pass_or_fail_string, input_file_name, url, OUTPUT_FILE_PATH, JOBID, \
                  LOT_NUMBER, GH_PART_NUMBER, INPUT_FILE_NAME, BP_MIN, BP_MAX):
    output_file_name = (OUTPUT_FILE_PATH + '{0}_G360EIOiQC_NTC_RESULT.csv'.format(JOBID))

    with open(output_file_name, "w") as output_file:

        output_file.write('Purpose:,G360 v2.7 EIO NTC Tool - QC Lot Release Data\n')
        output_file.write('Job ID:,{0}\n'.format(JOBID))
        output_file.write(
            'URL to Output File:,' + url + OUTPUT_FILE_PATH + '{0}_G360EIOiQC_NTC_RESULT.csv\n'.format(JOBID))
        output_file.write('Base Pair Range:,{0} - {1}\n'.format(BP_MIN, BP_MAX))
        output_file.write('Lot Number:,{0}\n'.format(LOT_NUMBER))
        output_file.write('Part Number:,{0}\n'.format(GH_PART_NUMBER))
        output_file.write('Script Name:,g360_eioqc_ntc.py\n')
        output_file.write('Script Version:,v1.0\n')
        output_file.write('Input File Name:,{0}\n'.format(input_file_name))
        output_file.write('Location of input files:,' + INPUT_FILE_NAME + '\n')
        output_file.write('Date Data Analyzed:,{0}\n'.format(TODAY_STRING))

        output_file.write('\nNTC Overall Qualification Summary:\n')
        output_file.write('Base Pair Range,Passing Criteria,Comments,Status\n')
        if pass_or_fail_string == 'Pass':
            output_file.write('{0} - {1},No unwanted product found in NTC area,NA,PASS\n\n'.format(BP_MIN, BP_MAX))
        else:
            output_file.write('{0} - {1},No unwanted product found in NTC area,{2},FAIL\n\n'.format(BP_MIN, BP_MAX,
                                                                                                    pass_or_fail_string))

        output_file.write('Peak Value Metrics\n')

        output_df.to_csv(output_file, index=False, header=col_header_list, line_terminator='\n')

        output_file.write('\nPerformed By:,{0}\n'.format(PERFORMED_BY))
        output_file.write('Performed Date:,{0}\n'.format(PERFORMED_DATE))

    output_file.close()

    # print('\n\nURL to NTC Output CSV: https://bifs.ghdna.io' + '{0}\n'.format(output_file_name))


# Checks over all wells after heatmap set is built to determine overall QC status
def passing_criteria_check(hdf):
    set_not_full = 0
    ntc_area = 0
    positive_control_col = 0
    pass_or_fail_string = str()

    # no pass if set has nan values in it
    for col in range(1, 8):  # all 7 cols
        for row in range(0, 8):  # all 8 rows
            try:
                if np.isnan(hdf[col][row]) == True:
                    set_not_full = set_not_full + 1
            except KeyError:
                set_not_full = set_not_full + 1

    # no pass if ntc area has a non-zero value
    for col in range(2, 8):
        for row in range(0, 8):
            try:
                if hdf[col][row] != 0.00:
                    ntc_area = ntc_area + 1
            except KeyError:
                ntc_area = ntc_area + 1

    # no pass if positive control row has nan or a zero value
    # always col 1
    for row in range(0, 8):
        try:
            if hdf[1][row] == 0.00 or np.isnan(hdf[1][row]) == True:
                positive_control_col = positive_control_col + 1
        except KeyError:
            positive_control_col = positive_control_col + 1

    if set_not_full == 0:
        if ntc_area == 0:
            if positive_control_col == 0:
                pass_or_fail_string = 'Pass'

    if set_not_full != 0:
        pass_or_fail_string = pass_or_fail_string + ' Data set not fully built.'
    if ntc_area != 0:
        pass_or_fail_string = pass_or_fail_string + ' NTC area has non-zero value(s) present.'
    if positive_control_col != 0:
        pass_or_fail_string = pass_or_fail_string + ' Positive Control columns has ' + \
                              'zero or NaN value present.'

    return pass_or_fail_string


def filter_peaks(ts_data, BP_MIN, BP_MAX):
    # sets peaks found to false
    peaks_found = 0
    peak_dict = pd.DataFrame()  # will contain only peaks of interest
    # puts value into 'peak_presence' column, based on user input for OUTPUT_UNITS
    for peak in range(0, len(ts_data.Well)):
        if ts_data.sizebp[peak] >= BP_MIN and ts_data.sizebp[peak] <= BP_MAX and ts_data.:
            # populates a dictionary with peak ts_data if within bp range
            peak_dict = peak_dict.append(ts_data.loc[peak])
            peaks_found = 1
            ts_data.peak_presence[peak] = ts_data.peak_molarity[peak]
        else:
            ts_data.peak_presence[peak] = 0.00
    ts_data.reset_index(inplace=True, drop=True)  # resets indexes to use in loops later

    return (ts_data, peaks_found, peak_dict)


def calc_sum(input_values):
    return np.sum(input_values)


def main():
    # user inputs TS file name including .csv extension
    parser = argparse.ArgumentParser()
    parser.add_argument('-LOT_NUMBER', type=str, required=True)
    parser.add_argument('-ntc_input_file', type=str, required=True)
    parser.add_argument('-FCID', type=str, required=True)
    parser.add_argument('-output_path', type=str, required=False,
                        default='/ghds/groups/reagent_qc/g360/eio_qc/output_files/')
    parser.add_argument('-GHPN', type=str, required=False, default='')
    parser.add_argument('-folder_id', type=str, required=False)
    parser.add_argument('-URL', type=str, required=False, default='https://bifs.ghdna.io')

    args = parser.parse_args()
    LOT_NUMBER = str(args.LOT_NUMBER)
    JOBID = str(args.FCID + '_' + args.LOT_NUMBER)
    input_file_name_w_ext = args.ntc_input_file
    GH_PART_NUMBER = args.GHPN
    folder_id = args.folder_id
    url = args.URL
    OUTPUT_FILE_PATH = args.output_path + folder_id + '/' + JOBID + '/'  # args.output_file_path
    INPUT_FILE_NAME = input_file_name_w_ext

    if not os.path.exists(OUTPUT_FILE_PATH):
        os.mkdir(OUTPUT_FILE_PATH)

    BP_MIN = int(100)  # args.user_input_bp_min
    BP_MAX = int(500)  # args.user_input_bp_max

    ts_data = pd.DataFrame()  # contains input csv file information
    hdf = pd.DataFrame()  # will contain filtered tapestation ts_data of all wells
    csv_input_file = pd.DataFrame()

    csv_input_file = pd.read_csv(INPUT_FILE_NAME, encoding='latin1')

    input_header_list = list(csv_input_file.columns)

    check_vital(str('compactPeakTable.csv') in str(input_file_name_w_ext),
                'Input file name must contain \'compactPeakTable.csv\' and be a compact peak table csv file.',
                input_header_list)
    check_vital('Well' in input_header_list, 'Input File Format Incorrect - Well Column Missing.', input_header_list)
    check_vital('Size [bp]' in input_header_list, 'Input File Format Incorrect - Size [bp] Column Missing.',
                input_header_list)
    check_vital('Peak Molarity [nmol/l]' in input_header_list,
                'Input File Format Incorrect - Peak Molarity [nmol/l] Column Missing.', input_header_list)

    # save ts_data file w new header names for easier use downstream
    ts_data = csv_input_file.rename(columns={"Size [bp]": "sizebp", "Peak Molarity [nmol/l]": "peak_molarity"})
    ts_data[["peak_presence"]] = ''

    ts_data, peaks_found, peak_dict = filter_peaks(ts_data, BP_MIN, BP_MAX)

    if peaks_found == 0:
        peak_oi_wells = list('NA')  # without this the program will error out if no POI'S are found
    else:
        peak_oi_wells = list(peak_dict.Well)  # list of wells where peaks of interest (POI) were found

    # trim ts_data frame so we only have the ts_data we want to use
    ts_data = ts_data[['Well', 'peak_presence']].drop_duplicates()
    ts_data.reset_index(inplace=True, drop=True)

    for i in range(0, len(ts_data.Well)):  # filters out the electronic ladder rows instnaces
        if ts_data.Well[i].startswith('EL'):
            ts_data = ts_data.drop(i, axis=0)
    ts_data.reset_index(inplace=True, drop=True)

    # POI wells have a duplicate instance at this point in the ts_data df the duplicates
    # are a false representation of the well and are trimmed from the ts_data df
    for i in range(0, len(ts_data.Well)):
        if ts_data.Well[i] in peak_oi_wells and ts_data.peak_presence[i] == 0.00:
            ts_data = ts_data.drop(i, axis=0)
    ts_data.reset_index(inplace=True, drop=True)

    # HEATMAP_DF sums up all the peak values for the heatmap output
    # heatmap_df = ts_data.groupby(['Well']).agg({'peak_presence':np.sum}).reset_index()
    heatmap_df = ts_data.groupby(["Well"])["peak_presence"].apply(calc_sum).reset_index()

    # creates a 2D array from ts_data df
    for index in range(0, len(heatmap_df.Well)):
        individual_row_char = str(heatmap_df.Well[index][0])
        individual_col_int = int(heatmap_df.Well[index][1:])
        hdf.loc[individual_row_char, individual_col_int] = float(heatmap_df.peak_presence[index])

    hdf.reset_index(inplace=True, drop=True)

    # asserts matrix ranges
    assert len(hdf) != 0
    assert len(hdf.count(axis='rows')) * len(hdf.count(axis='columns')) <= 96

    pass_or_fail_string = str(passing_criteria_check(hdf))

    # generate heatmap figure
    ax = plt.subplots(figsize=(12, 8))
    ax = sns.heatmap(hdf, annot=True, fmt=".2f", annot_kws={"size": 10}, cmap='Reds', vmin=0, vmax=20,
                     yticklabels=ROW_LABELS)  # heatmap uses hdf 2D dataframe to populate
    ax.set_ylim(top=8, bottom=.01)
    ax.invert_yaxis()  # makes sure well A1 is top left corner

    if pass_or_fail_string == 'Pass':
        plt.title('G360 EIO NTC Tool - RESULT - ' \
                  'QC Lot Release Data \n\n Job ID: {0}\nPeak Presence Summary:\n'.format(JOBID))
        plt.text(s='\
            Job ID: ' + str(JOBID) + '\n\
            URL to the Output File:' + url + OUTPUT_FILE_PATH \
                   + '{0}_G360EIOiQC_NTC_RESULT.csv'.format(JOBID) + '\n\
            Script Name: g360_eioqc_ntc.py \n\
            Script Version: v1.0' + '\n\
            Base Pair Range: ' + str(BP_MIN) + ' - ' + str(BP_MAX) + '\n\
            Lot Number: ' + str(LOT_NUMBER) + '\n\
            Part Number: ' + str(GH_PART_NUMBER) + '\n\
            Output Units: Peak Presence Value [nmol/L]' + '\n\
            Overall QC Status [Pass/Fail]: PASS', ha='left', va='bottom', fontdict={'fontsize': 9}, y=10.7, x=0)
    else:
        plt.title('G360 EIO NTC Tool - RESULT - ' \
                  'QC Lot Release Data \n\n Job ID: {0}\nPeak Presence Summary:\n'.format(JOBID))
        plt.text(s='\
            Job ID: ' + str(JOBID) + '\n\
            URL to the Output File:' + url + OUTPUT_FILE_PATH \
                   + '{0}_G360EIOiQC_NTC_RESULT.csv'.format(JOBID) + '\n\
            Script Name: g360_eioqc_ntc.py \n\
            Script Version: v1.0 ' + '\n\
            Base Pair Range: ' + str(BP_MIN) + ' - ' + str(BP_MAX) + '\n\
            Lot Number: ' + str(LOT_NUMBER) + '\n\
            Part Number: ' + str(GH_PART_NUMBER) + '\n\
            Output Units: Peak Presence Value [nmol/L]' + '\n\
            Overall QC Status [Pass/Fail]: FAIL', ha='left', va='bottom', fontdict={'fontsize': 9}, y=10.7, x=0)

    plt.savefig(OUTPUT_FILE_PATH + '{0}_G360EIOiQC_NTC_RESULT.pdf'.format(JOBID), bbox_inches='tight',
                pad_inches=0.75)  # saves figure
    # print('URL to the Output PDF: https://bifs.ghdna.io'\
    #     + OUTPUT_FILE_PATH +'{0}_G360EIOiQC_NTC_RESULT.pdf'.format(JOBID))

    output_df = pd.DataFrame(index=range(0, len(ts_data)), columns=['Well', 'Peak Presence', 'Peak Molarity [nM]'])
    for peak in range(0, len(ts_data)):
        if ts_data.peak_presence[peak] != 0.00:
            output_df['Well'][peak] = ts_data.Well[peak]
            output_df['Peak Presence'][peak] = 'Fail'
            output_df['Peak Molarity [nM]'][peak] = "{:.2f}".format(ts_data.peak_presence[peak])
        else:
            output_df['Well'][peak] = ts_data.Well[peak]
            output_df['Peak Presence'][peak] = 'Pass'
            output_df['Peak Molarity [nM]'][peak] = 0.00
    print_results(output_df, ['Well', 'Peak Presence Status', 'Peak Presence Value [nmol/L]'], pass_or_fail_string,
                  str(input_file_name_w_ext), url, OUTPUT_FILE_PATH, JOBID, \
                  LOT_NUMBER, GH_PART_NUMBER, INPUT_FILE_NAME, BP_MIN, BP_MAX)


if __name__ == "__main__":
    main()