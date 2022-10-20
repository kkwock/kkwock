# G360 EIO PCR Prod Tool SG3 02/24/2021
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

# user inputs TS file name including .csv extension


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


def pcrprod_passing_criteria_check(heatmap_df):
    output_df = pd.DataFrame()
    test_df = pd.DataFrame()
    ref_df = pd.DataFrame()

    for i in range(0, len(heatmap_df)):
        if int(heatmap_df.Well[i][1:]) in range(1, 7):
            test_df = test_df.append(heatmap_df.loc[i])
        if int(heatmap_df.Well[i][1:]) in range(7, 9):
            ref_df = ref_df.append(heatmap_df.loc[i])

    # to make the wells match up between reps (A6 becomes A1 because they share the same EIO)
    ref_df.reset_index(inplace=True)
    for i in range(0, len(ref_df)):
        ref_df.Well[i] = str(ref_df.Well[i][:1]) + str(int(ref_df.Well[i][1:]) - 6)

    test_qc_df = test_df.peak_presence.quantile([.25, .5, .75]).round(4)
    test_qr_value = float(test_qc_df[0.25])

    output_df = test_df.merge(ref_df, on='Well', sort=True, how='outer')

    output_df = output_df.rename(
        columns={'Well': 'EIO_Well', 'peak_presence_x': 'peak_presence_test', 'peak_presence_y': 'peak_presence_ref'})

    # orders wells for EIO Naming in the next code chunk to work
    orderedwelllist = sorted(list(output_df.EIO_Well), key=lambda x: int(x[1]))
    output_df['EIO_Well'] = pd.Categorical(output_df['EIO_Well'], categories=orderedwelllist)
    output_df = output_df.sort_values('EIO_Well')

    output_df[['EIO_Name']] = ''
    for i in range(0, 66):
        output_df['EIO_Name'][i] = 'EIO' + str(i + 1)
    for i in range(48, 66):
        output_df['EIO_Name'][i] = 'EIO' + str(i - 47)

    output_df['test_well_qc_check'] = np.where(output_df["peak_presence_test"] >= 90, "pass", "fail")
    output_df['ref_well_qc_check'] = np.where(output_df["peak_presence_ref"] >= 90, "pass", "fail")
    output_df.reset_index(inplace=True)
    output_df['ref_well_qc_check'][17:49] = 'EIO_Ignored'
    output_df['peak_presence_ref'][17:49] = 'EIO_Ignored'

    output_df = output_df[
        ['EIO_Well', 'peak_presence_test', 'test_well_qc_check', 'peak_presence_ref', 'ref_well_qc_check']]

    return (output_df, test_qr_value)


def pcrprod_print_results(qc_df, test_qr_value, url, OUTPUT_FILE_PATH, JOBID, LOT_NUMBER, GH_PART_NUMBER,
                          PCRPROD_INPUT_FILE_NAME):
    output_file_name = (OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_RESULT.csv'.format(JOBID))

    with open(output_file_name, "w") as output_file:

        output_file.write('Purpose:,G360 v2.7 EIO IQC PCR Productivity Tool - QC Lot Release Data\n')
        output_file.write('Job ID:,{0}\n'.format(JOBID))
        output_file.write(
            'URL to Output File:,' + url + OUTPUT_FILE_PATH + '{0}_G360EIOiQC_pcrProd_RESULT.csv'.format(JOBID))
        output_file.write('Base Pair Range:,{0} - {1}\n'.format(130, 165))
        output_file.write('Lot Number:,{0}\n'.format(LOT_NUMBER))
        output_file.write('Part Number:,{0}\n'.format(GH_PART_NUMBER))
        output_file.write('Script Name:,g360_eioqc_pcrprod_std.py\n')
        output_file.write('Script Version:,v1.0\n')
        output_file.write('Input File Name:,{0}\n'.format(PCRPROD_INPUT_FILE_NAME))
        output_file.write('Location of input files:,{0}\n'.format(PCRPROD_INPUT_FILE_NAME))
        output_file.write('Date Data Analyzed:,{0}\n'.format(TODAY_STRING))

        output_file.write('\nPCR Productivity Standard Overall Qualification Summary:\n')
        output_file.write('Base Pair Range,Passing Criteria,Lower IQR Value,Status\n')

        if test_qr_value >= 100:
            output_file.write('{0} - {1},Lower IQR >= 100nM, {2},PASS\n'.format(130, 165, test_qr_value))
        else:
            output_file.write('{0} - {1},Lower IQR >= 100nM, {2},FAIL\n'.format(130, 165, test_qr_value))

        output_file.write('\nQC Release Data\n')

        qc_df.to_csv(output_file, index=False, line_terminator='\n')

        output_file.write('\nPerformed By:,{0}\n'.format(PERFORMED_BY))
        output_file.write('Performed Date:,{0}\n'.format(PERFORMED_DATE))

    output_file.close()
    # print('\n\nURL to PCR Productivity [Standard] Output CSV: https://bifs.ghdna.io' + '{0}\n'.format(output_file_name))


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


def filter_peaks(ts_data, BP_MIN, BP_MAX):
    # populates a dictionary with peak ts_data if within bp range
    # populates a dictionary with "peaks_of_interest" (poi) to later refrence against to rid of duplicate instances when
    #   a poi value and a 0.00 value are found in the same well
    peaks_found = 0
    peak_dict = pd.DataFrame()  # will contain only peaks of interest
    for peak in range(0, len(ts_data.Well)):
        if ts_data.sizebp[peak] >= BP_MIN and ts_data.sizebp[peak] <= BP_MAX:
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
    parser = argparse.ArgumentParser()
    parser.add_argument('-LOT_NUMBER', type=str, required=True)
    parser.add_argument('-PCR_PRODstd_input_file', type=str, required=True)
    parser.add_argument('-FCID', type=str, required=True)
    parser.add_argument('-GHPN', type=str, required=False, default='')
    parser.add_argument('-output_path', type=str, required=True,
                        default='/ghds/groups/reagent_qc/g360/eio_qc/output_files/')
    parser.add_argument('-folder_id', type=str, required=False)
    parser.add_argument('-URL', type=str, required=False, default='https://bifs.ghdna.io')

    args = parser.parse_args()

    LOT_NUMBER = args.LOT_NUMBER
    OUTPUT_UNITS = 'molarity_output_units'
    PCRPROD_INPUT_FILE_NAME = args.PCR_PRODstd_input_file
    JOBID = args.FCID + '_' + args.LOT_NUMBER
    GH_PART_NUMBER = args.GHPN
    url = args.URL
    BP_MIN = int(130)  # args.user_input_bp_min
    BP_MAX = int(165)

    OUTPUT_FILE_PATH = args.output_path + args.folder_id + '/' + JOBID + '/'
    if not os.path.exists(OUTPUT_FILE_PATH):
        os.mkdir(OUTPUT_FILE_PATH)

    wells_list = list()

    for j in range(0, 8):
        for i in range(1, 13):
            wells_list.append(str(ROW_LABELS[j]) + str(i))

    PCRPROD_INPUT_FILE_NAME = args.PCR_PRODstd_input_file

    csv_data = pd.DataFrame

    hdf = pd.DataFrame()  # will contain filtered tapestation ts_data of all wells
    csv_input_file = pd.DataFrame()

    csv_input_file = pd.read_csv(PCRPROD_INPUT_FILE_NAME, encoding='latin1')

    input_header_list = list(csv_input_file.columns)
    check_vital(str('compactPeakTable.csv') in str(args.PCR_PRODstd_input_file),
                'Input file name must contain \'compactPeakTable.csv\' and be a compact peak table csv file.',
                input_header_list)
    check_vital('Well' in input_header_list, 'Input File Format Incorrect - Well Column Missing.', input_header_list)
    check_vital('Size [bp]' in input_header_list, 'Input File Format Incorrect - Size [bp] Column Missing.',
                input_header_list)
    check_vital('Peak Molarity [nmol/l]' in input_header_list,
                'Input File Format Incorrect - Peak Molarity [nmol/l] Column Missing.', input_header_list)

    # save ts_data file w new header names for easier use downstream
    csv_data = csv_input_file.rename(columns={"Size [bp]": "sizebp", "Peak Molarity [nmol/l]": "peak_molarity"})

    ts_data = pd.DataFrame(
        columns={"sizebp", "peak_molarity", "peak_presence", 'Well'})  # contains input csv file information

    ts_data.Well = csv_data.Well
    ts_data.sizebp = csv_data.sizebp
    ts_data.peak_molarity = csv_data.peak_molarity

    ts_data, peaks_found, peak_dict = filter_peaks(ts_data, BP_MIN, BP_MAX)

    if peaks_found == 0:
        peak_oi_wells = list('NA')  # without this the program will error out if no POI'S are found
    else:
        peak_oi_wells = list(peak_dict.Well)  # list of wells where peaks of interest (POI) were found

    # trim ts_data frame so we only have the ts_data we want to use
    # (drops a bunch of duplicate zero instances bc peak was out of range)
    ts_data = ts_data[['Well', 'peak_presence']].drop_duplicates()
    ts_data.reset_index(inplace=True, drop=True)

    # filters out the electronic ladder rows instnaces
    for i in range(0, len(ts_data.Well)):
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

    # check_vital(len(heatmap_df) == 64, 'Incorrect number of wells/plate layout in compactPeakTable.csv (Wells A1-H8 Required)', list())
    # no pass if set has nan values in it

    # creates a 2D array from ts_data df
    for index in range(0, len(heatmap_df.Well)):
        individual_row_char = str(heatmap_df.Well[index][0])
        individual_col_int = int(heatmap_df.Well[index][1:])
        hdf.loc[individual_row_char, individual_col_int] = float(heatmap_df.peak_presence[index])

    # hdf.reset_index(inplace=True, drop=True)
    # hdf.sort(columns=list(range(0,12)))
    hdf = hdf.sort_index(axis=1, ascending=True)

    # asserts matrix ranges
    assert len(hdf) != 0

    # correct plate layout error catch
    set_not_full = 0
    for col in range(1, 9):  # all 8 cols
        for row in range(0, 8):  # all 8 rows
            try:
                if np.isnan(hdf[col][row]) == True:
                    set_not_full = set_not_full + 1
            except KeyError:
                set_not_full = set_not_full + 1
    if set_not_full != 0:
        print('Incorrect number of wells/plate layout in compactPeakTable.csv (Wells A1-H8 Required)')
        sys.exit()

    qc_df, test_qr_value = pcrprod_passing_criteria_check(heatmap_df)

    qc_df = qc_df[['EIO_Well', 'peak_presence_test', 'test_well_qc_check', 'peak_presence_ref', 'ref_well_qc_check']]

    pcrprod_print_results(qc_df, test_qr_value, url, OUTPUT_FILE_PATH, JOBID, LOT_NUMBER, GH_PART_NUMBER,
                          PCRPROD_INPUT_FILE_NAME)

    heatmap_df.reset_index(inplace=True)
    orderedwelllist_hdf = sorted(list(heatmap_df.Well), key=lambda x: int(x[1]))
    heatmap_df['Well'] = pd.Categorical(heatmap_df['Well'], categories=orderedwelllist_hdf)
    heatmap_df = heatmap_df.sort_values('Well')

    heatmap_df.reset_index(inplace=True, drop=True)
    # heatmap_df[['EIO_Name']] = ''
    heatmap_df[['type']] = ''
    for j in range(0, 48):
        # heatmap_df['EIO_Name'][j] = 'EIO' + str(j + 1)
        heatmap_df['type'][j] = 'test'
    for j in range(48, 66):
        # heatmap_df['EIO_Name'][j] = 'EIO' + str(j - 47)
        heatmap_df['type'][j] = 'ref'

    pcrprod_generate_figures(heatmap_df, test_qr_value, OUTPUT_FILE_PATH, JOBID)


if __name__ == "__main__":
    main()