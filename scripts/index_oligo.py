script_name = 'index_oligo_qc'
version = 'v1.0'
author = 'Tracy Nguyen'

import os
import pandas as pd
from multiprocessing import Pool
import subprocess
from subprocess import Popen
import argparse
import numpy as np
from functools import reduce
import json
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="FCID", required=True)
    parser.add_argument("--output", help="Directory for writing output (will also contain fastq files)", required=True)
    parser.add_argument("--samples", help="CSV formatted sample sheet containing indexes and oligo sequences",
                        required=True)
    parser.add_argument("--config",
                        help="CSV formatted file containing acceptance criteria cutoffs for each assay (ie contamination cutoff, cluster density cutoff, etc)",
                        required=True)
    parser.add_argument("--threads", default=6, type=int, help="Number of processors to use for analysis")
    parser.add_argument("--no-demux",
                        help="For testing purposes can skip demux (set to True), input folder doesn't need to exist",
                        action='store_true')
    parser.add_argument("--no-fastq", help="For testing purposes can skip using fastq files (set to True)",
                        action='store_true')
    parser.add_argument("--paths", help="List of directories to search for the flowcell data", required=False, type=str,
                        default="/ghds/ivd/raw, /ghds/raw")
    parser.add_argument('--folder_id', type=str, required=False)
    parser.add_argument('--URL', type=str, required=False, default='https://bifs.ghdna.io')
    parser.add_argument('--LOT_NUMBER', type=str, required=True)

    args = parser.parse_args()
    url = args.URL

    name = args.input.split('/')[-1].split('_')[-1]
    if len(name) == 0:
        name = args.input.split('/')[-2].split('_')[-1]

    if not args.output.endswith('/'):
        args.output = args.output + '/'

    args.output = args.output + args.folder_id + '/' + args.input + '_' + args.LOT_NUMBER + '/'

    try:
        os.listdir(args.output)
    except FileNotFoundError:
        os.mkdir(args.output)

    paths = list(args.paths.split(','))

    # reading in list of oligos used for this assay
    oligo_file, assay = read_oligo(args.samples, args.output, name)
    # read in configuration file with the cutoffs and correction factor
    contam_cutoff, perf_cutoff, cd_min_cutoff, cd_max_cutoff, pf_cutoff, slope, intercept = read_config(args.config)
    # find the flowcell's path
    fcid_path = find_flowcell(args.input, args.output, assay, paths)

    if fcid_path == 0:
        print("FCID " + args.input + " could not be found")
        return None
    if fcid_path == "Incorrect":
        print("Incorrect flowcell run configuration")
        return None

    seq_qual_status, report_path = initial_report(args.output, oligo_file, args.input, fcid_path, cd_min_cutoff,
                                                  cd_max_cutoff, pf_cutoff, url)
    sample_sheet, project_name, mode = generate_sample_sheet(args.output, fcid_path, oligo_file, name)

    ######## SHORTCUTS: ########
    # 1) START AFTER FASTQ FILES ARE MADE --no-demux
    if not args.no_demux:
        code = call_bcl2fastq(sample_sheet, fcid_path, args.output, threads=args.threads)
        if code != 0:
            report_file = args.output + project_name + '_SEQ_RESULT.csv'
            f = open(report_file)
            for line in f:
                print(line)
            f.close()
            print("bcl2fastq exited with non-zero exit code. Check traceback for error.")
            return None
    # 2) START AFTER READ_COUNTS.CSV FILE IS MADE --no-fastq
    if not args.no_fastq:
        read_counts_df = count_matches(args.output, oligo_file, mode, threads=args.threads)
    else:
        read_counts_df = pd.read_csv(args.output + project_name + "_read_counts.csv")
    ############################

    cont_df = calculate_contamination(args.output, read_counts_df, project_name, assay, slope, intercept)
    performance_df = calculate_performance(args.output, project_name, read_counts_df)

    seq_qual_status, contam_fail, overall_status = final_report(assay, report_path, cont_df, performance_df,
                                                                seq_qual_status, contam_cutoff, perf_cutoff)


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


def read_config(config_file):
    config_df = pd.read_csv(config_file, dtype={'Value': object})
    config_df = config_df.set_index(['Metric'])
    perf_cutoff = config_df.loc['perf_cutoff']
    cd_min_cutoff = config_df.loc['cd_min_cutoff']
    cd_max_cutoff = config_df.loc['cd_max_cutoff']
    pf_cutoff = config_df.loc['pf_cutoff']
    slope = config_df.loc['slope']['Value']
    intercept = config_df.loc['intercept']['Value']
    assay = config_df.loc['assay']['Value']

    if assay != 'G':
        contam_cutoff = config_df.loc[['contam_cutoff', 'index_cutoff']]
    else:
        contam_cutoff = config_df.loc[['contam_cutoff']]
    return contam_cutoff, perf_cutoff, cd_min_cutoff, cd_max_cutoff, pf_cutoff, slope, intercept


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


def initial_report(directory, oligo_file, fcid, fcid_path, cd_min_cutoff, cd_max_cutoff, pf_cutoff, url):
    working_directory = os.getcwd()
    if directory.startswith(working_directory):
        output_directory = directory
    else:
        output_directory = os.path.join(working_directory, directory)

    output_directory = url + output_directory

    ########## Get acceptance criteria ############
    cd_min_cutoff_val = str(cd_min_cutoff['Value'])
    cd_min_cutoff_op = str(cd_min_cutoff['Operator'])
    cd_max_cutoff_val = str(cd_max_cutoff['Value'])
    cd_max_cutoff_op = str(cd_max_cutoff['Operator'])
    pf_cutoff_val = str(pf_cutoff['Value'])
    pf_cutoff_op = str(pf_cutoff['Operator'])
    ###### Finish getting acceptance criteria #####

    config = pd.read_csv(oligo_file)
    with open(directory + config['Sample_Project'].iloc[0] + '_SEQ_RESULT.csv', 'w') as fout:
        fout.write('Purpose:,' + ' '.join(config['Sample_Project'].iloc[0].split('_')[1:]) + ' Validation\n')
        fout.write('URL Output Files:,' + output_directory + '\n')
        fout.write('Flowcell ID:,' + fcid + '\n')
        fout.write('Script Name:,' + script_name + '\n')
        fout.write('Script Version:,' + version + '\n')
        fout.write('Location of input files:,' + fcid_path + '\n')
        fout.write(
            'Location of report:,' + output_directory + config['Sample_Project'].iloc[0] + '_SEQ_RESULT.csv' + '\n')
        fout.write('Date of analysis:,' + str(datetime.today().strftime("%m/%d/%Y")) + "\n")

        ### Sequencing Quality Metrics Section
        cluster_density, cluster_pf = seq_metrics(fcid_path)
        fout.write('\nSequencing Quality Metrics:\n')
        fout.write('Metrics,Passing Criteria,Value      ,Pass/Fail\n')

        # This evaluates if cluster density is >= 120 and then if cluster density is <= 220
        if eval(str(cluster_density) + " " + cd_min_cutoff_op + " " + cd_min_cutoff_val):
            # if eval((str(cluster_density)  + " " + cd_max_cutoff_op + " " + cd_max_cutoff_val )):
            cluster_density_status = 'Pass'
            # else:
            # cluster_density_status = 'Fail'
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

    return (seq_qual_status, directory + config['Sample_Project'].iloc[0] + '_SEQ_RESULT.csv')


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
    file_out = out_dir + "bcl2fastq.log"

    '''Calls bcl2fastq with no mismatches and enforcing a read length of 20 (default)'''
    args = "bcl2fastq -R {0} -o {1} --sample-sheet {2} -p {3} --barcode-mismatches {4} --ignore-missing-bcls".format(
        in_dir, out_dir, sample_csv, threads, mismatches)
    # --with-failed-reads --ignore-missing-filter --mask-short-adapter-reads {5} --no-lane-splitting --interop-dir {6}
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
    '''Count exact matches to each oligo in Undertermined fastq.gz file. These are the reads that didn't demux bc not perfect macth to both i7 and i5 index.
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


def calculate_contamination(directory, read_counts_df, project_name, assay, slope, intercept):
    m = float(slope)
    b = float(intercept)
    df = read_counts_df.copy()
    # Get the sum of reads for each oligo = OLIGO_SUM column
    df['OLIGO_SUM'] = df.sum(axis=1, skipna=True)

    # Split the Sample_ID column
    df['TYPE'] = df['Sample_ID'].str.split('_').str[0]
    df['SET'] = df['Sample_ID'].str.split('_').str[2]
    df['OLIGO'] = df['Sample_ID'].str.split('_').str[3]

    # Melt dataframe and get percentage of reads
    df = df.drop(columns=["Sample_ID"])
    df = pd.melt(df, id_vars=["TYPE", "SET", "OLIGO", "OLIGO_SUM"])
    df['i7'] = df['variable'].str.split(':').str[0].str.split('_').str[1]
    df['perc_reads'] = 100 * df['value'] / df['OLIGO_SUM']

    # Only keep testing oligos (TYPE == T)
    df = df[df['TYPE'] == "T"]

    # Remove if both i7 and i5 are equal to the oligo
    if assay != 'G':
        df['i5'] = df['variable'].str.split(':').str[1].str.split('_').str[1]
        df.loc[df['i5'].isnull(), 'i5'] = df['i7']
        df = df[(df['OLIGO'] != df['i7']) | (df['OLIGO'] != df['i5'])]
        df = df[['TYPE', 'SET', 'OLIGO', 'i7', 'i5', 'value', 'OLIGO_SUM', 'perc_reads']]
    else:
        df = df[(df['OLIGO'] != df['i7'])]
        df = df[['TYPE', 'SET', 'OLIGO', 'i7', 'value', 'OLIGO_SUM', 'perc_reads']]

    df.perc_reads = df.perc_reads.round(2)

    # do not apply correction factor if the perc_reads is 0
    df['predicted_perc_contam'] = np.where(df['perc_reads'] > 0, df['perc_reads'].subtract(b).divide(m),
                                           df['perc_reads'])

    # round the percentage of predicted reads to two decimal places
    df.predicted_perc_contam = df.predicted_perc_contam.round(2)

    # write out contamination csv file
    df.rename(columns={'value': 'oligo_reads', 'OLIGO_SUM': 'total_oligo_reads'}, inplace=True)
    df.to_csv(directory + project_name + '_contamination.csv', index=False)

    # Use only the testing oligos
    df = df.drop(columns=["oligo_reads", "total_oligo_reads", "perc_reads", "TYPE"])

    if assay == 'G':
        # flag known index pairs that will have eio bleeding for G360:
        df = df.pivot_table(index=['i7', 'OLIGO'], columns='SET', values='predicted_perc_contam', fill_value=0)
        df['Median'] = df.median(axis=1)
        df = df.reset_index(drop=False)
        df['eio_bleeding_pair'] = np.where((((df['i7'] == '26') & (df['OLIGO'] == '14')) | (
                    (df['i7'] == '09') & (df['OLIGO'] == "41")) | ((df['i7'] == '16') & (df['OLIGO'] == "23"))), 'TRUE',
                                           'FALSE')
        eio_bleeding_df = df[df['eio_bleeding_pair'] == "TRUE"]
        df = df[df['eio_bleeding_pair'] != "TRUE"]
        df = df.append(eio_bleeding_df)
        df = df.set_index(['i7', 'OLIGO'])
    else:
        if assay != 'O':
            df = df.pivot_table(index=['i7', 'i5', 'OLIGO'], columns='SET', values='predicted_perc_contam',
                                fill_value=0)
            df['Median'] = df.median(axis=1)
        else:
            # REP will be either A, B, C, D
            df['REP'] = df['SET'].apply(lambda x: x[-1])
            # SET will be either Set1, Set2, Set3
            df['SET'] = df['SET'].apply(lambda x: x[:-1])
            # Pivot table
            df = df.pivot_table(index=['i7', 'i5', 'OLIGO', 'REP'], columns='SET', values='predicted_perc_contam',
                                fill_value=0)
            df['Median'] = df.median(axis=1)

    df = df.reset_index(drop=False)
    df = df.rename_axis(None, axis=1)
    df.index.names = ['ix']
    df.to_csv(directory + project_name + '_contamination_set.csv', index=True, index_label='ix')

    return df


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


def final_report(assay, report_path, cont_df, performance_df, seq_qual_status, contam_cutoff, perf_cutoff):
    ######## TYPE OF INDEX OLIGO  #########
    if assay == 'O':
        io_type = 'LIO'
    else:
        io_type = 'EIO'

    ######## TYPE OF CONTAMINATION TO REPORT OUT #########
    contam_mode = len(contam_cutoff['Value'])

    # if mode != 1, then also evaluate at i7-level and i5-level contamination.
    if contam_mode != 1:
        index_cutoff = contam_cutoff.loc['index_cutoff']
        index_cutoff_val = str(index_cutoff['Value'])
        index_cutoff_op = str(index_cutoff['Operator'])

    ######### GET ACCEPTANCE CRITERIA ##########
    contam_cutoff = contam_cutoff.loc['contam_cutoff']
    contam_cutoff_val = str(contam_cutoff['Value'])
    contam_cutoff_op = str(contam_cutoff['Operator'])
    perf_cutoff_val = str(perf_cutoff['Value'])
    perf_cutoff_op = str(perf_cutoff['Operator'])
    #### FINISH GETTING ACCEPTANCE CRITERIA ####

    if assay == 'G':
        # Separate known eio bleeding pairs from the other pairs
        eio_bleed = cont_df[cont_df['eio_bleeding_pair'] == "TRUE"]
        cont_df = cont_df[cont_df['eio_bleeding_pair'] != "TRUE"]

    cont_contam_fail = 0
    cont_sample_fail = 0
    failed_sample_list = list()
    failed_i7_list = list()

    # Writing out report
    with open(report_path, 'a') as fout:
        ##################################################################################################
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

            # i5-level contamination:
            fout.write('\nCross contamination quality metrics (' + io_type + ' Contamination: i5 Contamination):\n')
            df_oligo_i5 = cont_df[cont_df['i5'] != cont_df['OLIGO']]
            df_oligo_i5 = df_oligo_i5.set_index('OLIGO')
            if assay == 'O':
                rep = cont_df['REP'].unique()
            else:
                rep = ['A']
            for x in rep:
                if len(rep) == 1:
                    df_oligo_i5_rep = df_oligo_i5
                else:
                    df_oligo_i5_rep = df_oligo_i5[df_oligo_i5['REP'] == x]
                    fout.write(io_type + ' Replicate ' + x + ' i5 Contamination\n')
                fout.write('Index,Passing Criteria,Value,Pass/Fail,Source of Contamination\n')
                for index, new_df in df_oligo_i5_rep.groupby(level=0):
                    max_row = new_df[new_df['Median'] == new_df['Median'].max()]
                    oligo = io_type + '_' + str(index)
                    i5 = io_type + '_' + str((max_row.iloc[0]['i5']))
                    median = max_row.iloc[0]['Median']
                    if eval(str(median) + " " + index_cutoff_op + " " + index_cutoff_val):
                        fout.write(
                            oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Pass\n")
                    else:
                        fout.write(
                            oligo + "," + index_cutoff_op + index_cutoff_val + "%," + "%0.2f" % median + "%,Fail," + i5 + "\n")
                        cont_contam_fail += 1
        '''
        ################################################################################################################################
        ### Concentration Metrics Section
        fout.write('\nPerformance quality metrics:\n')
        fout.write('Performance quality metrics (Concentration):\n')
        if len(performance_df) != 0:
            fout.write('Index,Passing Criteria,Value,Pass/Fail\n')
            performance_fail = 0
            for index, row in performance_df.iterrows():
                i7  = io_type + '_' + row['Index']
                median = row['Median']
                if eval(str(median) + " " + perf_cutoff_op + " " + perf_cutoff_val) : 
                    fout.write(i7 + "," + perf_cutoff_op +  perf_cutoff_val + '%,' + "%0.2f" %row['Median'] + "%,Pass\n")
                else:
                    fout.write(i7 + "," + perf_cutoff_op +  perf_cutoff_val + '%,' + "%0.2f" %row['Median'] + "%,Fail\n")
                    performance_fail +=1           
        else:
            performance_fail = 'NA'
            fout.write('Reference oligos not detected in run\n')
        '''
        ################################################################################################################################################
        ### Summary Section
        fout.write('\nCross-Contamination Qualification Summary:')
        fout.write('\nQC Metric,Passing Criteria,Comments,Status\n')

        ## Overall Status: Pass/Fail/Repeat
        if seq_qual_status != "Pass":
            overall_status = "FAIL"
        elif cont_contam_fail != 0:
            overall_status = "FAIL"
        else:
            overall_status = "PASS"

        if seq_qual_status != "Pass":
            seq_qual_status = "FAIL"
        else:
            seq_qual_status = 'PASS'

        if overall_status != "PASS" and seq_qual_status == "FAIL":
            fout.write('Overall Cross-Contamination QC,All sections pass,Sequence Qualification Fail,' + str(
                overall_status) + '\n')
        elif overall_status != "PASS" and cont_contam_fail != 0:
            fout.write('Overall Cross-Contamination QC,All sections pass,Cross-Contamination Qualification Fail,' + str(
                overall_status) + '\n')
        else:
            fout.write('Overall Cross-Contamination QC,All sections pass,NA,' + str(overall_status) + '\n')
        ################################################################################################################################################
        # fout.write('Sequencing Qualification Status:,Both Cluster Density and Cluster Passing Filter pass,,' + seq_qual_status + '\n')
        fout.write('Sequencing Qualification,Both Cluster Density and Cluster Passing Filter pass,NA,' + str(
            seq_qual_status) + '\n')
        if cont_contam_fail == 0:
            fout.write('Cross Contamination Qualification,All ' + io_type + 's pass,NA,PASS\n')
        else:
            fout.write('Cross Contamination Qualification,All ' + io_type + 's pass,' + io_type + 's ' + str(
                failed_i7_list).replace(',', ' ').replace(']', '').replace('[', '') + ' failed,FAIL\n')
        '''
        if performance_fail == 0:
            fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,Pass\n')
        else:
            if performance_fail == 'NA':
                fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,NA\n')
            else:
                fout.write('\nPerformance Qualification Status:,All ' + io_type + 's pass,,Fail\n')

        '''
        fout.write('\nPerformed By:\n')
        fout.write('Performed Date:\n')

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


if __name__ == '__main__':
    main()