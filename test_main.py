# import modules
import os
import sys
import numpy as np
import pandas as pd
import string
import itertools
import configparser

repo_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split('/'))
script_path = os.path.join(repo_path, "scripts")
sys.path.append(script_path)
import utility as test

# Variables
config_file = 'configfile.ini'
config = test.read_config(config_file)
# pcr = config['pcr']
# print(pcr['bp_min'])
# print(pcr['bp_max'])

# mode = 'seq'
# job_id = 12345

# ts_parse test
# if mode == 'ntc':
#     input_file = "/Users/kkwock/Documents/GitHub/kkwock/input/20220407_Post-6806-NTC_kk10 - 2022-04-07 - 09-06-01-D1000_compactPeakTable.csv"
# elif mode == 'pcr':
#     input_file = "/Users/kkwock/Documents/GitHub/kkwock/input/20220407_Post-6806-REF_kk10 - 2022-04-07 - 09-05-47-D1000_compactPeakTable.csv"
# elif mode =='seq':
#     input_file = '/Users/kkwock/Documents/GitHub/kkwock/input/HJ5WFBGXL_G360_EIO_QC_read_counts.csv'

# print(input_file)
# filtered_ts, hdf, pass_fail, test_iqr = test.ts_parse(input_file, config, mode)

# read_counts_df = pd.read_csv(input_file)

# df, pass_fail, list = test.contam_pass_fail('./', read_counts_df, job_id, config, False)

# test.results_config(repo_path, 'seq', 'PASS',
#                     80, 80, 'PASS', ['A','B','C'],
#                     'NA', 'NA', 'NA', 'NA')
#
# print(repo_path)
# # pdf = test.calculate_performance('./', str(job_id), read_counts_df)
#
# results_path = f'{repo_path}/results.ini'
# os.path.isfile(results_path)

# df test.calculate_contamination(directory, read_counts_df, job_id, config)

# if(mode in config):
#     print("TRUE")
# else:
#     print("FALSE")
#
# print(config['results']['config_list'])

df = pd.read_csv('/Users/kkwock/Documents/GitHub/kkwock/input/HJ5WFBGXL_G360_EIO_QC_read_counts.csv')

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
df['Perc_Reads'] = pd.eval(round((df['Count'] / df['Total_Reads'])*100, 2))

directory = repo_path + '/output/'
job_id = 'ABCD_1234'
seq = config['seq']
readfrac_op = str(seq['readfrac_operator'])
readfrac_min = seq['readfrac_min']
df.to_csv(directory + str(job_id) + '_read-fractions.csv', index=True, index_label='ix')

def op(config_op):
    import operator
    op_dict = {
        ">": lambda a,b: operator.gt(a,b),
        "<": lambda a,b:operator.lt(a,b),
        ">=": lambda a,b:operator.ge(a,b),
        "<=": lambda a,b:operator.le(a,b),
        "!=": lambda a,b:operator.ne(a,b),
        "==": lambda a,b:operator.eq(a,b)
    }

    return op_dict[config_op]

rf_op = op(readfrac_op)
seq_df = pd.DataFrame({"Lot": '1234',
                       "Metric": ["Cluster Density", "Cluster Passing Filter", "Cross-Contamination", "Read Fraction"],
                       "Outcome": "NA",
                       "Comment": ""})
seq_df.loc[0]['Outcome'] = "Fail" if not rf_op(df['Perc_Reads'][1], float(readfrac_min)) else "Pass"

print(list(df.loc[np.invert(rf_op(df['Perc_Reads'], float(readfrac_min)))]['Sample_ID']))

# print(list(df.loc[[not x for x in rf_op(df['Perc_Reads'], float(readfrac_min))]]['Sample_ID']))
# pd.eval(df.loc[df['Read_Frac'], readfrac_op, readfrac_min])
# print(list)