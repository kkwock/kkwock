import os
import configparser
from pathlib import Path

# path variables
file = os.path.dirname(os.path.realpath(__file__)).replace("/scripts", "/configs")
print(file)
config_dir = file

# Read config.ini file
config = configparser.ConfigParser()
config.read(f"{config_dir}/configfile.ini")
contam = config['cross_contamination']
pcr = config['pcr']
ntc = config['ntc']
conc_operator = config.get('ntc', 'conc_operator')
conc = ntc['conc']

# Print out configs to test
print(pcr['iQR_min'])
print(pcr['bp_min'])
print(ntc['bp_min'])
print(pcr['bp_max'])
print(conc_operator)

