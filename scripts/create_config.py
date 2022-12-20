# create config.INI file
import os
import configparser
from pathlib import Path
config = configparser.ConfigParser()

# Create Structure

## Cross Contamination
config.add_section('cross_contamination')
config.set('cross_contamination', 'contam_cutoff', '0.10')
config.set('cross_contamination', 'contam_operator', '<=')
config.set('cross_contamination', 'perf_cutoff', '40.00')
config.set('cross_contamination', 'perf_operator', '>=')
config.set('cross_contamination', 'cd_min', '120')
config.set('cross_contamination', 'cd_min_operator', '>=')
config.set('cross_contamination', 'cd_max', '220')
config.set('cross_contamination', 'cd_max_operator', '<=')
config.set('cross_contamination', 'pf_cutoff', '80')
config.set('cross_contamination', 'pf_operator', '>=')
config.set('cross_contamination', 'slope', '1')
config.set('cross_contamination', 'intercept', '0')
config.set('cross_contamination', 'assay', 'G')

## PCR
config.add_section('pcr')
config.set('pcr', 'iQR_min', '60')
config.set('pcr', 'iQR_operator', '<=')
config.set('pcr', 'bp_min', '130')
config.set('pcr', 'bp_max', '165')

## NTC
config.add_section('ntc')
config.set('ntc', 'conc', '0.1')
config.set('ntc', 'conc_operator', '>=')
config.set('ntc', 'bp_min', '100')
config.set('ntc', 'bp_max', '500')

# Write the new structure to the new file
file = Path(os.path.realpath(__file__)).parent
path = Path(file).parent
config_dir = os.path.join(path, 'configs')


with open(f"{config_dir}/configfile.ini", 'w') as configfile:
    config.write(configfile)