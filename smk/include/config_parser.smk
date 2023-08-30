#!/usr/bin/env python

'''
Config file parser

Authors: Carmen Saenz, Mani Arumugam
'''

################################
# config file name:
# Getting the configfile's name is simple but just not obvious.
# See: https://stackoverflow.com/a/72096363
################################
config_path = 'configuration yaml file' #args[args_idx+1]
if '--configfile' in sys.argv:
    i = sys.argv.index('--configfile')
    config_path = sys.argv[i+1]
elif '--configfiles' in sys.argv:
    i = sys.argv.index('--configfiles')
    config_path = sys.argv[i+1]

print(" *******************************")
print(" Reading configuration yaml file: ", config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file
if config['PROJECT'] is None:
    print('ERROR in ', config_path, ': PROJECT variable is empty. Please, complete ', config_path)
else:
    project_id = config['PROJECT']

if config['working_dir'] is None:
    print('ERROR in ', config_path, ': working_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['working_dir']) is False:
    print('WARNING in ', config_path, ': working_dir path does not exit. The directory will be created by MIntO')
    working_dir = config['working_dir']
else:
    working_dir = config['working_dir']

if config['omics'] in ('metaG', 'metaT', 'metaG_metaT'):
    omics = config['omics']
else:
    print('ERROR in ', config_path, ': omics variable is not correct. "omics" variable should be metaG, metaT or metaG_metaT.')

if config['minto_dir'] is None:
    print('ERROR in ', config_path, ': minto_dir variable in configuration yaml file is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]
    script_dir=config["minto_dir"]+"/scripts"

if config['METADATA'] is None:
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata=config["METADATA"]
    #metadata="None"
elif config['METADATA'] == "None":
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata=config["METADATA"]
if config['METADATA'] is not None:
    if path.exists(config['METADATA']) is False:
        print('ERROR in ', config_path, ': METADATA variable path does not exit. Please, complete ', config_path)
    else:
        metadata=config["METADATA"]

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
else:
    print('WARNING in ', config_path, ': ILLUMINA list of samples is missing. Proceed with caution.')
