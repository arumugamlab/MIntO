#!/usr/bin/python

'''
Merge gene abundances raw profiles and exclude not present and not expressed genes

Authors: Carmen Saenz
'''

import sys
import pandas as pd
import csv

threads_n = sys.argv[1]
memory_lim = sys.argv[2]
wd = sys.argv[3]
omics = sys.argv[4] #'metaG_metaT'
map_reference = sys.argv[5]
normalization = sys.argv[6]
identity = sys.argv[7]

#soft, hard = 10**7, 10**7

tpm_profile_all = pd.read_csv(wd + '/output/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p' + identity + '.' + normalization + '.csv',sep='\t')

if omics == 'metaG_metaT':
    metaG_file_raw = pd.read_csv(wd + '/metaG/9-mapping-profiles/' + map_reference + '/genes_abundances.p' + identity + '.bed',sep='\t')
    metaT_file_raw = pd.read_csv(wd + '/metaT/9-mapping-profiles/' + map_reference + '/genes_abundances.p' + identity + '.bed',sep='\t')
    #################################
    tpm_profile_all_sub = tpm_profile_all[['chr','start','stop','name','score','strand','source','feature','frame','info']]
    metaG_new_names = [(i,'metaG.'+i) for i in metaG_file_raw.iloc[:, 10:].columns.values]
    metaG_file_raw.rename(columns = dict(metaG_new_names), inplace=True)
    metaG_file_raw_sub = pd.merge(metaG_file_raw, tpm_profile_all_sub, on=['chr','start','stop','name','score','strand','source','feature','frame','info'], how = 'right')
    metaT_new_names = [(i,'metaT.'+i) for i in metaT_file_raw.iloc[:, 10:].columns.values]
    metaT_file_raw.rename(columns = dict(metaT_new_names), inplace=True)
    tpm_profile_raw_all = pd.merge(metaG_file_raw_sub, metaT_file_raw, on=['chr','start','stop','name','score','strand','source','feature','frame','info'], how = 'left')

    tpm_profile_raw_all.to_csv(wd + '/output/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p' + identity + '.bed', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)
else:
    omics_file_raw = pd.read_csv(wd + '/' + omics +'/9-mapping-profiles/' + map_reference + '/genes_abundances.p' + identity + '.bed',sep='\t')
    #################################
    tpm_profile_all_sub = tpm_profile_all[['chr','start','stop','name','score','strand','source','feature','frame','info']]
    omics_new_names = [(i, omics+'.'+i) for i in omics_file_raw.iloc[:, 10:].columns.values]
    omics_file_raw.rename(columns = dict(omics_new_names), inplace=True)
    omics_file_raw_sub = pd.merge(omics_file_raw, tpm_profile_all_sub, on=['chr','start','stop','name','score','strand','source','feature','frame','info'], how = 'right')
    omics_file_raw_sub.to_csv(wd + '/output/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p' + identity + '.bed', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)
