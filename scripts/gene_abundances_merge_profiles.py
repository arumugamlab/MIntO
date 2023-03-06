#!/usr/bin/python

'''
Merge gene abundances profiles and exclude not present and not expressed genes

Authors: Carmen Saenz, Mani Arumugam
'''

import sys
import pandas as pd
import csv

map_reference = sys.argv[1]
out_csv = sys.argv[2]
metaG_csv = sys.argv[3]
metaT_csv = sys.argv[4]

#soft, hard = 10**7, 10**7

metaG_file = pd.read_csv(metaG_csv, sep='\t')
metaT_file = pd.read_csv(metaT_csv, sep='\t')

if map_reference == 'db_genes':
    metaG_file.rename(columns = {list(metaG_file)[0]:'Genes'}, inplace=True)
    metaT_file.rename(columns = {list(metaT_file)[0]:'Genes'}, inplace=True)
    metaG_file_sub = metaG_file.loc[~(metaG_file.iloc[:, 1:]==0).all(axis=1)]
    metaT_file_sub = metaT_file.loc[~(metaT_file.iloc[:, 1:]==0).all(axis=1)]
    tpm_profile_all = pd.merge(metaG_file_sub, metaT_file_sub, on='Genes', how = 'inner')
    #tpm_profile_all_last=tpm_profile_all.fillna(0)
    #tpm_profile_all = pd.merge(metaG_file, metaT_file, on='Genes', how = 'inner')
    #tpm_profile_all_sub = tpm_profile_all.loc[~(tpm_profile_all.iloc[:, 1:]==0).all(axis=1)]

elif map_reference == 'MAGs_genes' or map_reference == 'reference_genes':
    metaG_file_sub = metaG_file.loc[~(metaG_file.iloc[:, 12:]==0).all(axis=1)]
    metaT_file_sub = metaT_file.loc[~(metaT_file.iloc[:, 12:]==0).all(axis=1)]
    tpm_profile_all = pd.merge(metaG_file_sub, metaT_file_sub, on=['coord','chr','start','stop','name','score','strand','source','feature','frame','info','gene_length'], how = 'outer')
    #tpm_profile_all_last=tpm_profile_all.fillna(0)
    #tpm_profile_all = pd.merge(metaG_file, metaT_file, on=['coord','chr','start','stop','name','score','strand','source','feature','frame','info','gene_length'], how = 'outer')
    #tpm_profile_all_sub = tpm_profile_all.loc[~(tpm_profile_all.iloc[:, 12:]==0).all(axis=1)]

#tpm_profile_all = metaT_file.append(metaT_file, ignore_index=True, sort=False)

tpm_profile_all.to_csv(out_csv, index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)
#tpm_profile_all_sub.to_csv(wd + '/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p95.TPM.csv', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)

