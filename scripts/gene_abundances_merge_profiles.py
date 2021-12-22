#!/usr/bin/python

'''
Merge gene abundances profiles and exclude not present and not expressed genes

Authors: Carmen Saenz
'''

#from Bio import SeqIO
#import os, sys
#from glob import glob
import sys
import pandas as pd
import csv
#import resource

threads_n = sys.argv[1]
memory_lim = sys.argv[2]
wd = sys.argv[3]
omics = sys.argv[4] #'metaG_metaT'
map_reference = sys.argv[5]
normalization = sys.argv[6]
identity = sys.argv[7]

#soft, hard = 10**7, 10**7

# wd = '/emc/cbmr/users/rzv923/ibdmdb'
# omics = 'metaG_metaT'
# map_reference = 'MAGs_genes'
# normalization = 'TPM'
# identity = '95'

#metaG_file = open(wd + '/metaG/6-mapping-profiles/BWA_reads-' + map_reference + '/genes_abundances.p' + identity + '.' + normalization + '.csv', 'r')
metaG_file = pd.read_csv(wd + '/metaG/6-mapping-profiles/BWA_reads-' + map_reference + '/genes_abundances.p' + identity + '.' + normalization + '.csv',sep='\t')
#metaT_file = open(wd + '/metaT/6-mapping-profiles/BWA_reads-' + map_reference + '/genes_abundances.p' + identity + '.' + normalization + '.csv', 'r')
metaT_file = pd.read_csv(wd + '/metaT/6-mapping-profiles/BWA_reads-' + map_reference + '/genes_abundances.p' + identity + '.' + normalization + '.csv',sep='\t')

if map_reference == 'db_genes':
    metaG_file.rename(columns = {list(metaG_file)[0]:'Genes'}, inplace=True)
    metaT_file.rename(columns = {list(metaT_file)[0]:'Genes'}, inplace=True)
    metaG_file_sub = metaG_file.loc[~(metaG_file.iloc[:, 1:]==0).all(axis=1)]
    metaT_file_sub = metaT_file.loc[~(metaT_file.iloc[:, 1:]==0).all(axis=1)]
    tpm_profile_all = pd.merge(metaG_file_sub, metaT_file_sub, on='Genes', how = 'inner')
    #tpm_profile_all = pd.merge(metaG_file, metaT_file, on='Genes', how = 'inner')
    #tpm_profile_all_sub = tpm_profile_all.loc[~(tpm_profile_all.iloc[:, 1:]==0).all(axis=1)]

elif map_reference == 'MAGs_genes' or map_reference == 'reference_genes':
    metaG_file_sub = metaG_file.loc[~(metaG_file.iloc[:, 12:]==0).all(axis=1)]
    metaT_file_sub = metaT_file.loc[~(metaT_file.iloc[:, 12:]==0).all(axis=1)]
    tpm_profile_all = pd.merge(metaG_file_sub, metaT_file_sub, on=['coord','chr','start','stop','name','score','strand','source','feature','frame','info','gene_lenght'], how = 'outer')

    #tpm_profile_all = pd.merge(metaG_file, metaT_file, on=['coord','chr','start','stop','name','score','strand','source','feature','frame','info','gene_lenght'], how = 'outer')
    #tpm_profile_all_sub = tpm_profile_all.loc[~(tpm_profile_all.iloc[:, 12:]==0).all(axis=1)]

#tpm_profile_all = metaT_file.append(metaT_file, ignore_index=True, sort=False)

#tpm_profile_all.to_csv(wd + '/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p95.TPM.original.csv', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)
#tpm_profile_all_sub.to_csv(wd + '/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p95.TPM.csv', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)

tpm_profile_all.to_csv(wd + '/output/data_integration/' + map_reference + '/' + omics + '.genes_abundances.p' + identity + '.' + normalization + '.csv', index=False, encoding='utf-8', sep="\t", quoting=csv.QUOTE_NONE)
