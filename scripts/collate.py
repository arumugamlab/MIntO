#!/usr/bin/python

# '''
# integrates the different gene annotations from eggnog-mapper, dbcan and kofamscan

# Authors: Vithiagaran Gunalan
# '''

import pandas as pd
import numpy as np
import sys
import os
import sys

mydir = sys.argv[1]
name = sys.argv[2]

eggnog = pd.DataFrame()
dbcan  = pd.DataFrame()
kofam  = pd.DataFrame()

eggnogfile = mydir+"/"+name+"_eggNOG.tsv"
if os.path.exists(eggnogfile):
    eggnog = pd.read_csv(eggnogfile, sep="\t", index_col=False)
    # Prefix all columns with 'eggNOG.'
    eggnog = eggnog.rename(lambda x: f'eggNOG.{x}', axis='columns')
    # Fix some names and remove prefix from ID
    eggnog = eggnog.rename(columns={"eggNOG.KEGG_ko"    : "eggNOG.KEGG_KO",
                                    "eggNOG.ID"         : "ID",
                                    "eggNOG.eggNOG_OGs" : "eggNOG.OGs"})
    eggnog.set_index("ID", inplace=True)
    df = eggnog 

dbcanfile = mydir+"/"+name+"_dbCAN.tsv"
if os.path.exists(dbcanfile):
    dbcan = pd.read_csv(dbcanfile, sep="\t", index_col=False)
    # Prefix eCAMI columns with 'dbCAN.'
    dbcan = dbcan.rename(columns={"eCAMI.submodule" : "dbCAN.eCAMI_submodule",
                                  "eCAMI.subfamily" : "dbCAN.eCAMI_subfamily"})
    dbcan.set_index("ID", inplace=True)
    df = dbcan 

kofamfile = mydir+"/"+name+"_kofam.tsv"
if os.path.exists(kofamfile):
    kofam = pd.read_csv(kofamfile, sep="\t", index_col=False)
    # Prefix all columns with 'kofam.'
    kofam = kofam.rename(columns={"kofam_KO"      : "kofam.KEGG_KO",
                                  "kofam_Module"  : "kofam.KEGG_Module",
                                  "kofam_Pathway" : "kofam.KEGG_Pathway"})
    kofam.set_index("ID", inplace=True)
    df = kofam 

# eggnog and kofam
if (not eggnog.empty and not kofam.empty):
    KOdf = pd.merge(eggnog,kofam, on="ID", how="outer").fillna("-")
    eggnoglist = KOdf['eggNOG.KEGG_KO'].tolist()
    kofamlist = KOdf['kofam.KEGG_KO'].tolist()
    newlist = []
    final = []
    for i in range(len(kofamlist)):
        mydict = {}
        myline = eggnoglist[i]+","+kofamlist[i]
        mylist = myline.split(",")
        for i in mylist:
            mydict[i] = ""
        newlist = list(mydict.keys())
        mystring = ""
        if len(newlist) > 1:
            if "-" in newlist:
                newlist.remove("-")
            mystring = ",".join(newlist)
            final.append(mystring)
        else:
            mystring = "".join(newlist)
            final.append(mystring)
    KOdf["merged.KEGG_KO"] = final
    df = KOdf
# eggnog and dbcan, no kofam
elif (not eggnog.empty and not dbcan.empty and kofam.empty):
    del df
    df = pd.merge(eggnog,dbcan, on="ID", how="outer").fillna("-")
# kofam and dbcan, no eggnog
elif (not kofam.empty and not dbcan.empty and eggnog.empty):
    del df
    df = pd.merge(dbcan,kofam, on="ID", how="outer").fillna("-")

# eggnog, kofam and dbcan
if (not eggnog.empty and not kofam.empty and not dbcan.empty):
    del df
    df = pd.merge(KOdf,dbcan, on="ID", how="outer").fillna("-")

df.to_csv(sys.stdout, sep="\t")
