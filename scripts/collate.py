#!/usr/bin/python

# '''
# integrates the different gene annotations from eggnog-mapper, dbcan and kofamscan

# Authors: Vithiagaran Gunalan
# '''

import pandas as pd
import numpy as np
import sys

mydir = sys.argv[1]
name = sys.argv[2]

eggnogfile = mydir+"/"+name+".eggNOG5.annotations"
eggnog = pd.read_csv(eggnogfile, sep="\t", index_col=False)
eggnog.set_index("ID", inplace=True)

dbcanfile = mydir+"/"+name+"_dbCAN.tsv"
dbcan = pd.read_csv(dbcanfile, sep="\t", index_col=False)
dbcan.set_index("ID", inplace=True)

kofamfile = mydir+"/"+name+"_kofam.tsv"
kofam = pd.read_csv(kofamfile, sep="\t", index_col=False)
kofam.set_index("ID", inplace=True)

KOdf = pd.merge(eggnog,kofam, on="ID", how="outer").fillna("-")
kegglist = KOdf['KEGG_ko'].tolist()
kofamlist = KOdf['kofam_KO'].tolist()
newlist = []
final = []

for i in range(len(kofamlist)):
    mydict = {}
    myline = kegglist[i]+","+kofamlist[i]
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
        
KOdf["merged_KO"] = final
df = pd.merge(KOdf,dbcan, on="ID", how="outer").fillna("-")

myoutfile = mydir+"/"+name+".annotations.tsv"

df.to_csv(myoutfile, sep="\t")
