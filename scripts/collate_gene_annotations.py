#!/usr/bin/env python

# '''
# integrates the different gene annotations from eggnog-mapper, dbcan and kofamscan

# Authors: Vithiagaran Gunalan, Mani Arumugam, Judit Szarvas
# '''

import pandas as pd
import sys
import os

eggnog = pd.DataFrame()
dbcan  = pd.DataFrame()
kofam  = pd.DataFrame()
df  = pd.DataFrame()

input_annotations = sys.argv[1:]
invalid_inputs = []
for filename in input_annotations:

    # Make sure file exists
    if os.path.exists(filename) == False or os.path.getsize(filename) == 0:
        raise Exception("File {} not found!".format(filename))

    # eggNOG annotations
    if filename.endswith('_eggNOG.tsv'):
        eggnog = pd.read_csv(filename, sep="\t", index_col=False)
        # Prefix all columns with 'eggNOG.'
        eggnog = eggnog.rename(lambda x: f'eggNOG.{x}', axis='columns')
        # Fix some names and remove prefix from ID
        eggnog = eggnog.rename(columns={"eggNOG.KEGG_ko"    : "eggNOG.KEGG_KO",
                                        "eggNOG.ID"         : "ID",
                                        "eggNOG.eggNOG_OGs" : "eggNOG.OGs"})
        eggnog.set_index("ID", inplace=True)
        df = eggnog

    # dbCAN annotations
    elif filename.endswith('_dbCAN.tsv'):
        dbcan = pd.read_csv(filename, sep="\t", index_col=False)
        # Prefix eCAMI columns with 'dbCAN.'
        dbcan = dbcan.rename(columns={"eCAMI.submodule" : "dbCAN.eCAMI_submodule",
                                      "eCAMI.subfamily" : "dbCAN.eCAMI_subfamily"})
        dbcan.set_index("ID", inplace=True)
        df = dbcan

    # kofam annotations
    elif filename.endswith('_kofam.tsv'):
        kofam = pd.read_csv(filename, sep="\t", index_col=False)
        # Prefix all columns with 'kofam.'
        kofam = kofam.rename(columns={"kofam_KO"      : "kofam.KEGG_KO",
                                      "kofam_Module"  : "kofam.KEGG_Module",
                                      "kofam_Pathway" : "kofam.KEGG_Pathway"})
        kofam.set_index("ID", inplace=True)
        df = kofam
    else:
        invalid_inputs.append(filename)
        print("WARNING:", filename, "not supported", file=sys.stderr)

for fn in invalid_inputs:
    input_annotations.remove(fn)

if len(input_annotations) > 1:
    # eggnog and kofam
    if (not eggnog.empty and not kofam.empty):
        KOdf = pd.merge(eggnog,kofam, on="ID", how="outer").fillna("-")
        eggnoglist = KOdf['eggNOG.KEGG_KO'].tolist()
        kofamlist = KOdf['kofam.KEGG_KO'].tolist()
        final_kos = []
        for i in range(len(kofamlist)):
            uniqmergedlist = eggnoglist[i].split(",")
            for kofam_ko in kofamlist[i].split(","):
                if kofam_ko not in uniqmergedlist:
                    uniqmergedlist.append(kofam_ko)
            mystring = ""
            if len(uniqmergedlist) > 1:
                if "-" in uniqmergedlist:
                    uniqmergedlist.remove("-")
                mystring = ",".join(uniqmergedlist)
            else:
                mystring = uniqmergedlist[0]
            final_kos.append(mystring)
        KOdf["merged.KEGG_KO"] = final_kos
        df = KOdf
        # eggnog, kofam and dbcan
        if not dbcan.empty:
            df = pd.merge(KOdf,dbcan, on="ID", how="outer").fillna("-")
    # eggnog and dbcan, no kofam
    elif (not eggnog.empty and not dbcan.empty and kofam.empty):
        df = pd.merge(eggnog,dbcan, on="ID", how="outer").fillna("-")
    # kofam and dbcan, no eggnog
    elif (not kofam.empty and not dbcan.empty and eggnog.empty):
        df = pd.merge(dbcan,kofam, on="ID", how="outer").fillna("-")
if not df.empty:
    df.to_csv(sys.stdout, sep="\t")
else:
    print("WARNING: no supported input", file=sys.stderr)
