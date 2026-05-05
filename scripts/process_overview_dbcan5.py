#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv) < 2:
    print("program <input overview.tsv> <output prefix>")
    exit()

input_file = sys.argv[1]
output_prefix = sys.argv[2]

colnames = ['ID', 'EC_nums', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', 'NumberofTools', 'RecommendResults', 'Substrate']
df = pd.read_csv(input_file, header=0, names=colnames, sep = "\t")

# process only those with min 2 tools agree
df = df[df.NumberofTools > 1]

# separate the multi-domain hits into list-in-cells
df["dbCAN.hit"]=df["RecommendResults"].str.split("|")
df["dbCAN.hit_EC"]=df["EC_nums"].str.split("|")

# separate list-in-cell into rows for dbCAN.hit
df1 = df.explode("dbCAN.hit").reset_index(drop=True)
df1["rownum"] = df1.groupby("ID")["dbCAN.hit"].cumcount().add(1)
df1.drop(columns = ['EC_nums','RecommendResults','dbCAN.hit_EC'], inplace=True)

# separate list-in-cell into rows for EC numbers
df2 = df.explode("dbCAN.hit_EC").reset_index(drop=True)
df2["rownum"] = df2.groupby("ID")["dbCAN.hit_EC"].cumcount().add(1)
df2 = df2[["ID", "dbCAN.hit_EC", "rownum"]]

# in theory each domain separated by | have an EC number so df1 and df2 should have corresponding amount of rows
df = df1.merge(df2, on = ['ID', 'rownum'])

# CBM will be bindingmodule, rest dbCAN_sub
df["dbCAN.binding_module"] = '-'
df.loc[df['dbCAN.hit'].str.contains("CBM"), 'dbCAN.binding_module'] = df['dbCAN.hit']
df["dbCAN.dbCAN_sub"] = '-'
df.loc[df['dbCAN.hit'].str.contains("CBM") == False, 'dbCAN.dbCAN_sub'] = df['dbCAN.hit']
# remove the eCAMI sub-clustering numbers
df['dbCAN.dbCAN_sub.subfamily'] = df['dbCAN.dbCAN_sub'].str.replace(r'_e\d+', '', regex=True)

if (output_prefix != "-"):
    # output subfamily file and substrate separately
    output_file = f"{output_prefix}.cazymes.tsv"
    df.to_csv(output_file, columns = ["ID", "dbCAN.hit", "dbCAN.hit_EC", "dbCAN.binding_module", "dbCAN.dbCAN_sub", "dbCAN.dbCAN_sub.subfamily"], sep = "\t", index=False)

    output_substr = f"{output_prefix}.substrate.tsv"
    df = df[['ID', 'Substrate']]
    df[df.Substrate != "-"].drop_duplicates().to_csv(output_substr, sep = "\t", index=False)
else:
    # stdout and only cazymes
    df.to_csv(sys.stdout, columns = ["ID", "dbCAN.hit", "dbCAN.hit_EC", "dbCAN.binding_module", "dbCAN.dbCAN_sub", "dbCAN.dbCAN_sub.subfamily"], sep = "\t", index=False)
