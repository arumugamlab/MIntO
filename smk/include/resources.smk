#!/usr/bin/env python

'''
Helper functions to estimate workflow resources

Authors: Mani Arumugam
'''

# If file exists, return its size.
# If not, return 1GB - this is placeholder for allowing --dry-run to work
def get_file_size(f):
    import os.path
    if os.path.exists(f):
        return(os.path.getsize(f))
    else:
        return(1<<9)

# TSV file dimensions
def get_tsv_dimensions(filename):
    import pandas as pd
    num_columns = len(pd.read_csv(filename, sep='\t', nrows=0).columns)
    with open(filename) as f:
        num_rows = sum(1 for line in f)
    return([num_rows, num_columns])

def get_tsv_cells(filename):
    dim = get_tsv_dimensions(filename)
    return(dim[0]*dim[1])
