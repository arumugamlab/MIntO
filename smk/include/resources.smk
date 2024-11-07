#!/usr/bin/env python

'''
Helper functions to estimate workflow resources

Authors: Mani Arumugam
'''

#############################################
# File size
#############################################

# If file exists, return its size.
# If not, return 1GB - this is placeholder for allowing --dry-run to work
def get_file_size(f):
    import os.path
    if os.path.exists(f):
        return(os.path.getsize(f))
    else:
        return(1<<9)

#############################################
# TSV table sizes
#############################################

# If file exists, return TSV table dimensions
# If not, return 2^27 so that memsize will be 1GB at 8 bytes per cell - this is placeholder for allowing --dry-run to work
def get_tsv_dimensions(filename):
    import os.path
    if not os.path.exists(filename):
        return([1, 1<<27])

    import pandas as pd
    num_columns = len(pd.read_csv(filename, sep='\t', nrows=0).columns)
    with open(filename) as f:
        num_rows = sum(1 for line in f)
    return([num_rows, num_columns])

def get_tsv_cells(filename):
    dim = get_tsv_dimensions(filename)
    return(dim[0]*dim[1])
