#!/usr/bin/env python

'''
Ensure that key fields or variables are defined.

Authors: Mani Arumugam
'''

# Make sure that shadow-prefix is specified.

if '--shadow-prefix' in sys.argv:
    print("NOTE: Using SHADOW_PREFIX={}".format(sys.argv[1 + sys.argv.index('--shadow-prefix')]))
else:
    raise Exception("MIntO error: shadow_prefix needs to be defined in Snakemake commandline using '--shadow-prefix'")
