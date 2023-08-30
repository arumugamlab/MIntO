#!/usr/bin/env python

'''
Ensure that key fields or variables are defined.

Authors: Mani Arumugam
'''

# Make sure that shadow-prefix is specified.

if workflow.shadow_prefix is None:
    raise Exception("MIntO error: shadow_prefix needs to be defined using '--shadow-prefix'")
