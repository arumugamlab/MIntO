###################################################
# List the nodes in your cluster as a python array.
# For example, if you have four nodes:
#   leanserver1, leanserver2, fatserver1, fatserver2
# define it as:
#
# CLUSTER_NODES = ['leanserver1', 'leanserver2', 'fatserver1', 'fatserver2']
#
# This is a python variable definition, so you can define patterns.
# If you have 5 nodes with uniform naming scheme like this:
#
# CLUSTER_NODES = ['bladex01', 'bladex02', 'bladex03', 'bladex04', 'bladex05']
#
# you can also define it as:
#
# CLUSTER_NODES = [f"bladex{id}" for id in range(1,6)]
#
###################################################
CLUSTER_NODES = None

###################################################
# List the directory that is located on the local disk on each node.
# This directory MUST be present in all nodes.
# User should have write permissions in this directory.
# Example:
# CLUSTER_LOCAL_DIR = "/scratch/MIntO/mirror"
###################################################
CLUSTER_LOCAL_DIR = None

###################################################
# Name the cluster workload manager (or) queuing
# system. This will ensure that file distribution
# to the nodes can take place effectively using
# the right argument to the batch submission
# commmand.
# Currently handles (case-insensitively):
#  slurm, torque, PBS, SGE, LSF
###################################################
CLUSTER_WORKLOAD_MANAGER = None
