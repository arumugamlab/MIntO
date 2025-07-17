"""
Step 7:
    This script will replace 7-run-retrieving-score-genomes.py. 
    The script calcuate a set of statistics for the genomes in a folder
        
        input: 
            a folder containing fasta files (or maybe just a fasta file)
        
        output:
            a table with the statistics
    gaps : TGTTGCTTCTCNNNNNNNNNNTTCAATACTTTCTCTA
    
    I will still calculate the old score.
    But a new score will be calculated with Mani's formula:
"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import os
from os import listdir
from os.path import isfile, join
import warnings
import argparse as ap
import sys
import numpy as np
from itertools import groupby
import math
import pandas as pd
import glob
import time
import math


def read_params():
	p = ap.ArgumentParser( description = ( 
		"""Given a folder of genomes and a checkm output it calculates the score for the genomes.
		The score can be calculated in two ways : checkm and genome [--score_method].
		The option checkm calculate the score using checkm output, while genome only on the genome quality.
		"""),
		formatter_class = ap.ArgumentDefaultsHelpFormatter )


	p.add_argument("--checkm_output", 
		help = "The checkm output",
		)
	
	p.add_argument("--fasta_folder", 
		help = "The folder with all the genomes")

	p.add_argument("--output_file",
		help = "Where to put the output file")

	p.add_argument("--score_method", 
		help = "The score can be calculated based on checkm results [checkm] [defualt] or with a more complex formula based on the assembly [genome].\n checkm is suggested when assembly has been done only with illumina. \n genome when assembly takes care also nanopore",
		default= "checkm",
		choices= ["checkm", "genome"])

	p.add_argument("--gap_penalty", 
		help = "Penalty for the number of gaps in the sequence [only if --score_method genome]", 
		default = 0.15)

	p.add_argument('--circularity_threshold',
		help = "Proportion of the genome that should be captured from the a contig to be considered circular [only if --score_method genome]",
		default= 0.9)

	p.add_argument("--favor_circular", 
		help = "The score in --score_method will be updated by adding the circularity fracion * number of circular bases, meaning an higher score will be given to circular genomes",
		choices= ["yes", "no"], 
		default= "yes")


	return(p.parse_args())



def printtime():
	import time
	return(time.strftime("%H:%M:%S", time.localtime()))


def is_fasta(filename):
	"""
	Check if the file is a fasta
	"""
	with open(filename, "r") as f:
		fasta = SeqIO.parse(f, "fasta")
		return (any(fasta))


def get_number_of_gaps(sequence):
	"""
	Given a DNA string, it parse it to find out how many gaps are inside
	"""
	count_of_gaps = 0
	# find the first occurrence
	find_n = sequence.find("N")
	if find_n != -1 :
		count_of_gaps = count_of_gaps + 1
		for i in range(find_n, len(sequence)):
			if sequence[i] != "N":
				break
		
		if i + 1 == len(sequence):
			return(count_of_gaps)
		count_of_gaps  = count_of_gaps +  get_number_of_gaps(sequence[i:])
	return(count_of_gaps)




def calculate_N_L_statistics(list_of_contigs_length, threshold = 50):
	"""
	Given a list of lenght of contigs, it calculates the N50 etc statistics
	"""

	# total length
	total_length = sum(list_of_contigs_length)

	# sort the list descending
	list_of_contigs_length.sort(reverse = True)

	# deciding the of the genomes
	perc_length = total_length * threshold /100

	# starting from the end of the list, 
	count  = 0
	for i in range(len(list_of_contigs_length)):
		
		# summing the bp
		count = count + list_of_contigs_length[i]

		if count >= perc_length:
			return(list_of_contigs_length[i], i + 1 ) #N50, L50 




def calculate_fasta_statistics_checkm_based(fasta_file):

	"""
	The function takes a fastafiles in input and calculate the statistics needed 
	Bin_Id, Contigs, Bases, L50, L75, L90, N50_bp, N75_bp, N90_bp, Longest_bp, SeqScore, seq_1M, seq_2M
	It returns a dictionary with all these information (computationally not really wise since I can open a file and make it write it at the moment)
	"""

	# get the file name
	bin_id = os.path.basename(fasta_file)

	print("[{}] Calculating statistics for {}".format(printtime(), bin_id))
	#

	total = 0
	contigs = 0
	length_of_contigs = [] # useful to calculate entropy later 

	gc_content_total = 0 # this should be divided by the total

	longest_contig_name = ""
	longest_contig_size = 0
	shortest_contig_name = ""
	shortest_contig_size = len(list(SeqIO.parse(fasta_file, "fasta"))[0]) # taking the first contig for reference

	N_statistics_L_statistics_calculation = [] # list of length of contigs

	# defining sequences with length > 1,000,000 and 2,000,000 for calculating the score 
	seq_1M = 0
	seq_2M = 0

	with open(fasta_file) as handle:

		for record in SeqIO.parse(handle, "fasta"):
			
			contigs = contigs + 1
			id = record.id

			seq = record.seq
			seq_len = len(record.seq)
			length_of_contigs.append(seq_len)
			total = total + seq_len
			
			# updating the number of seqeunce higher than 1,000,000 and 2,000,000
			if seq_len > 1000000 :
				seq_1M = seq_1M + 1
			if seq_len > 2000000: 
				seq_2M = seq_2M + 1
			

			# adding the contig length in the N and L 
			N_statistics_L_statistics_calculation.append(seq_len)

			# updating the longest contig
			if seq_len >= longest_contig_size: 
				longest_contig_name = id
				longest_contig_size = seq_len 

			# updating the shortest contig
			if seq_len <= shortest_contig_size:
				shortest_contig_name = id
				shortest_contig_size = seq_len

			gc_content_total = gc_content_total + seq_len*GC(seq) # this should be divided by the length of genome

	statistics_dictionary = {"L50" : 0,
		"L75" : 0, 
		"L90" : 0, 
		"N50" : 0, 
		"N75" : 0 , 
		"N90" : 0 }

	for stat in [50, 75, 90]:
		statistics = calculate_N_L_statistics(N_statistics_L_statistics_calculation, threshold=stat)
		statistics_dictionary["N{}".format(stat)] = statistics[0]
		statistics_dictionary["L{}".format(stat)] = statistics[1]

	info_dict = {
		"Bin_Id" : bin_id,
		"Contigs" : contigs,
		"Bases" : total, 
		"GC" : gc_content_total/total,
		"L50" : statistics_dictionary["L50"], 
		"L75" : statistics_dictionary["L75"], 
		"L90" : statistics_dictionary["L90"], 
		"N50_bp" : statistics_dictionary["N50"], 
		"N75_bp" : statistics_dictionary["N75"], 
		"N90_bp" : statistics_dictionary["N90"], 
		"Longest_bp" : longest_contig_size,
		"Longest_contig_name": longest_contig_name,
		"Shortest_bp" : shortest_contig_size, 
		"Shortest_contig_name" : shortest_contig_name,
		"seq_1M": seq_1M, 
		"seq_2M" : seq_2M
		}

	return(info_dict)




def calculate_fasta_statistics_genome_based(fasta_file, circularity_threshold = 0.9, gap_penalty = 0.15 , favor_circular = "yes"):

	"""
	The function takes a fastafiles in input and calculate the statistics needed 
	"Bin_Id"
	"Contigs"
	"Bases"
	"L50"
	L75"
	L90"
	"N50_bp"
	"N75_bp"
	"N90_bp"
	"Longest_bp"
	"SeqScore"
	"Gaps"
	"GapLengths"
	"Circular"
	"CircularFraction"
	"N1Mb"
	"N2Mb"
	"Circular_chromosome"
	It returns a dictionary with all these information (computationally not really wise since I can open a file and make it write it at the moment)
	"""


	# get the file name
	bin_id = os.path.basename(fasta_file)

	print("[{}] Calculating statistics for {}".format(printtime(), bin_id))

	total = 0
	contigs = 0
	length_of_contigs = [] # useful to calculate entropy later 

	gc_content_total = 0 # this should be divided by the total

	# gaps information
	number_of_gaps = 0
	total_length_of_gaps = 0

	# circular information
	how_many_contigs_considered_circular = 0
	name_of_circular_contigs = []
	circular_contigs_length = []

	longest_contig_name = ""
	longest_contig_size = 0
	shortest_contig_name = ""
	shortest_contig_size = len(list(SeqIO.parse(fasta_file, "fasta"))[0]) # taking the first contig for reference

	N_statistics_L_statistics_calculation = [] # list of length of contigs

	# defining sequences with length > 1,000,000 and 2,000,000 for calculating the score 
	seq_1M = 0
	seq_2M = 0

	# defining if a genome is circular if 90% of the genomes is inside one contig
	circular = "no"


	with open(fasta_file) as handle:

		for record in SeqIO.parse(handle, "fasta"):
			
			contigs = contigs + 1
			id = record.id

			seq = record.seq
			seq_len = len(record.seq)
			length_of_contigs.append(seq_len)
			total = total + seq_len
			
			if "circular" in id:
				how_many_contigs_considered_circular += 1
				circular_contigs_length.append(seq_len)
				name_of_circular_contigs.append(id)
			
			# updating the number of seqeunce higher than 1,000,000 and 2,000,000
			if seq_len > 1000000 :
				seq_1M = seq_1M + 1
			if seq_len > 2000000: 
				seq_2M = seq_2M + 1
			

			# adding the contig length in the N and L 
			N_statistics_L_statistics_calculation.append(seq_len)

			# updating the longest contig
			if seq_len >= longest_contig_size: 
				longest_contig_name = id
				longest_contig_size = seq_len 

			# updating the shortest contig
			if seq_len <= shortest_contig_size:
				shortest_contig_name = id
				shortest_contig_size = seq_len

			gc_content_total = gc_content_total + seq_len*GC(seq) # this should be divided by the number of contigs

			# check if there are gaps
			if "N" in seq:
				# count how long are the gaps 
				total_length_gaps_contigs = seq.count("N")
				total_length_of_gaps = total_length_of_gaps + total_length_gaps_contigs
				number_of_gaps_contig =  get_number_of_gaps(seq)
				number_of_gaps = number_of_gaps + number_of_gaps_contig

	statistics_dictionary = {"L50" : 0,
		"L75" : 0, 
		"L90" : 0, 
		"N50" : 0, 
		"N75" : 0 , 
		"N90" : 0 }

	for stat in [50, 75, 90]:
		statistics = calculate_N_L_statistics(N_statistics_L_statistics_calculation, threshold=stat)
		statistics_dictionary["N{}".format(stat)] = statistics[0]
		statistics_dictionary["L{}".format(stat)] = statistics[1]


	# defining if the genomes is circular
	# check if the longest contig is also in the list of circular
	longest_contig_is_circular = "no"
	if (longest_contig_name  in  name_of_circular_contigs ):
		longest_contig_is_circular = "yes"
	
	if sum(circular_contigs_length) >= ( circularity_threshold * total) :
		circular = "yes"
	

	# calculating the score
	entropy = 0
	for slen in length_of_contigs:
		p = slen/total
		p_log = math.log10(slen/ total)
		entropy = entropy - (p * p_log)
	

	#############################################################
	# This probably can be written inside reading the fasta in order to avoid the loop here
	initial_length = 0
	circular_chr = 0
	# finding the best circular chromosome
	for i in range(len(name_of_circular_contigs)):
		c_contigs_name = name_of_circular_contigs[i]
		c_contigs_length = circular_contigs_length[i]

		if c_contigs_length > initial_length: # only working on the longest one
			initial_length = c_contigs_length

			if c_contigs_length > circularity_threshold * total:
				if "_circularU" in c_contigs_name:
					circular_chr = 3.2
				elif "_circularT" in c_contigs_name:
					circular_chr = 3.15
				elif "_circularA_circularL" in c_contigs_name:
					circular_chr = 3
				elif "_circularL" in c_contigs_name:
					circular_chr = 2
				else:
					circular_chr = 1 
	############################################################
	
	circular_fraction = sum(circular_contigs_length)/total

	# L90 and N90 works for MAGs with long contigs overall, e.g. short+long-read hybrid assemblies.
	# L75 and N75 works for MAGs with shorter contigs overall, e.g. short-read-only assemblies.
	L_stat = statistics_dictionary["L75"]
	N_stat = statistics_dictionary["N75"]

	score = math.log10(longest_contig_size) + math.log10(N_stat) - math.log10(total) - math.log10(L_stat) - entropy
	score = score + (seq_1M * 0.1)
	score = score + (seq_2M * 0.1)
	score = score + (circular_chr * circular_fraction)
	score = score - (gap_penalty * number_of_gaps)
	score = round(score, 6)

	# if favor-circular was given, then adjust score using fraction-circular


	if favor_circular == "yes":
		score = score + (circular_fraction*sum(circular_contigs_length))



	info_dict = {
		"Bin_Id" : bin_id,
		"Contigs" : contigs,
		"Bases" : total, 
		"GC" : gc_content_total/total,
		"L50" : statistics_dictionary["L50"], 
		"L75" : statistics_dictionary["L75"], 
		"L90" : statistics_dictionary["L90"], 
		"N50_bp" : statistics_dictionary["N50"], 
		"N75_bp" : statistics_dictionary["N75"], 
		"N90_bp" : statistics_dictionary["N90"], 
		"Longest_bp" : longest_contig_size,
		"Longest_contig_name": longest_contig_name,
		"Longest_contig_circular": longest_contig_is_circular, 
		"Shortest_bp" : shortest_contig_size, 
		"Shorterst_contig_name" : shortest_contig_name, 
		"Gaps" : number_of_gaps, 
		"GapLengths" : total_length_of_gaps, 
		"CircularContigs" : name_of_circular_contigs, 
		"CircularFraction" : circular_fraction,
		"Circular": circular_contigs_length,
		"ConsideredCircular" : circular,
		"Score" : score , # this is the sequence score , not the real score! We have to use it in the calculation for the score with checkm
		"seq_1M" : seq_1M, 
		"seq_2M" : seq_2M
		}

	return(info_dict)

############################# Reading the parameters and performing the functions

args = read_params()
folder_fasta_file = args.fasta_folder
score_method = args.score_method
gap_penalty = args.gap_penalty
circularity_threshold = args.circularity_threshold
favor_circular = args.favor_circular
path_to_checkm_table = args.checkm_output
checkm = pd.read_csv(path_to_checkm_table, sep = "\t", index_col = "Name") #index_col = "Bin Id"
output_file = args.output_file




fna_files = glob.glob(folder_fasta_file + "/*.fna")
fa_files = glob.glob(folder_fasta_file + "/*.fa")
folder_fasta_file = fna_files + fa_files


print("Folder Fasta files: {}".format(folder_fasta_file) )
print("Score Method: {}".format(score_method) )

# if checkm score is decided:
if score_method == "checkm":

	print("[calculate_score_genomes] --score_method: checkm")

	# open the output file
	with open(output_file, "w") as fh:
		fh.write("#--score_method {}\n".format(score_method))
		fh.write("Bin_id\tScore\tSeq_score\tQual_score\tCompleteness\tContamination\tContigs\tBases\tGC\tL50\tL75\tL90\tN50_bp\tN75_bp\tN90_bp\tLongest_bp\tShortest_bp\tStrain_heterogeneity\tseq_1M\tseq_2M\n")
	
		# calculate the statistics 
		for file in folder_fasta_file:
			if is_fasta(file): # if the file is a fasta
				genome_info = calculate_fasta_statistics_checkm_based(file) # information of the genome
				# split the genomes name in order to have the same in the checkm table as a key
				genome_name = file.split("/")[-1]
				if genome_name[-4:] == ".fna":
					genome_name = genome_name.replace(genome_name[-4:], "")
				elif genome_name[-3:] == ".fa":
					genome_name = genome_name.replace(genome_name[-3:], "")
				
				print("[calculate_score_genomes] Score for {}".format(genome_name))
				
				# found information in checkm file
				if genome_name in checkm.index:
					completeness = float(checkm["Completeness"].loc[genome_name])
					contamination = float(checkm["Contamination"].loc[genome_name])

					# calculate the score (based on genome)
					seq_score = math.log10(genome_info["Longest_bp"]/int(genome_info["Contigs"])) + math.log10(genome_info["N50_bp"]/int(genome_info["L50"]))
					qual_score = completeness - (2 * contamination)
					score = (0.1 * qual_score) + seq_score

					fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(genome_name, 
						score, 
						seq_score, 
						qual_score, 
						completeness, 
						contamination, 
						genome_info["Contigs"],
						genome_info["Bases"], 
						genome_info["GC"], 
						genome_info["L50"], 
						genome_info["L75"], 
						genome_info["L90"],
						genome_info["N50_bp"], 
						genome_info["N75_bp"], 	
						genome_info["N90_bp"], 
						genome_info["Longest_bp"],
						genome_info["Shortest_bp"], 
						genome_info["seq_1M"], 
						genome_info["seq_2M"]
						))
				else:
					print("There is a problem in calculate_genomes_score.py, check if the the genomes name are identical between checkm and the dictionary in memory!") # probably this can be done in an only passage (open the fasta and calcuate and write immediately)
					sys.exit()

				

# if genome score is decided:				
elif score_method == "genome":

	print("[calculate_score_genomes] --score_method: genome")

	# open the output file
	with open(output_file, "w") as fh:
		fh.write("#--score_method {}\n".format(score_method))
		fh.write("Bin_id\tScore\tSeqScore\tQualScore\tCompleteness\tContamination\tContigs\tBases\tGC\tL50\tL75\tL90\tN50_bp\tN75_bp\tN90_bp\tLongest_bp\tShortest_bp\tGaps\tGapLengths\tCircularContigs\tCircularFraction\tCircular\tConsideredCircular\tStrain_heterogeneity\tseq_1M\tseq_2M\n")
	
		# calculate the statistics 
		for file in folder_fasta_file:
			if is_fasta(file): # if the file is a fasta
				genome_info = calculate_fasta_statistics_genome_based(file) # information of the genome
				# split the genomes name in order to have the same in the checkm table as a key
				genome_name = file.split("/")[-1]
				if genome_name[-4:] == ".fna":
					genome_name = genome_name.replace(genome_name[-4:], "")
				elif genome_name[-3:] == ".fa":
					genome_name = genome_name.replace(genome_name[-3:], "")

				# found information in checkm file
				if genome_name in checkm.index:
					completeness = float(checkm["Completeness"].loc[genome_name])
					contamination = float(checkm["Contamination"].loc[genome_name])

					seq_score = float(genome_info["Score"]) # the seqeunce score has been calculated before 
					qual_score = completeness - (2 * contamination)
					score_checkm = (0.1 * qual_score) + seq_score



					fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
					genome_name, 
					score_checkm, # this is updated with the checkm  (probably I can make just one function now that I understand)
					seq_score, 
					qual_score,
					completeness, 
					contamination, 
					genome_info["Contigs"],
					genome_info["Bases"], 
					genome_info["GC"], 
					genome_info["L50"], 
					genome_info["L75"], 
					genome_info["L90"],
					genome_info["N50_bp"], 
					genome_info["N75_bp"], 	
					genome_info["N90_bp"], 
					genome_info["Longest_bp"],
					genome_info["Shortest_bp"], 
					genome_info["Gaps"], #column added in genome modality
					genome_info["GapLengths"],  #column added in genome modality
					genome_info["CircularContigs"],  #column added in genome modality
					genome_info["CircularFraction"], #column added in genome modality
					genome_info["Circular"],  #column added in genome modality
					genome_info["ConsideredCircular"],  #column added in genome modality
					genome_info["seq_1M"], 
					genome_info["seq_2M"]
					  ))
				else:
					print("There is a problem in calculate_genomes_score.py, check if the the genomes name are identical between checkm and the dictionary in memory!") # probably this can be done in an only passage (open the fasta and calcuate and write immediately)
					sys.exit()

				# calculate the score
				#seq_score = math.log10(genome_info["Longest_bp"])/int(genome_info["Contigs"])) + math.log10(genome_info["N50_bp"])/int(genome_info["L50"]))
				#qual_score = completeness - (2 * contamination)
				#score (0.1 * qual_score) + seq_score

				
else:
	print("Which one did you mean? checkm or genome? Check if you type them right!")
	sys.exit()
