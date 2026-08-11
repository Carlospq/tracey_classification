from Bio import AlignIO
from prody import *
from pylab import *
from sklearn.metrics.cluster import adjusted_mutual_info_score

import numpy as np
import os, subprocess, sys, math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inMSA', dest='MSA', required=True, help='Input multisequence alignment file')
parser.add_argument('-t', '--MIthreshold', dest='cut', default=0.6, required=False, help='MI threshold')
args = parser.parse_args()

MIpath = "coevolution/MI/"
MSA = args.MSA
cut = args.cut

################################################################################
# FUNCTIONS
################################################################################
# Jensen-Shannon Divergence (MSA conservation)
################################################################################
def js_divergence(col, sim_matrix, bg_distr, seq_weights, pseudocount, amino_acids, gap_penalty=1):
	""" Return the Jensen-Shannon Divergence for the column with the background
	distribution bg_distr. sim_matrix is ignored. JSD is the default method."""
	distr = bg_distr[:]
	fc = weighted_freq_count_pseudocount(col, seq_weights, pseudocount, amino_acids)
	# if background distrubtion lacks a gap count, remove fc gap count
	if len(distr) == 20:
		new_fc = fc[:-1]
		s = sum(new_fc)
		for i in range(len(new_fc)):
			new_fc[i] = new_fc[i] / s
		fc = new_fc
	if len(fc) != len(distr): return -1
	# make r distriubtion
	r = []
	for i in range(len(fc)):
		r.append(.5 * fc[i] + .5 * distr[i])
	d = 0.
	for i in range(len(fc)):
		if r[i] != 0.0:
			if fc[i] == 0.0:
				d += distr[i] * math.log(distr[i]/r[i], 2)
			elif distr[i] == 0.0:
				d += fc[i] * math.log(fc[i]/r[i], 2)
			else:
				d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)
	# d /= 2 * math.log(len(fc))
	d /= 2
	if gap_penalty == 1:
		return d * weighted_gap_penalty(col, seq_weights)
	else:
		return d


def jaccard(a, b):
	a = set(a)
	b = set(b)
	intersection = len(a.intersection(b))
	union = len(a) + len(b) - intersection
	return (intersection/union)
################################################################################
# Read/Load/Calculate input files
################################################################################
amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"]
aa_to_index = {}
for i, aa in enumerate(amino_acids):
    aa_to_index[aa] = i


def read_scoring_matrix(sm_file, aa_to_index = aa_to_index):
	""" Read in a scoring matrix from a file, e.g., blosum80.bla, and return it as an array. """
	aa_index = 0
	first_line = 1
	row = []
	list_sm = []
	try:
		matrix_file = open(sm_file, 'r')
		for line in matrix_file:
			if line.startswith("#"): continue
			if first_line:
				first_line = 0
				for c in line.split():
					aa_to_index[c.lower()] = aa_index
					#amino_acids.append(c.lower())
					aa_index += 1
			else:
				row = line.split()
				list_sm.append(row)
	except IOError:
		print("Could not load similarity matrix: %s. Using identity matrix..." % sm_file)
		return identity(20)
	# if matrix is stored in lower tri form, copy to upper
	if len(list_sm[0]) < 20:
		for i in range(0,19):
			for j in range(i+1, 20):
				list_sm[i].append(list_sm[j][i])
	for i in range(len(list_sm)):
		for j in range(len(list_sm[i])):
			list_sm[i][j] = float(list_sm[i][j])
	return list_sm


def read_fasta_alignment(filename):
	""" Read in the alignment stored in the FASTA file, filename. Return two lists: the identifiers and sequences. """
	f = open(filename)
	names = []
	alignment = []
	cur_seq = ''
	for line in f:
		line = line[:-1]
		if len(line) == 0: continue
		if line[0] == ';': continue
		if line[0] == '>':
			names.append(line[1:].replace('\r', ''))
			if cur_seq != '':
				cur_seq = cur_seq.upper()
				for i, aa in enumerate(cur_seq):
					if aa not in iupac_alphabet:
						cur_seq = cur_seq.replace(aa, '-')
				alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
				cur_seq = ''
		elif line[0] in iupac_alphabet:
			cur_seq += line.replace('\r', '')
	# add the last sequence
	cur_seq = cur_seq.upper()
	for i, aa in enumerate(cur_seq):
		if aa not in iupac_alphabet:
			cur_seq = cur_seq.replace(aa, '-')
	alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
	# Check if alignment format is correct
	if len(alignment) != len(names) or alignment == []:
		print("Unable to parse alignment.\n")
		sys.exit(1)
	# Check length of sequences
	seq_len = len(alignment[0])
	for i, seq in enumerate(alignment):
		if len(seq) != seq_len:
			print("ERROR: Sequences of different lengths: %s (%d) != %s (%d).\n" % (names[0], seq_len, names[i], len(seq)))
			sys.exit(1)
	return names, alignment


def load_sequence_weights(fname):
	""""Read in a sequence weight file f and create sequence weight list. The weights are in the same order as the sequences each on a new line. """
	seq_weights = []
	try:
		f = open(fname)
		for line in f:
			l = line.split()
			if line[0] != '#' and len(l) == 2:
				seq_weights.append(float(l[1]))
	except IOError:
		pass
	return seq_weights


def calculate_sequence_weights(msa, amino_acids):
	""" Calculate the sequence weights using the Henikoff '94 method for the given msa. """
	seq_weights = [0.] * len(msa)
	for i in range(len(msa[0])):
		freq_counts = [0] * len(amino_acids)
		col = []
		for j in range(len(msa)):
			if msa[j][i] != '-': # ignore gaps
				freq_counts[aa_to_index[msa[j][i]]] += 1
		num_observed_types = 0
		for j in range(len(freq_counts)):
			if freq_counts[j] > 0: num_observed_types += 1
		for j in range(len(msa)):
			d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
			if d > 0:
				seq_weights[j] += 1. / d
	for w in range(len(seq_weights)):
		seq_weights[w] /= len(msa[0])
	return seq_weights


def get_column(col_num, alignment):
	"""Return the col_num column of alignment as a list."""
	col = []
	for seq in alignment:
		if col_num < len(seq): col.append(seq[col_num])
	return col
################################################################################
# Frequency Count and Gap Penalty
################################################################################
def gap_percentage(col):
	"""Return the percentage of gaps in col."""
	num_gaps = 0.
	for aa in col:
		if aa == '-': num_gaps += 1
	return num_gaps / len(col)


def weighted_freq_count_pseudocount(col, seq_weights, pseudocount, amino_acids):
	""" Return the weighted frequency count for a column--with pseudocount."""
	# if the weights do not match, use equal weight
	if len(seq_weights) != len(col):
		seq_weights = [1.] * len(col)
	aa_num = 0
	freq_counts = len(amino_acids)*[pseudocount] # in order defined by amino_acids
	for aa in amino_acids:
		for j in range(len(col)):
			if col[j] == aa:
				freq_counts[aa_num] += 1 * seq_weights[j]
		aa_num += 1
	for j in range(len(freq_counts)):
		freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pseudocount)
	return freq_counts


def weighted_gap_penalty(col, seq_weights):
	""" Calculate the simple gap penalty multiplier for the column. If the
	sequences are weighted, the gaps, when penalized, are weighted
	accordingly. """
	# if the weights do not match, use equal weight
	if len(seq_weights) != len(col):
		seq_weights = [1.] * len(col)
	gap_sum = 0.
	for i in range(len(col)):
		if col[i] == '-':
			gap_sum += seq_weights[i]
	return 1 - (gap_sum / sum(seq_weights))
################################################################################
# Window Score
################################################################################
def window_score(scores, window_len, lam=.5):
	""" This function takes a list of scores and a length and transforms them
	so that each position is a weighted average of the surrounding positions.
	Positions with scores less than zero are not changed and are ignored in the
	calculation. Here window_len is interpreted to mean window_len residues on
	either side of the current residue. """
	w_scores = scores[:]
	for i in range(window_len, len(scores) - window_len):
		if scores[i] < 0:
			continue
		sum = 0.
		num_terms = 0.
		for j in range(i - window_len, i + window_len + 1):
			if i != j and scores[j] >= 0:
				num_terms += 1
				sum += scores[j]
		if num_terms > 0:
			w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]
	return w_scores


def calc_z_scores(scores, score_cutoff):
	"""Calculates the z-scores for a set of scores. Scores below
	score_cutoff are not included."""
	average = 0.
	std_dev = 0.
	z_scores = []
	num_scores = 0
	for s in scores:
		if s > score_cutoff:
			average += s
			num_scores += 1
	if num_scores != 0:
		average /= num_scores
	for s in scores:
		if s > score_cutoff:
			std_dev += ((s - average)**2) / num_scores
			std_dev = math.sqrt(std_dev)
	for s in scores:
		if s > score_cutoff and std_dev != 0:
			z_scores.append((s-average)/std_dev)
		else:
			z_scores.append(-1000.0)
	return z_scores
######################################################

# Calculate occupancy per residue
occupancy = calcMSAOccupancy(parseMSA(MSA))
writeArray(os.path.join(MIpath, "occupancy.txt"), occupancy)

# Calculate entropy per residue
entropy = calcShannonEntropy(parseMSA(MSA))
writeArray(os.path.join(MIpath, "entropy.txt"), entropy)

# Calculate conservation per residue
# Set default parameters
# AMI: READ msa alignment in fasta format
alignment = AlignIO.read(open(MSA), "fasta")

# TODO? -> set this parameters as options for the script. customizables by user
window_size = 3	# 0 = no window
win_lam = .5	# for window method linear combination
s_matrix_file = "scripts/blosum62.bla"
bg_distribution = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
scoring_function = js_divergence
use_seq_weights = True
background_name = 'blosum62'
gap_cutoff = .3
use_gap_penalty = 1

PSEUDOCOUNT = .0000001

# Calculate conservation per residues
# Read substitution matrix scores
s_matrix = read_scoring_matrix(s_matrix_file, aa_to_index)

# Read fasta alignment
names, alignment = read_fasta_alignment(MSA)

# Calculate sequence weights
seq_weights = calculate_sequence_weights(alignment, amino_acids)
if len(seq_weights) != len(alignment): seq_weights = [1.] * len(alignment)

# Calculate scores
scores = []
for i in range(len(alignment[0])):
	col = get_column(i, alignment)
	if len(col) == len(alignment) and gap_percentage(col) <= gap_cutoff:
		scores.append(scoring_function(col, s_matrix, bg_distribution, seq_weights, PSEUDOCOUNT, amino_acids, use_gap_penalty))
	else:
		scores.append(-.1)

# window_size scores
wscores = window_score(scores, window_size, win_lam)

# normalize_scores
zscores = calc_z_scores(scores, -999)

# print Conservation to file/stdout
outfile_name = 'coevolution/MI/conservation.txt'
outfile = open(outfile_name, 'w')
outfile.write("# %s -- %s - window_size: %d - window lambda: %.2f - background: %s - seq. weighting: %s - gap penalty: %d\n" % (MSA, scoring_function.__name__, window_size, win_lam, background_name, use_seq_weights, use_gap_penalty))
outfile.write("# align_column_number\tscore\tcolumn\n")

for i, score in enumerate(scores):	# Choose: scores, wscores or zscores
	outfile.write("%d\t%5f\t%s\n" % (i, score, "".join(get_column(i, alignment))))
outfile.close()

# List with columns from alignment - convert each column into string and add to list "columns"
l = len(alignment[0])	# length of sequences in alignment
columns = []
for i in range(0, l):
	columns.append([x[i] for x in alignment])

# Calculate Adjusted Mutual Information (AMI)
# Initializate empty array
AMI = np.empty((l,l), dtype="f")
AMI_jaccard = np.empty((l,l), dtype="f")

# Iterate over each pair of columns and calculate AMI
significant_interactions = []
for i in range(0, l):
	print(i)
	for j in range(0, l):
		ami_score = adjusted_mutual_info_score(columns[i], columns[j])
		AMI[i,j] = ami_score
		AMI_jaccard[i,j] = jaccard(columns[i], columns[j])
		# Selecte significant pairs of columns | AMI > 0.6 as default
		if abs(ami_score) > float(cut):
			if j >= i: continue
			significant_interactions.append([min(i,j)+1, max(i,j)+1, ami_score])

# Write matrixes into files
np.savetxt(os.path.join(MIpath, "AMI.txt"), AMI, fmt='%.4f', delimiter=',', newline='\n')
np.savetxt(os.path.join(MIpath, "AMIjaccard.txt"), AMI_jaccard, fmt='%.4f', delimiter=',', newline='\n')

# Write significant interaction into files
with open(MIpath+"/significant_interactions.txt", 'w') as fo:
	fo.write("\t".join(["residue1", "residue2", "AMI\n"]))
	for i in significant_interactions:
		i = [str(x) for x in i]
		fo.write("\t".join(i+["\n"]))

# Run Rscript to plot matrixes
log = MIpath+"/AMI.log"
script = os.path.join("scripts", "MutualInformation.R")
Rargs = " ".join(["Rscript", script, MSA, MIpath, str(cut), ">", log])
subprocess.call(Rargs, shell=True, stderr=open(log, 'wb'))