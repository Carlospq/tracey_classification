### Shared helpers for the SNARE/Habc classification pipeline ###
# Consolidates functions that were previously copy-pasted (with minor drift) across
# scripts/func.py, scripts/generateStatisticalHmmModels.py, scripts/domainPredict.py,
# results/hmmscan_to_table.py and initialSetSelection/setSelection.py.
import os
import subprocess

import numpy as np
import pandas as pd
from Bio import AlignIO
import seqlogo


def readFasta(fasta_file):
	# Returns a dictionary with the fasta file content {seq_name: sequence}
	fasta = {}
	with open(fasta_file) as f:
		for line in f:
			if line.startswith('>'):
				seq_name = line[1:].strip()
				fasta[seq_name] = ''
			else:
				fasta[seq_name] += line.strip()
	return fasta


def cleanFasta(fasta_file, replace=True):
	fasta = readFasta(fasta_file)
	new_fasta = {}
	# Remove empty sequences and clean dashes
	for seq in list(fasta.keys()):
		fasta[seq] = fasta[seq].replace("-", "")
		if not len(fasta[seq]): continue
		new_fasta[seq.replace(".", "_")] = fasta[seq].replace('-', '')

	if replace:
		writeFasta(new_fasta, fasta_file, replace=replace)

	return new_fasta


def writeFasta(fasta_dict, fasta_file, replace=False):
	if replace:
		fasta_file = fasta_file.replace(".fasta", "_clean.fasta")
	with open(fasta_file, 'w') as f:
		for seq in fasta_dict:
			f.write(">" + seq + "\n" + fasta_dict[seq] + "\n")
	if replace:
		os.rename(fasta_file, fasta_file.replace("_clean.fasta", ".fasta"))


def alignFasta(sequences_fasta, alignment_file):
	# Requires the `muscle` binary to be available on PATH
	cmd = 'muscle'
	bashCommand = "%s -super5 %s -output %s" % (cmd, sequences_fasta, alignment_file)
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()
	return output, error


def trimAlignment(alignment_file, trimmed_alignment_file, col_threshold=0.8, row_threshold=0.8):
	if type(alignment_file) is dict:
		alignment = alignment_file
	else:
		alignment = readFasta(alignment_file)
	alignment_trimm = {}
	alignment_trimm_ = {}
	elements = [[k, v] for k, v in alignment.items()]			# [[seq_id, seq in alignment], ... ]
	# iterate over columns in alignment
	for n in range(len(elements[0][1])):
		coln = [x[1][n] for x in elements]
		# Keep column if column has less gaps than col_threshold (default: 80%)
		if coln.count('-') / len(coln) < col_threshold:
			for name in alignment:
				alignment_trimm_[name] = alignment_trimm_.get(name, '') + alignment[name][n]
	# Iterate over rows in alignment_trimm_ and skip if row has more than row_threshold (default: 80%) of gaps
	for name in alignment_trimm_:
		if alignment_trimm_[name].count('-') / len(alignment_trimm_[name]) < row_threshold:
			alignment_trimm[name] = alignment_trimm_[name].replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('*', '-').replace('x', '-')
			length_seq = len(alignment_trimm[name])
	alignment_dims = [len(alignment_trimm), length_seq]
	# Write alignment into file
	with open(trimmed_alignment_file, 'w') as af:
		for name in alignment_trimm:
			af.write(">" + name + "\n" + alignment_trimm[name] + "\n")
	return alignment_trimm, alignment_dims


def searchHMMs(sequences_file, hmm_file, domtbloutPath):
	# Search sequences against hmm DB and store results in domtblout files
	cmd = "hmmsearch --domtblout %s --cpu 4 %s %s" % (domtbloutPath, hmm_file, sequences_file)
	out, error = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print(error)


def retriveEvalues(domtbloutFile, motifs):
	"""Best (lowest) raw e-value per sequence/motif, with a 0 -> 0.01 floor.

	Used by the training pipeline (generateStatisticalHmmModels.py, results/hmmscan_to_table.py).
	NOTE: scripts/domainPredict.py's inference path uses retriveEvaluesLog instead, a
	log-transformed variant. That divergence predates this cleanup and is preserved as-is
	rather than silently unified, since merging them would change either training or
	inference numerics -- flagged here for the author to double check.
	"""
	evalues = {}
	for line in open(domtbloutFile, 'r').readlines():
		if line.startswith('#'): continue
		line = line.split()
		line[12] = float(line[12])
		seq_name = line[0]
		motif_name = line[3]
		if not motif_name in motifs: continue
		if seq_name not in evalues:
			evalues[seq_name] = {}
			for m in motifs:
				evalues[seq_name][m] = 10e10
		if float(line[12]) < evalues[seq_name][motif_name]:
			if float(line[12]) == 0:
				line[12] = 0.01
			evalues[seq_name][motif_name] = float(line[12])
	return evalues


def retriveEvaluesLog(domtbloutFile, all_motifs):
	"""Best (lowest) log-transformed e-value per sequence/motif.

	Used by the inference/comparison path (domainPredict.py, HMMcomparison.py).
	See retriveEvalues for the raw/clamped variant used by the training pipeline.
	"""
	results = {}
	with open(domtbloutFile, 'r') as f:
		for line in f:
			if line.startswith("#"): continue
			line = line.strip().split()
			seq_name = line[0]
			hmm_name = line[3]
			evalue = np.log(float(line[12]))
			if hmm_name not in all_motifs: continue
			if seq_name not in results:
				results[seq_name] = {}
				for m in all_motifs:
					results[seq_name][m] = np.log(10e10)
			if evalue < results[seq_name][hmm_name]:
				results[seq_name][hmm_name] = evalue
	return results


def nestedDictValues(d, values=None, rank=0):
	if values is None:
		values = {}
	if not rank in values:
		values[rank] = []
	for k in d:
		if type(d[k]) is dict:
			values = nestedDictValues(d[k], values, rank+1 if rank < 3 else rank)
		if not k in values[rank]:
			values[rank].append(k)
	return values


def plot_web_logo(aln, pdf_file, title=''):

	def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
		alnRows = aln.get_alignment_length()
		compDict = {char: [0] * alnRows for char in characters}
		for record in aln:
			header = record.id
			seq = record.seq
			for aaPos in range(len(seq)):
				aa = seq[aaPos]
				if aa in characters:
					compDict[aa][aaPos] += 1
		return pd.DataFrame.from_dict(compDict)

	# Load alignment data
	t7_alignmentFile = aln
	t7_alignment = AlignIO.read(t7_alignmentFile, format="fasta")

	# Calculate amino acid frequency
	t7_alignmentSiteCompDF = alnSiteCompositionDF(t7_alignment)
	t7_alignmentSiteFreqDF = t7_alignmentSiteCompDF.div(t7_alignmentSiteCompDF.sum(axis=1), axis=0)

	# Create sequence logo
	t7_alignmentSiteFreqSeqLogo = seqlogo.Ppm(t7_alignmentSiteFreqDF, alphabet_type="AA")
	seqlogo.seqlogo(t7_alignmentSiteFreqSeqLogo, filename=pdf_file, format='pdf',
					title=title, size='xlarge', ic_scale=True,
					stacks_per_line=400, number_interval=50, number_fontsize=3,
					color_scheme='chemistry', stack_aspect_ratio=10, stack_margin=0.3)
