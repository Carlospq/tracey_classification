import os, sys, subprocess
import numpy as np
from sklearn.mixture import GaussianMixture
from statistics import NormalDist
import matplotlib.pyplot as plt
from scipy.stats import norm

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from scripts.pipeline_utils import alignFasta, trimAlignment

os.chdir("initialSetSelection")

####################################################################################################
### FUNCTIONS
def readFasta(fastaFile, cleanHeaders=True):
	sequences = {}
	with open(fastaFile, "r") as f:
		for line in f:
			if line.startswith(">"):
				line = line.strip().replace(">", "").replace(".", "_")
				if cleanHeaders:
					line = line.split("/")[0]
				header = line
				sequences[header] = ""
			else:
				sequences[header] += line.strip()
	return sequences

def writeFasta(fasta_dict, fasta_file, replace=False):
	if replace:
		fasta_file = fasta_file.replace(".fasta", "_clean.fasta")
	with open(fasta_file, 'w') as f:
		for seq in fasta_dict:
			f.write(">" + seq + "\n" + fasta_dict[seq] + "\n")
	if replace:
		os.rename(fasta_file, fasta_file.replace("_clean.fasta", ".fasta"))

def cleanFasta(fasta_file, replace=True, cleanHeaders=False):

	new_fasta = {}
	fasta = readFasta(fasta_file, cleanHeaders=cleanHeaders)
	averageLength = sum([len(fasta[x]) for x in fasta])/ len(fasta)

	# Remove empty sequences and clean dashes
	for seq in list(fasta.keys()):

		fasta[seq] = fasta[seq].replace("-", "")
		if not len(fasta[seq]) or len(fasta[seq]) < averageLength/2: continue
		new_fasta[seq.replace(".", "_")] = fasta[seq].replace('-', '')

	if replace:
		writeFasta(new_fasta, fasta_file, replace=replace)

	return new_fasta

####################################################################################################

group = "Qb.III"

# Get sequences for initialSet using Seq shortname from core sequences
subgroupsList = list(set([x.replace(group+".", "").replace(".initialSet.txt", "") for x in os.listdir("initialSets") if x.startswith(group+".") and x.endswith(".initialSet.txt")]))

for subgroup in subgroupsList:

	# Get initialSet shortnames
	initialSetShortnames = []
	with open(f"initialSets/{group}.{subgroup}.initialSet.txt", "r") as f:
		for line in f:
			initialSetShortnames.append(line.strip())
		print(subgroup, len(initialSetShortnames))
	# Read fasta file with all seuqences
	fastaFile = f'../fastas/{group}.fasta'
	sequences = readFasta(fastaFile, cleanHeaders=False)

	# Write intialSet sequences to file
	initialFastaFile = f'initialSets/{group}.{subgroup}.initialSet.fasta'
	with open(initialFastaFile, "w") as f:
		for shortname in initialSetShortnames:
			header = f'{shortname}'.replace(".", "_")
			if header not in sequences:
				print(header, "not found")
				continue
			else:
				sequence = sequences[header]
			f.write(f'>{header}\n{sequence}\n')

# Collect all initial sequences
allInitialSeqs = []
for subgroup in subgroupsList:
	allInitialSeqs += [x for x in readFasta(f'initialSets/{group}.{subgroup}.initialSet.fasta')]

# Iterate over selection steps to form a strict group of sequences for each subgroup
for subgroup in subgroupsList:

	print(f'Processing {group}.{subgroup}...')
	subgroup = f'{group}.{subgroup}'

	# empty files
	for f in [f'iterations/{subgroup}.data.csv', f'iterations/{subgroup}.fasta', f'iterations/{subgroup}.fasta']:
		if os.path.isfile(f):
			with open(f, 'w') as f: pass

	# Initialize variables
	maxIterations = 11
	earlyExitMax = 30
	iteration = 0
	earlyExit = 0
	mean_i0 = 0
	median_i0 = 0
	meanTreshold = 0.15
	thresholdDiff = 0
	incrementThresholdDiff = 0.075

	# Number of groups in the distribution
	groups=2

	while iteration < maxIterations:

		mainFastaFile = f'../fastas/{group}.fasta'
		iterationFastaFile = f'iterations/{subgroup}.fasta'
		alignment = f'iterations/{subgroup}.aln'
		alignmentTrimmed = f'iterations/{subgroup}.trimmed.aln'
		hmm = f'iterations/{subgroup}.hmm'
		results = f'iterations/{subgroup}.results'

		if iteration == 0:
			with open(iterationFastaFile, 'w') as f:
				for line in open(f'initialSets/{subgroup}.initialSet.fasta', 'r'):
					f.write(line)

		iterationSeqs = readFasta(iterationFastaFile)
		mainSeqs = readFasta(mainFastaFile)

		# 1. Align sequences
		alignmentOut, alignmentErr = alignFasta(iterationFastaFile, alignment)

		# 2. Trim alignment
		alnTrimmed, alignmentDims = trimAlignment(alignment, alignmentTrimmed, col_threshold=0.8, row_threshold=0.8)

		# 3. Generate new HMM from
		cmd = f"hmmbuild -n {subgroup} iterations/{subgroup}.hmm {alignmentTrimmed}"
		process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, error = process.communicate()

		# 4. Search for sequences in HMM
		cmd = f"hmmsearch --tblout {results} {hmm} {mainFastaFile}"
		process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, error = process.communicate()

		# 5. Get sequences belonging to the correct distribution
		resultsDict = {}

		# 5.1 Read results file
		for line in open(results, 'r'):
			if line.startswith("#") or not line: continue
			line = line.strip().split()
			if not line[0].split("/")[0] in resultsDict:
			#if not line[0] in resultsDict:
				resultsDict[line[0].split("/")[0]] = float(line[7])
				#resultsDict[line[0]] = float(line[7])
			else:
				if float(line[7]) < resultsDict[line[0].split("/")[0]]:
				#if float(line[7]) < resultsDict[line[0]]:
					resultsDict[line[0].split("/")[0]] = float(line[7])
					#resultsDict[line[0]] = float(line[7])

		# 5.2 Fit a Gaussian Mixture Model to the distribution
		multimodal_dist = np.array([-np.log2(x) for x in resultsDict.values()])
		gmm = GaussianMixture(n_components=groups)
		gmm.fit(multimodal_dist.reshape(-1, 1))

		means = gmm.means_
		weights = gmm.weights_

		# Convert covariance into Standard Deviation
		standard_deviations = gmm.covariances_ ** 0.5

		# 6. Select sequences with higher probability for correct distribution
		distMaxIndex = np.argmax(means)
		threshold = means[distMaxIndex][0] - standard_deviations[distMaxIndex][0] - thresholdDiff
		selectedSeqs = [k for k, v in resultsDict.items() if -np.log2(v) > threshold and k not in allInitialSeqs and not k in iterationSeqs] + [k for k in iterationSeqs]

		if len(iterationSeqs) >= len(selectedSeqs):
			thresholdDiff += threshold*incrementThresholdDiff if thresholdDiff == 0 else thresholdDiff*incrementThresholdDiff
			earlyExit += 1
			print(f'Early exit: {earlyExit} - Threshold: {threshold} - Selected: {len(selectedSeqs)} - Initial: {len(iterationSeqs)}')
			if earlyExit >= earlyExitMax: break
			continue
		else:
			thresholdDiff = 0
			earlyExit = 0

		with open(f'iterations/{subgroup}.data.csv', 'a') as f:
			f.write(f'{iteration} {threshold} {len(selectedSeqs)} {len(iterationSeqs)}\n')

		if earlyExit >= earlyExitMax: break

		# 7. Write selected sequences to file
		selectedFastaFile = f'iterations/{subgroup}.fasta'
		with open(selectedFastaFile, 'a') as f:
			for seq in selectedSeqs:
				if seq not in iterationSeqs:
					if not seq in mainSeqs: continue
					f.write(f'>{seq}\n{mainSeqs[seq]}\n')

		# 8. Print iteration information
		print(f'Iteration: {iteration} - Threshold: {threshold} - Selected: {len(selectedSeqs)} - Initial: {len(iterationSeqs)}')

		iteration += 1

# 9. Plot the distributions
multimodal_dist = np.array([-np.log2(x) for x in resultsDict.values()])
gmm = GaussianMixture(n_components=3)
gmm.fit(multimodal_dist.reshape(-1, 1))
means = gmm.means_
weights = gmm.weights_
standard_deviations = gmm.covariances_ ** 0.5
fig, axes = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=(6.4, 7))
axes[0].hist(multimodal_dist, bins=50, alpha=0.5)
x = np.linspace(min(multimodal_dist), max(multimodal_dist), 100)
for mean, std, weight in zip(means, standard_deviations, weights):
	pdf = weight * norm.pdf(x, mean, std)
	plt.plot(x.reshape(-1, 1), pdf.reshape(-1, 1), alpha=0.5)
plt.axvline(x=threshold, color='red', label='axvline - full height')
plt.show()






qaseqs = []
for f in os.listdir('fastas'):
	if not f.startswith("Qa"): continue
	with open(f'fastas/{f}', 'r') as fi:
		for line in fi:
			if line.startswith(">"):
				seqname = line.strip().split("/")[0]
				if not seqname in qaseqs:
					qaseqs.append(seqname)


main_domains = ["Qa", "Qb", "Qc", "R", "SNAPbc"]
for domain in main_domains:

	print( domain )
	# Read fasta file with all seuqences
	fastaFile = f'../fastas/{domain}.fasta'
	sequences = readFasta(fastaFile, cleanHeaders=False)

	# Write intialSet sequences to file
	initialFastaFile = f'initialSets/{domain}.initialSet.fasta'
	with open(initialFastaFile, "w") as f:
		for shortname in sequences:
			header = f'{shortname}'.replace(".", "_")
			sequence = sequences[shortname]
			f.write(f'>{header}\n{sequence}\n')

	cleanFasta(initialFastaFile, replace=True)

	alignment = f'iterations/{domain}.aln'
	alignmentTrimmed = f'iterations/{domain}.trimmed.aln'
	hmm = f'iterations/{domain}.hmm'
	results = f'iterations/{domain}.results'
	mainFastaFile = f'../fastas/{domain}.fasta'


	# 1. Align sequences
	alignmentOut, alignmentErr = alignFasta(initialFastaFile, alignment)

	# 2. Trim alignment
	alnTrimmed, alignmentDims = trimAlignment(alignment, alignmentTrimmed, col_threshold=0.8, row_threshold=0.8)

	# 3. Generate new HMM from
	cmd = f"hmmbuild -n {domain} {hmm} {alignmentTrimmed}"
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()

	# 4. Search for sequences in HMM
	cmd = f"hmmsearch --tblout {results} {hmm} {mainFastaFile}"
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()








#####
alns = [aln for aln in os.listdir('../alignments') if aln.startswith("Hb") and aln.endswith("trimmed.aln")]
for aln in alns:
	print(aln)
	alignment_file = f'../alignments/{aln}'
	clean_aln = cleanFasta(alignment_file, replace=False, cleanHeaders=True)
	writeFasta(clean_aln, f'../fastas/{aln.replace("trimmed.aln", "fasta")}', replace=False)