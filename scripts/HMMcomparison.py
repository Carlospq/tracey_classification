# Comparison between old and new SNARE HMMs
import os, tempfile
import numpy as np

# searchHMM, predictFromFasta, snareMotifPrediction, formatPrediction, domainPredict, nestedDictValues, retriveEvalues, seqToFasta
from scripts.domainPredict import *


#################################################################
def retriveScores(domtbloutFile, all_motifs):
	""""
	Retrieves the e-values from a domtblout file and returns them in a dictionary.
	Args:
		domtbloutFile (str): Path to the domtblout file.
		all_motifs (list): List of all motifs to be considered.
		column (int): Column number in the domtblout file that contains the e-values (default is 12, for p-val. Use 13 for bit-score. Use 14 for bias).
	"""""
	results = {}
	with open(domtbloutFile, 'r') as f:
		for line in f:

			if line.startswith("#"): continue

			# Read values
			line = line.strip().split()
			seq_name = line[0]
			hmm_name = line[3]
			cEval = -np.log(float(line[11]))
			iEval = -np.log(float(line[12]))
			score = float(line[13])
			bias = float(line[14])
			length = int(line[20]) - int(line[19]) + 1

			# Initialize results for this sequence if not already done
			if hmm_name not in all_motifs: continue
			if seq_name not in results:
				results[seq_name] = {}

			for m in all_motifs:

				if not m in results[seq_name]:
					results[seq_name][m] = {}
					results[seq_name][m]['cEval'] = -np.log(10e10)
					results[seq_name][m]['iEval'] = -np.log(10e10)
					results[seq_name][m]['score'] = 0
					results[seq_name][m]['bias'] = 0
					results[seq_name][m]['length'] = 0

			# Update results for this motif if the values are lower than the current one
			if iEval > results[seq_name][hmm_name]['iEval'] or cEval > results[seq_name][hmm_name]['cEval']:
				results[seq_name][hmm_name]['cEval'] = cEval
				results[seq_name][hmm_name]['iEval'] = iEval
				results[seq_name][hmm_name]['score'] = score
				results[seq_name][hmm_name]['bias'] = bias
				results[seq_name][hmm_name]['length'] = length

	return results

#################################################################
fastaDir = 'fastas'
newHmmsDb = 'hmms/SNAREDb.hmm'
traceyHmmsDb = 'traceyHmms/traceyHmmDb.hmm'
motifNames = nestedDictValues(SNARETree)[0] + nestedDictValues(SNARETree)[1] + nestedDictValues(SNARETree)[2] + nestedDictValues(SNARETree)[3]

results = {}
n = 1
for db in [newHmmsDb, traceyHmmsDb]:

	print(db)
	results[db] = {}
	for fastaFile in os.listdir(fastaDir):

		results[db][fastaFile] = {}
		if fastaFile == "notSNARE.fasta": continue

		print('%s/%s: %s' % (n, len(os.listdir(fastaDir)), fastaFile))
		fastaPath = os.path.join(fastaDir, fastaFile)
		new_file, domtblout = tempfile.mkstemp(dir="tmp")
		searchHMMs(fastaPath, db, domtblout)
		values = retriveScores(domtblout, motifNames)
		os.close(new_file)
		os.remove(domtblout)

		for seqId in values:

			if not seqId in results[db][fastaFile]:
				results[db][fastaFile][seqId] = []

			results[db][fastaFile][seqId] += [values[seqId][m]['iEval'] for m in motifNames]
			results[db][fastaFile][seqId] += [values[seqId][m]['cEval'] for m in motifNames]
			results[db][fastaFile][seqId] += [values[seqId][m]['score'] for m in motifNames]
			results[db][fastaFile][seqId] += [values[seqId][m]['bias'] for m in motifNames]
			results[db][fastaFile][seqId] += [values[seqId][m]['length'] for m in motifNames]

		n += 1
	n = 1



# Wide format
with open('tmp/results_wide.csv', 'w') as of:

	headers = ['id', 'fasta'] + ['newiEval_%s' % x for x in motifNames] + ['newcEval_%s' % x for x in motifNames] + [
								 'newScore_%s' % x for x in motifNames] + ['newBias_%s' % x for x in motifNames] + [
								 'newLength_%s' % x for x in motifNames] + ['oldiEval_%s' % x for x in motifNames] + [
								 'oldcEval_%s' % x for x in motifNames] + ['oldScore_%s' % x for x in motifNames] + [
								 'oldBias_%s' % x for x in motifNames] + ['oldLength_%s' % x for x in motifNames]

	of.write(",".join(headers) + "\n")

	for fastaFile in os.listdir(fastaDir):

		if fastaFile == "notSNARE.fasta": continue

		# Data
		newData = results[newHmmsDb][fastaFile]
		oldData = results[traceyHmmsDb][fastaFile]

		# Combine results
		allIds = list(set(newData.keys()).union(set(oldData.keys())))

		# Write results for each sequence ID
		for seqId in allIds:
			# Initializate data if not present ;iEval + cEval                             ;Score + Bias               ;Length
			seqIdNewData = newData.get(seqId, [-25.3284360229345] * len(motifNames) * 2 + [0] * len(motifNames) * 2 + [0] * len(motifNames))
			seqIdOldData = oldData.get(seqId, [-25.3284360229345] * len(motifNames) * 2 + [0] * len(motifNames) * 2 + [0] * len(motifNames))
			fastaFileName = fastaFile.replace(".fasta", "")
			of.write(",".join([seqId, fastaFileName] + [str(x) for x in seqIdNewData] + [str(x) for x in seqIdOldData]) + "\n")


# Long format
with open('tmp/results.csv', 'w') as of:

	headers = ['id', 'fasta', 'motif', 'newiEval', 'oldiEval', 'newcEval', 'oldcEval', 'newScore', 'oldScore', 'newBias', 'oldBias', 'newLength', 'oldLength']
	of.write(",".join(headers) + "\n")

	for fastaFile in os.listdir(fastaDir):
		if fastaFile == "notSNARE.fasta": continue

		# Data
		newData = results[newHmmsDb][fastaFile]
		oldData = results[traceyHmmsDb][fastaFile]

		# Combine results
		allIds = list(set(newData.keys()).union(set(oldData.keys())))

		# Added length
		mLength = len(motifNames)

		for seqId in allIds:
			seqIdNewData = newData.get(seqId, [25.3284360229345] * len(motifNames) * 2 + [0] * len(motifNames) * 2 + [0] * len(motifNames))
			seqIdOldData = oldData.get(seqId, [25.3284360229345] * len(motifNames) * 2 + [0] * len(motifNames) * 2 + [0] * len(motifNames))
			for i in range(len(motifNames)):
				fastaFileName = fastaFile.replace(".fasta", "")
				of.write(",".join([seqId, fastaFileName, motifNames[i], str(seqIdNewData[i]), str(seqIdOldData[i]),
								   										str(seqIdNewData[i + mLength]), str(seqIdOldData[i + mLength]),
								   										str(seqIdNewData[i + mLength * 2]), str(seqIdOldData[i + mLength * 2]),
								   										str(seqIdNewData[i + mLength * 3]), str(seqIdOldData[i + mLength * 3]),
								   										str(seqIdNewData[i + mLength * 4]), str(seqIdOldData[i + mLength * 4])
								   ]) + "\n")
