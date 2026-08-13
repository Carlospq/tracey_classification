### Imports ###
import os
import subprocess
import pickle

import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import train_test_split
from scripts.motifNames import *
from scripts.pipeline_utils import (
	readFasta, cleanFasta, writeFasta, alignFasta, trimAlignment,
	searchHMMs, retriveEvalues, nestedDictValues,
)
####################################################################################################
### Global variables ###
fastasPath = "fastas/"
alignmentsPath = "alignments/"
hmmsPath = "hmms/"
treesPath = "trees/"
regressionPath = "regressionModel/"
initialSetSelectionPath = "initialSetSelection/"

####################################################################################################
### Functions ###
def trainModel(model, groupsList, initialdata, modelName, initialColumn = "Qa", groupColumn = "Group"):
	print("Saving model for %s (%s)..." % (modelName, model))
	data = initialdata[initialdata[groupColumn].isin(groupsList)]
	motifsData_X = []
	motifsData_Y = []
	for index, row in data.iterrows():
		motifsData_X.append(np.log(row.loc[initialColumn:].values.flatten().tolist()))
		motifsData_Y.append(row.loc[groupColumn])
	#regModelsData[model][group + "Data_X"] = motifsData_X
	#regModelsData[model][group + "Data_Y"] = motifsData_Y
	regmodel = LogisticRegression(max_iter=10000)
	regmodel.fit(motifsData_X, motifsData_Y)
	# regModelsData[model][modelName + "Model"] = regmodel
	pickle.dump(regmodel, open('regressionModel/lm_%s_%s.sav' % (model, modelName), 'wb'))
	return [regmodel, motifsData_X, motifsData_Y]


def getTrainingData(initialdata, groupsList, initialColumn = "Qa", groupColumn = "Group"):
	data = initialdata[initialdata[groupColumn].isin(groupsList)]
	motifsData_X, motifsData_Y = [], []
	for index, row in data.iterrows():
		motifsData_X.append(np.log(row.loc[initialColumn:].values.flatten().tolist()))
		motifsData_Y.append(row.loc[groupColumn])
	return motifsData_X, motifsData_Y
####################################################################################################

### Get all domains in lists
domainSNARE = nestedDictValues(SNARETree)[0]
groupSNARE = nestedDictValues(SNARETree)[1]
subgroupSNARE = [hmm for rank in nestedDictValues(SNARETree) for hmm in nestedDictValues(SNARETree)[rank] if rank >= 2]
allSNARE = domainSNARE + groupSNARE + subgroupSNARE

### Generate fasta files for each main group with sequences from tracey (step done manually) ###

### Clean fasta files
print("# Cleaning fasta files")
for fastaFile in os.listdir(fastasPath):
	if not fastaFile.endswith(".fasta"): continue
	if "notSNARE" in fastaFile: continue
	if fastaFile.replace(".fasta", ".aln") in os.listdir(alignmentsPath): continue
	print("Cleaning %s" % fastaFile)
	cleanFasta(fastasPath + fastaFile)

### Align full SNARE sequences (step done on an HPC cluster)
print("# Aligning sequences")
for fastaFile in os.listdir(fastasPath):
	if not fastaFile.endswith(".fasta"): continue
	if "notSNARE" in fastaFile: continue
	if not os.path.isfile(alignmentsPath + fastaFile.replace(".fasta", ".aln")):
		print("Aligning %s" % (fastaFile))
		alignFasta(fastasPath + fastaFile, alignmentsPath + fastaFile.replace(".fasta", ".aln"))

### Trim alignments
print("# Trimming alignments")
for alignmentFile in os.listdir(alignmentsPath):
	if alignmentFile.endswith("trimmed.aln"): continue
	trimmedAln = alignmentsPath + alignmentFile.replace(".aln", ".trimmed.aln")
	if not os.path.isfile(trimmedAln):
		print("Trimming alignment %s" % (alignmentFile))
		trimAlignment(alignmentsPath + alignmentFile, trimmedAln)

### Build Trees  (recommended: run on an HPC cluster) from these alignments and use this information to generate subgroups by selecting initial sequences for each group
print("# Building trees")
for trimmedAln in os.listdir(alignmentsPath):
	if not trimmedAln.endswith("trimmed.aln"): continue
	treeFile = treesPath + trimmedAln + ".treefile"
	if not os.path.isfile(treeFile):
		print("Building tree for %s" % (trimmedAln))
		cmd = "iqtree2 -s %s -B 1000 --prefix %s" % (alignmentsPath + trimmedAln, treesPath + trimmedAln)
		process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, error = process.communicate()

### Plot trees PDF - For some reason running the command from python may not work. Run directly in Rstudio if necessary
print("# Plotting trees")
cmd = "Rscript scripts/PlotTree.R"
process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, error = process.communicate()


### Initial sequence set selection for subgroups and main domains (initialSetSelection/setSelection.py)
# Seed sequence lists for each subgroup (initialSetSelection/initialSets/*.initialSet.txt) are still
# curated manually from the trees above; this just runs the selection/refinement itself per fasta.
# new fastas, alignments and hmms are generated in the iterations folder
print("# Selecting initial sequence sets")
from initialSetSelection.setSelection import runInitialSetSelection
for fastaFile in os.listdir(fastasPath):
	if not fastaFile.endswith(".fasta"): continue
	prefix = fastaFile.replace(".fasta", "")
	if prefix not in domainSNARE and prefix not in groupSNARE: continue
	print("Selecting initial set for %s" % (prefix))
	runInitialSetSelection(prefix)

### Copy fastas/alignments/hmms from iterations to the corresponding folders
print("# Copying files")
for fastaFile in os.listdir(initialSetSelectionPath + "iterations/"):
	if fastaFile.endswith(".fasta"):
		cmd = "cp %s %s" % (initialSetSelectionPath + "iterations/" + fastaFile, fastasPath)
		## Generate idsForSubgroups.txt files from fastas in iterations
		ids = [x.split("/")[0][1:] for x in open(initialSetSelectionPath + "iterations/" + fastaFile, 'r').readlines() if x.startswith(">")]
		with open("idsForSubgroups/%s"%(fastaFile.replace("fasta", "txt")), "w") as f:
			for i in ids:
				f.write(i)
	elif fastaFile.endswith(".aln"):
		cmd = "cp %s %s" % (initialSetSelectionPath + "iterations/" + fastaFile, alignmentsPath)
	elif fastaFile.endswith(".hmm"):
		cmd = "cp %s %s" % (initialSetSelectionPath + "iterations/" + fastaFile, hmmsPath)
	else:
		continue
	print("Copying %s" % (fastaFile))
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()


### Generate HMMs for main groups (Qa, Qb, ...) and generate HMM db concatenating all HMMs
# No domain-wide fasta is ever supplied directly (a domain-wide tree/alignment is too costly to
# build locally, hence the group/subgroup split above), so this redoes clean/align/trim/hmmbuild
# just for the main domains -- skipping tree-building -- using a fasta merged from every
# fastas/{domain}.*.fasta already on disk (group- and subgroup-level, including the ones the
# "Copying files" step above just produced).
# Some domains were exported under a legacy raw filename before the group/subgroup naming
# convention was settled (SNAP's combined domain fasta is fastas/SNAPbc.fasta, not fastas/SNAP.fasta).
domainFastaAliases = {"SNAP": "SNAPbc"}

print("# Generating main HMMs")
for domain in domainSNARE:
	hmmFile = hmmsPath + domain + ".hmm"
	if os.path.isfile(hmmFile): continue

	domainFastaFile = fastasPath + domain + ".fasta"
	if not os.path.isfile(domainFastaFile):
		domainSeqs = {}
		alias = domainFastaAliases.get(domain)
		for fastaFile in os.listdir(fastasPath):
			if not fastaFile.endswith(".fasta"): continue
			name = fastaFile[:-len(".fasta")]
			if not (name.startswith(domain + ".") or name == alias): continue
			domainSeqs.update(readFasta(fastasPath + fastaFile))
		if not domainSeqs:
			print("No fastas found for %s, skipping" % (domain))
			continue
		print("Merging fastas for %s (%d sequences)" % (domain, len(domainSeqs)))
		writeFasta(domainSeqs, domainFastaFile)

	alignment = alignmentsPath + domain + ".aln"
	trimmedAln = alignmentsPath + domain + ".trimmed.aln"
	if not os.path.isfile(alignment):
		print("Aligning %s" % (domainFastaFile))
		alignFasta(domainFastaFile, alignment)
	if not os.path.isfile(trimmedAln):
		print("Trimming alignment for %s" % (domain))
		trimAlignment(alignment, trimmedAln)

	print("Generating HMM for %s" % (domain))
	cmd = "hmmbuild -n %s %s %s" % (domain, hmmFile, trimmedAln)
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()

print("# Generating HMM db")
dbName = "SNAREDb"
hmmDbFile = hmmsPath + "SNAREDb.hmm"
hmmFiles = [f'%s%s'%(hmmsPath, x) for x in os.listdir(hmmsPath) if x.replace(".hmm", "") in allSNARE]
with open(hmmDbFile, 'w') as outfile:
	for hname in hmmFiles:
		with open(hname) as infile:
			for line in infile:
				outfile.write(line)
for fname in [x for x in os.listdir(hmmsPath) if x.startswith("SNAREDb.hmm.")]:
	os.remove(hmmsPath+fname)
cmd = "hmmpress %s" % (hmmDbFile)
process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, error = process.communicate()


### Run hmmsearch for notSNARE sequences
print("Generating notSNARE (negative dataset) sequences")
hmmDbFile = "hmms/SNAREDb.hmm"
notSNAREfasta = "fastas/notSNARE.fasta"
domtblout = regressionPath + "/notSNARE.domtblout"
if not os.path.isfile(domtblout):
	searchHMMs(notSNAREfasta, hmmDbFile, domtblout)
else:
	print("%s already exists" % (domtblout))


### Search for HMMs using all sequences
for fasta in os.listdir('fastas/'):
	# name = ".".join(fasta.split('.')[1:-1])
	name = fasta.replace(".fasta", "")
	if not name: continue
	print(name)
	domtblout = regressionPath + name + '.domtblout'
	cmd = "hmmsearch --domtblout %s %s %s" % (domtblout, hmmDbFile, 'fastas/' + fasta)
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()
	filteredLines = [x for x in open(domtblout, 'r').readlines() if not x.startswith("#")]
	# Get only best hit for each motif and sequence
	bestHits = {}
	for line in filteredLines:
		seqId = line.split()[0]
		motif = line.split()[3]
		if not seqId in bestHits:
			bestHits[seqId] = {}
		if not motif in bestHits[seqId]:
			bestHits[seqId][motif] = line
		else:
			if float(line.split()[11]) < float(bestHits[seqId][motif].split()[11]):
				bestHits[seqId][motif] = line
	# Write filtered domtblout
	with open(domtblout, 'w') as fo:
		#for line in filteredLines:
		#	nullhandler = fo.write(line)
		for b in bestHits:
			for m in bestHits[b]:
				nullhandler = fo.write(bestHits[b][m])


### Parse domtblout files and store results in a dictionary ###
print("Parsing domtblout files")
domtbloutFiles = [os.path.join(regressionPath, f) for f in os.listdir(regressionPath) if f.endswith(".domtblout")]
results = {}

# all_motifs =  [hmm for rank in nestedDictValues(SNARETree) for hmm in nestedDictValues(SNARETree)[rank]]
all_motifs = allSNARE

for domtblout in domtbloutFiles:
	print(domtblout)
	sequencesGroup = domtblout.split("/")[1].replace(".domtblout", "")
	results[sequencesGroup] = {}
	with open(domtblout, 'r') as f:
		for line in f:
			if line.startswith("#"): continue
			line = line.strip().split()
			seq_name = line[0]
			hmm_name = line[3]
			hmm_start = int(line[15])
			hmm_end = int(line[16])
			seq_start = int(line[17])
			seq_end = int(line[18])
			length = seq_end - seq_start
			evalue = float(line[12])
			score = float(line[7])
			if seq_name not in results[sequencesGroup]:
				results[sequencesGroup][seq_name] = {}
				for m in all_motifs:
					results[sequencesGroup][seq_name][m] = {"eval": 10e10, "len": 0}
			if not hmm_name in results[sequencesGroup][seq_name]:
				results[sequencesGroup][seq_name][hmm_name]["eval"] = 10e10
			if evalue < results[sequencesGroup][seq_name][hmm_name]["eval"]:
				results[sequencesGroup][seq_name][hmm_name]["eval"] = evalue
				results[sequencesGroup][seq_name][hmm_name]["len"]  = length


### Add set of negative sequences with random evals
n = 10000
for i in range(n):
	results["notSNARE"]["randomSeq_" + str(i)] = {}
	for m in all_motifs:
		results["notSNARE"]["randomSeq_" + str(i)][m] = {"eval": 10e10, "len": 0}
		if random.uniform(0, 1) > 0.99:
			eval = random.uniform(0, 3)
			length = random.uniform(5, 15)
			results["notSNARE"]["randomSeq_" + str(i)][m] = {"eval": eval, "len": length}


### Prepare classification pd ###
data = []
for group in results:
	for seqId in results[group]:
		#data.append([seqId, group] + [results[group][seqId][motif]["eval"] for motif in all_motifs] + [results[group][seqId][motif]["len"] for motif in all_motifs])
		data.append([seqId, group] + [results[group][seqId][motif]["eval"] for motif in all_motifs])

datapd =  pd.DataFrame(data, columns = ['ID', 'Group'] + all_motifs)


### Including length of hmm match in case is needed for a better classification
# datapd = pd.DataFrame(data, columns = ['ID', 'Group'] + all_motifs + ["len"+m for m in all_motifs])
# ...


### Generate regression models
redo = True

# Split data for each of the lr models -- Habc and SNARE Domain, Group and Subgroup models
print("Gather data for model testing...")
domainGroups = {"SNARE": {"mainMotifs": domainSNARE,
						  "groupMotifs": groupSNARE,
						  "subMotifs": subgroupSNARE},
				"Habc": {}
				}
lrModels = {}

if redo:
	for model in ["SNARE"]: #["Habc", "SNARE"]:

		lrModels[model] = {}
		tree = SNARETree if model == "SNARE" else "" #HabcTree

		# Models considering all motifs at once -- mainModel (Qa, Qb, Qc, R, SNAP), groupModel (Qa.I, Qa.II, ...) and subgroupModel (Qa.I.Syx18, Qa.I.Syx5, ...)
		for group in ["mainMotifs", "groupMotifs", "subMotifs"]:

			groupList = domainGroups[model][group]
			lrModels[model][group] =  {'model': '', 'dataX': '', 'dataY': ''}
			lrModels[model][group]["dataX"], lrModels[model][group]["dataY"] = getTrainingData(datapd, groupList)

		# Models grouped by main motifs and group motifs -- (Qa: Qa.I, Qa.II, Qa.III, Qa.IV; ...)
		for mainMotif in domainGroups[model]["mainMotifs"]:

			groupsList = [x for x in domainGroups[model]["groupMotifs"] if x.startswith(mainMotif)]
			lrModels[model][mainMotif] = {'model': '', 'dataX': '', 'dataY': ''}
			lrModels[model][mainMotif]["dataX"], lrModels[model][mainMotif]["dataY"] = getTrainingData(datapd, groupsList)

			for group in groupsList:

				if not tree[mainMotif][group]: continue
				rankedMotifs = nestedDictValues(tree[mainMotif][group])

				# Models for each group and rank motif considering only its submotifs (Qa.I: Qa.I.Syx18, Qa.I.Ufe1; Qa.I.Syx18: Qa.I.Syx18.Metazoa, Qa.I.Syx18.SAR, Qa.I.Syx18.Viridiplantae; ...)
				for rank in nestedDictValues(tree[mainMotif][group]):

					subgroupsList = rankedMotifs[rank]
					lrModels[model]["%s_rank%s"%(group, rank)] = {'model': '', 'dataX': '', 'dataY': ''}
					lrModels[model]["%s_rank%s"%(group, rank)]["dataX"], lrModels[model]["%s_rank%s"%(group, rank)]["dataY"] = getTrainingData(datapd, subgroupsList)


### TEST MODELS ###
if redo:

	n = 10
	lrScores = {}
	for domain in lrModels:
		lrScores[domain] = {}
		for group in lrModels[domain]:
			lrScores[domain][group] = []

	for domain in lrModels:
		for model in lrModels[domain]:
			print("Testing model %s %s..." % (domain, model))

			modelScore = 0.0
			modelData = [lrModels[domain][model]['dataX'], lrModels[domain][model]['dataY']]

			for i in range(n):
				# Prepare training and testing datasets
				x_train, x_test, y_train, y_test = train_test_split(modelData[0], modelData[1], test_size=0.2)
				lm = LogisticRegression(solver='liblinear', C=0.05, multi_class='ovr')
				fitted = lm.fit(x_train, y_train)
				y_pred = fitted.predict(x_test)
				# Scoring model
				score = fitted.score(x_test, y_test)
				lrScores[domain][model].append(score)
				print("\t\tIteration (%d/%d): "%(i+1, n) + "\tTest accuracy: %s" % (fitted.score(x_test, y_test)))
				if score > modelScore:
					modelScore = score
					lrModels[domain][model]['model'] = fitted
					lrModels[domain][model]['score'] = score
					plotData = [fitted, x_test, y_test, y_pred]

			### Report and confusion matrix from the best iteration ###
			fitted = plotData[0]
			x_test = plotData[1]
			y_test = plotData[2]
			y_pred = plotData[3]
			with open(regressionPath + "/report_%s_%s.txt" % (domain, model), 'w') as fo:
				fo.write("Best model score: %s\n" % (lrModels[domain][model]['score']))
				fo.write(classification_report(y_test, y_pred))
				titles_options = [
					("Confusion matrix (%s, %s), without normalization" % (domain, model), None),
					("Normalized confusion matrix (%s, %s)" % (domain, model), "true"),
				]
				for title, normalize in titles_options:
					disp = ConfusionMatrixDisplay.from_estimator(
						fitted,
						x_test,
						y_test,
						cmap=plt.cm.Blues,
						normalize=normalize,
					)
					disp.ax_.set_title(title)
					fig, ax = plt.subplots(figsize=(15, 15))
					disp.plot(ax=ax)
					plt.xticks(rotation=90)
					plt.savefig(regressionPath + "/confusionMatrix_%s_%s_%s.png" % (domain, model, normalize))

	### Save best iteration models ###
	for domain in lrModels:
		for model in lrModels[domain]:
			pickle.dump(lrModels[domain][model]["model"], open(regressionPath+'/lm_%s_%s.sav'%(domain, model), 'wb'))


### Report average scores ###
fmt = '{:<35}{:<20}'
print("Average scores for %s iterations:\t"%(n))
for domain in lrScores:
	for model in lrScores[domain]:
		print(fmt.format(  "Model %s (%s):"%(model, domain), str(round( np.mean(lrScores[domain][model]), 4))))













