####	IMPORTS		####
import subprocess
import os, math

from scripts.pipeline_utils import readFasta, alignFasta, trimAlignment, plot_web_logo

####	FUNCTIONS	####
def filterHMMoutput(tblout, threshold=0.01):
	hits = []
	lines = str(tblout).split("\\n")
	for line in lines[14:]:
		if line.startswith("#"): continue
		if line == '  ------ inclusion threshold ------': break
		evalue = float(line.split()[3])
		if evalue <= threshold:
			hits.append(line.split()[-1])
	return hits


def sampleAlignment(alignment_file, start, end, sample_alignment_file=False, trimm=False, col_threshold=0.8, row_threshold=0.8):
	# Returns a dictionary with the fasta file content
	# {seq_name: sequence}
	alignment = readFasta(alignment_file)
	alignment_sample = {}
	for name in alignment:
		alignment_sample[name] = alignment[name][start - 1:end]
	# Write alignment into file
	if sample_alignment_file:
		if trimm:
			alignment_sample, alignment_dims = trimAlignment(alignment_sample, sample_alignment_file, col_threshold=col_threshold, row_threshold=row_threshold)
		else:
			with open(sample_alignment_file, 'w') as af:
				for name in alignment_sample:
					af.write(">" + name + "\n" + alignment_sample[name] + "\n")
	return alignment_sample


def buildHMM(alignment_file, hmm_name, hmm_file):
	cmd = "hmmbuild -n %s --amino %s %s" % (hmm_name, hmm_file, alignment_file)
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()


def calculate_sequence_weights(trimmed_alignment, amino_acids, aa_to_index):
	""" Calculate the sequence weights using the Henikoff '94 method for the given msa. """
	seq_weights = {}
	for seq_id in trimmed_alignment:
		seq_weights[seq_id] = 0.
	for i in range(len(trimmed_alignment[list(trimmed_alignment.keys())[0]])):	# iterate over columns in alignment
		freq_counts = [0] * len(amino_acids)
		col = []
		for seq_id in trimmed_alignment:
			if trimmed_alignment[seq_id][i] != '-': # ignore gaps
				freq_counts[aa_to_index[trimmed_alignment[seq_id][i].upper()]] += 1
		num_observed_types = 0
		for j in range(len(freq_counts)):
			if freq_counts[j] > 0: num_observed_types += 1
		for seq_id in trimmed_alignment:
			d = freq_counts[aa_to_index[trimmed_alignment[seq_id][i].upper()]] * num_observed_types
			if d > 0:
				seq_weights[seq_id] += ((0. * len(trimmed_alignment)) + (1. / d)) / len(trimmed_alignment[seq_id])
			else:
				seq_weights[seq_id] += 0.
	return seq_weights


def weighted_freq_count_pseudocount(aa_weigth, pseudocount, amino_acids):
	""" Return the weighted frequency count for a column--with pseudocount."""
	# if the weights do not match, use equal weight
	aa_num = 0
	freq_counts = len(amino_acids)*[pseudocount] # in order defined by amino_acids
	for aa in amino_acids:
		for j in range(len(aa_weigth)):
			if aa_weigth[j][0] == aa:
				freq_counts[aa_num] += 1 * aa_weigth[j][1]
		aa_num += 1
	sum_weights = sum([x[1] for x in aa_weigth])
	for j in range(len(freq_counts)):
		freq_counts[j] = freq_counts[j] / (sum_weights + len(amino_acids) * pseudocount)
	return freq_counts


def weighted_gap_penalty(aa_weigth):
	""" Calculate the simple gap penalty multiplier for the column. If the
	sequences are weighted, the gaps, when penalized, are weighted
	accordingly. """
	# if the weights do not match, use equal weight
	gap_sum = 0.
	for aa, w in aa_weigth:
		if aa == '-':
			gap_sum += w
	return 1 - (gap_sum / sum([x[1] for x in aa_weigth]))


def js_divergence(trimmed_alignment, alignment_dims, sim_matrix, bg_distr, seq_weights, pseudocount, amino_acids, gap_penalty=1):
	""" Return the Jensen-Shannon Divergence for the column with the background
	distribution bg_distr. sim_matrix is ignored. JSD is the default method."""
	distr = bg_distr[:]
	conservation = []
	sequences_order = [seq_id for seq_id in trimmed_alignment]
	for i in range(alignment_dims[1]):	# for each column in alignment
		aa_weigth = [ [trimmed_alignment[seq_id][i], seq_weights[seq_id]] for seq_id in sequences_order ]
		fc = weighted_freq_count_pseudocount(aa_weigth, pseudocount, amino_acids)
		# if background distrubtion lacks a gap count, remove fc gap count
		if len(distr) == 20:
			new_fc = fc[:-1]
			s = sum(new_fc)
			for i in range(len(new_fc)):
				new_fc[i] = new_fc[i] / s
			fc = new_fc
		if len(fc) != len(distr): return -1
		# make r distribution
		r = []
		for x in range(len(fc)):
			r.append(.5 * fc[x] + .5 * distr[x])
		d = 0.
		for y in range(len(fc)):
			if r[y] != 0.0:
				if fc[y] == 0.0:
					d += distr[y] * math.log(distr[y]/r[y], 2)
				elif distr[y] == 0.0:
					d += fc[y] * math.log(fc[y]/r[y], 2)
				else:
					d += fc[y] * math.log(fc[y]/r[y], 2) + distr[y] * math.log(distr[y]/r[y], 2)
		# d /= 2 * math.log(len(fc))
		d /= 2
		if gap_penalty == 1:
			d = d * weighted_gap_penalty(aa_weigth)
		# else:
		# 	return d
		conservation.append(d)
	return conservation


####	MAIN		####
### Obtain motif sequences (SNARE/Habc) from TraceyDB for main groups (Qa.I, Qa.II, Qa.III, Qa.IV, ..., Ha.I, Ha.II, Ha.III, Ha.IV, ...)
### Combine into one fasta Habc + SNARE for those proteins containing both motifs
pairs = {'QI.SNAREa.I': ['QI.Ha.I.fasta', 'QI.Qa.I.fasta'],
		 'QI.SNAREa.I.Syx18': ['QI.Ha.I.Syx18.fasta', 'QI.Qa.I.Syx18.fasta'],
		 'QI.SNAREa.I.Ufe1': ['QI.Ha.I.Ufe1.fasta', 'QI.Qa.I.Ufe1.fasta'],
		 'QII.SNAREa.II': ['QII.Ha.II.fasta', 'QII.Qa.II.fasta'],
		 'QII.SNAREa.II.Sed5': ['QII.Ha.II.fasta', 'QII.Qa.II.Sed5.fasta'],
		 'QII.SNAREa.II.Syx5': ['QII.Ha.II.fasta', 'QII.Qa.II.Syx5.fasta'],
		 'QIII.SNAREa.III': ['QIII.Ha.III.fasta', 'QIII.Qa.III.fasta'],
		 'QIII.SNAREa.III.a': ['QIII.Ha.III.a.fasta', 'QIII.Qa.III.a.fasta'],
		 'QIII.SNAREa.III.b': ['QIII.Ha.III.b.fasta', 'QIII.Qa.III.b.fasta'],
		 'QIV.SNAREa.IV': ['QIV.Ha.IV.fasta', 'QIV.Qa.IV.fasta'],
		 'QIV.SNAREa.IV.Sso': ['QIV.Ha.IV.Sso.fasta', 'QIV.Qa.IV.Sso.fasta'],
		 'QIV.SNAREa.IV.Syx': ['QIV.Ha.IV.Syx.fasta', 'QIV.Qa.IV.Syx.fasta'],
		 
		 'QI.SNAREb.I': ['QI.Hb.I.fasta', 'QI.Qb.I.fasta'],
		 'QII.SNAREb.II': ['QII.Hb.II.fasta', 'QII.Qb.II.fasta'],
		 'QII.SNAREb.II.Bos1': ['QII.Hb.II.Bos1.fasta', 'QII.Qb.II.Bos1.fasta'],
		 'QII.SNAREb.II.Gos1': ['QII.Hb.II.Gos1.fasta', 'QII.Qb.II.Gos1.fasta'],
		 'QII.SNAREb.II.Membrin': ['QII.Hb.II.Membrin.fasta', 'QII.Qb.II.Membrin.fasta'],
		 'QIII.SNAREb.III': ['QIII.Hb.III.fasta', 'QIII.Qb.III.fasta'],
		 'QIII.SNAREb.III.b':['QIII.Hb.III.b.fasta', 'QIII.Qb.III.b.fasta'],
		 'QIII.SNAREb.III.d':['QIII.Hb.III.d.fasta', 'QIII.Qb.III.d.fasta'],
		 
		 'QI.SNAREc.I': ['QI.Hc.I.fasta', 'QI.Qc.I.fasta'],
		 'QIII.SNAREc.III': ['QIII.Hc.III.fasta', 'QIII.Qc.III.fasta'],
		 'QIII.SNAREc.III.b': ['QIII.Hc.III.b.fasta', 'QIII.Qc.III.b.fasta'],
		 'QIII.SNAREc.III.c': ['QIII.Hc.III.c.fasta', 'QIII.Qc.III.c.fasta']
		}


for group in pairs:
	hDict = readFasta('fastas/'+pairs[group][0])
	sDict = readFasta('fastas/'+pairs[group][1])
	if os.path.isfile('fastas/'+group+'.fasta'): continue
	with open('fastas/'+group+'.fasta', 'w') as fo:
		for hid in hDict:
			sid = [x for x in sDict if x.split("|")[0] == hid.split("|")[0]]
			if sid:
				header = ">" + hid.split("|")[0] + "|Habc+SNARE"
				seq = hDict[hid] + sDict[sid[0]]
				fo.write(header+"\n"+seq+"\n")


### Groups & Subgroups
fastaFiles = [x for x in os.listdir('fastas/') if x.endswith(".fasta")]
allGroups = sorted([".".join(g.split(".")[1:-1]) for g in fastaFiles if not g == "notSNARE.fasta"])
# 'Ha.I', 'Ha.II', 'Ha.III', 'Ha.IV', 'Hb.I', 'Hb.II', 'Hb.III', 'Hc.I', 'Hc.III', 'Qa.I', 'Qa.II', 'Qa.III', 'Qa.IV', 'Qb.I', 'Qb.II', 'Qb.III', 'Qc.I', 'Qc.II', 'Qc.III', 'R.I', 'R.II', 'R.III', 'R.IV', 'R.Reg', 'SNAP', 'SNAREa.I', 'SNAREa.II', 'SNAREa.III', 'SNAREa.IV', 'SNAREb.I', 'SNAREb.II', 'SNAREb.III', 'SNAREc.I', 'SNAREc.III'
allMainGroups = [x for x in allGroups if len(x.split(".")) <= 2 and x not in ["SNAPb", "SNAPc"]]
# 'Ha.I.Syx18.I', 'Ha.I.Ufe1.I', 'Ha.III.Ha.III.a', 'Ha.III.Ha.III.b', 'Ha.IV.Sso.IV', 'Ha.IV.Syx.IV', 'Hb.II.Bos1.II', 'Hb.II.Gos1.II', 'Hb.II.Membrin.II', 'Hb.III.Hb.III.b', 'Hb.III.Hb.III.d', 'Hc.III.Hc.III.b', 'Hc.III.Hc.III.c', 'Qa.I.Syx18.I', 'Qa.I.Ufe1.I', 'Qa.II.Sed5.II', 'Qa.II.Syx5.II', 'Qa.III.Qa.III.a', 'Qa.III.Qa.III.b', 'Qa.IV.Sso.IV', 'Qa.IV.Syx.IV', 'Qb.II.Bos1.II', 'Qb.II.Gos1.II', 'Qb.II.Membrin.II', 'Qb.III.Qb.III.b', 'Qb.III.Qb.III.d', 'Qc.III.Qc.III.b', 'Qc.III.Qc.III.c', 'SNAPb', 'SNAPc', 'SNAREa.I.Syx18.I', 'SNAREa.I.Ufe1.I', 'SNAREa.II.Sed5.II', 'SNAREa.II.Syx5.II', 'SNAREa.III.SNAREa.III.a', 'SNAREa.III.SNAREa.III.b', 'SNAREa.IV.Sso.IV', 'SNAREa.IV.Syx.IV', 'SNAREb.II.Bos1.II', 'SNAREb.II.Gos1.II', 'SNAREb.II.Membrin.II', 'SNAREb.III.SNAREb.III.b', 'SNAREb.III.SNAREb.III.d', 'SNAREc.III.SNAREc.III.b', 'SNAREc.III.SNAREc.III.c'
allSubGroups  = [x for x in allGroups if not x in allMainGroups]

# 'SNAREa.I', 'SNAREa.II', 'SNAREa.III', 'SNAREa.IV', 'SNAREb.I', 'SNAREb.II', 'SNAREb.III', 'SNAREc.I', 'SNAREc.III'
fullMainGroups = [x for x in allMainGroups if x.startswith("SNARE")]
# 'Ha.I', 'Ha.II', 'Ha.III', 'Ha.IV', 'Hb.I', 'Hb.II', 'Hb.III', 'Hc.I', 'Hc.III'
habcMainGroups = [x for x in allMainGroups if x.startswith("H")]
# 'Qa.I', 'Qa.II', 'Qa.III', 'Qa.IV', 'Qb.I', 'Qb.II', 'Qb.III', 'Qc.I', 'Qc.II', 'Qc.III', 'R.I', 'R.II', 'R.III', 'R.IV', 'R.Reg', 'SNAP'
snareMainGroups = [x for x in allMainGroups if x not in fullMainGroups + habcMainGroups]

# 'SNAREa.I.Syx18.I', 'SNAREa.I.Ufe1.I', 'SNAREa.II.Sed5.II', 'SNAREa.II.Syx5.II', 'SNAREa.III.SNAREa.III.a', 'SNAREa.III.SNAREa.III.b', 'SNAREa.IV.Sso.IV', 'SNAREa.IV.Syx.IV', 'SNAREb.II.Bos1.II', 'SNAREb.II.Gos1.II', 'SNAREb.II.Membrin.II', 'SNAREb.III.SNAREb.III.b', 'SNAREb.III.SNAREb.III.d', 'SNAREc.III.SNAREc.III.b', 'SNAREc.III.SNAREc.III.c'
fullSubGroups = [x for x in allSubGroups if x.startswith("SNARE")]
# 'Ha.I.Syx18.I', 'Ha.I.Ufe1.I', 'Ha.III.Ha.III.a', 'Ha.III.Ha.III.b', 'Ha.IV.Sso.IV', 'Ha.IV.Syx.IV', 'Hb.II.Bos1.II', 'Hb.II.Gos1.II', 'Hb.II.Membrin.II', 'Hb.III.Hb.III.b', 'Hb.III.Hb.III.d', 'Hc.III.Hc.III.b', 'Hc.III.Hc.III.c'
habcSubGroups = [x for x in allSubGroups if x.startswith("H")]
# 'Qa.I.Syx18.I', 'Qa.I.Ufe1.I', 'Qa.II.Sed5.II', 'Qa.II.Syx5.II', 'Qa.III.Qa.III.a', 'Qa.III.Qa.III.b', 'Qa.IV.Sso.IV', 'Qa.IV.Syx.IV', 'Qb.II.Bos1.II', 'Qb.II.Gos1.II', 'Qb.II.Membrin.II', 'Qb.III.Qb.III.b', 'Qb.III.Qb.III.d', 'Qc.III.Qc.III.b', 'Qc.III.Qc.III.c', 'SNAPb', 'SNAPc'
snareSubGroups = [x for x in allSubGroups if x not in fullSubGroups + habcSubGroups]


### Align and trim sequences
fastaFiles = os.listdir('fastas/')
for f in fastaFiles:
	if not f.endswith('.fasta'): continue
	if f == 'notSNARE.fasta': continue
	alignmentOut = os.path.join('alignments', f.replace(".fasta", ".aln"))
	alignmentTrimmed = os.path.join('alignments', f.replace(".fasta", ".trimmed.aln"))
	if os.path.isfile(alignmentOut): continue
	alignment = alignFasta(os.path.join('fastas', f), alignmentOut)
	trimmeAln = trimAlignment(alignmentOut, alignmentTrimmed, col_threshold=0.8, row_threshold=0.8)


### Build trees from alignments // Use this code or IQtree website: http://iqtree.cibiv.univie.ac.at/
for aln in [f for f in os.listdir('alignments') if f.endswith('.trimmed.aln')]:
	if not aln.replace(".trimmed.aln", "").split(".", 1)[1] in allMainGroups: continue
	# Check for finish line in log file
	try:
		log = open(('trees/%s'%(aln.replace(".trimmed.aln", ".log"))), 'r').readlines()[-1]
	except FileNotFoundError:
		log = "notFound"
	if log.startswith("Date"):
		print("\t%s: completed" % (aln.replace(".trimmed.aln", ".treefile")))
		continue
	print(aln)
	cmd = 'iqtree -s alignments/%s --seqtype AA -m TEST -bb 1000 -alrt 1000 --prefix trees/%s ' % (aln, aln.replace(".trimmed.aln", ""))
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = process.communicate()


### Replace SequenceIDs in trees to match IDs in fastas
# Names changed to XXX_XXX are replaces with XXX|XXX
for f in os.listdir('trees'):
	if f.endswith('.treefile'):
		fasta = readFasta('fastas/' + f.replace('.treefile', '.fasta'))
		tree = open('trees/' + f, 'r').readlines()[0]
		for k in fasta:
			ktree = k.replace('|', '_')
			tree = tree.replace(ktree, k)
		with open('trees/'+f, 'w') as t:
			t.write(tree)


# Names changed to XXX_XXX_XXX are replaces with XXX|XXX,XXX
#for f in os.listdir('trees'):
for f in ["QIV.SNAP.treefile"]:
	if f.endswith('.treefile'):
		fasta = readFasta('fastas/' + f.replace('.treefile', '.fasta'))
		tree = open('trees/' + f, 'r').readlines()[0]
		# Read log file to retrieve changes on names
		nameChanges = {}
		skipLine = 1
		for line in open('trees/' + f.replace('.treefile', '.log'), 'r').readlines():
			if line.startswith('WARNING: Some sequence names are changed as follows:'):
				skipLine = 0
				continue
			if line == "\n" and skipLine == 0: skipLine = 1
			if skipLine: continue
			changes = line.strip().split(' -> ')
			nameChanges[changes[1]] = changes[0].replace(',', '_')
		# Replace names in tree
		for k in nameChanges:
			tree = tree.replace(k, nameChanges[k])
		with open('trees/'+f, 'w') as t:
			t.write(tree)


### Use trees to extract SequenceIDs for subgroups (Syx18.I, Ufe1.I, Syx5.II, ...)
### Align and filter new fasta files
idFiles = os.listdir('idsForSubgroups/')
for f in idFiles:
	if f.endswith('.txt'):
		# Read fasta file and Ids for subgroup
		root = f.replace('.txt', '')
		mainGroup = root.split('.')[0]
		group = ".".join(root.split('.')[1:3])
		subgroup = ".".join(root.split('.')[3:])
		fasta = readFasta('fastas/' + ".".join([mainGroup, group]) + '.fasta')
		alignmentOut = os.path.join('alignments', root + '.aln')
		alignmentTrimmed = os.path.join('alignments', root + '.trimmed.aln')
		if os.path.isfile(alignmentTrimmed): continue
		# Extract sequences for subgroup
		subgroupIds = [x.strip() for x in open('idsForSubgroups/' + f, 'r').readlines()]
		subgroupFasta = {}
		for sgId in subgroupIds:
			subgroupFasta[sgId] = fasta[sgId]
		# Write subgroup fasta file
		subgroupFastaFile = 'fastas/' + root + '.fasta'
		with open(subgroupFastaFile, 'w') as fo:
			for k in subgroupFasta:
				nullhandler = fo.write('>%s\n%s\n' % (k, subgroupFasta[k]))
		# Align and trim subgroup fasta file
		print(root)
		alignment = alignFasta(subgroupFastaFile, alignmentOut)
		trimmeAln = trimAlignment(alignmentOut, alignmentTrimmed, col_threshold=0.8, row_threshold=0.8)


### Generate HMMs for each group and subgroup // and create hmmDB
alignments = [f for f in os.listdir('alignments') if f.endswith('.trimmed.aln')]
for aln in alignments:
	motif = aln.replace(".trimmed.aln", "").split(".", 1)[1]
	print(motif)
	buildHMM('alignments/%s'%(aln), motif, 'hmms/%s.hmm'%(motif))


input_files = ['hmms/'+x for x in os.listdir('hmms') if 'hmmDb.hmm' not in x]
cmd = ['cat'] + input_files
with open('hmms/hmmDb.hmm', "w") as outfile:
	subprocess.run(cmd, stdout=outfile)


for f in [x for x in os.listdir('hmms') if 'Db' in x and not x.endswith('.hmm')]:
	os.remove('hmms/%s'%(f))


cmdPress = 'hmmpress hmms/hmmDb.hmm'
process = subprocess.Popen(cmdPress.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, error = process.communicate()


## Search for HMMs using all sequences
hmmDb = 'hmms/hmmDb.hmm'
for fasta in os.listdir('fastas/'):
	name = ".".join(fasta.split('.')[1:-1])
	print(name)
	domtblout = 'regressionModel/' + name + '.domtblout'
	cmd = "hmmsearch --domtblout %s %s %s" % (domtblout, hmmDb, 'fastas/' + fasta)
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








