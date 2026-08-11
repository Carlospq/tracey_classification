import sys, os, subprocess
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from scripts.pipeline_utils import searchHMMs, retriveEvalues


hmm_groups = [["Ha.I", "Ha.I"],
			  ["Ha.II", "Ha.II"],
			  ["Ha.III", "Ha.III"],
			  ["Ha.III.a", "Ha.III.Ha.III.a"],
			  ["Ha.III.b", "Ha.III.Ha.III.b"],
			  ["Ha.IV", "Ha.IV"],
			  ["Ha.IV.Sso", "Ha.IV.Sso.IV"],
			  ["Ha.IV.Syx", "Ha.IV.Syx.IV"],
			  ["Hb.I", "Hb.I"],
			  ["Hb.II", "Hb.II"],
			  ["Hb.II.Gos", "Hb.II.Hb.II.a"],
			  ["Hb.II.Membos", "Hb.II.Hb.II.b"],
			  ["Hb.III", "Hb.III"],
			  ["Hc.I", "Hc.I"],
			  ["Hc.III", "Hc.III"],
			  ["Hc.III.b", "Hc.III.Hc.III.b"],
			  ["Hc.III.c", "Hc.III.Hc.III.c"],
			  ["Qa.I", "Qa.I"],
			  ["Qa.II", "Qa.II"],
			  ["Qa.III", "Qa.III"],
			  ["Qa.III.a", "Qa.III.Qa.III.a"],
			  ["Qa.III.b", "Qa.III.Qa.III.b"],
			  ["Qa.IV", "Qa.IV"],
			  ["Qb.I", "Qb.I"],
			  ["Qb.II", "Qb.II"],
			  ["Qb.III", "Qb.III"],
			  ["Qb.III.b", "Qb.III.Qb.III.b"],
			  ["Qb.III.d", "Qb.III.Qb.III.d"],
			  ["Qc.I", "Qc.I"],
			  ["Qc.II", "Qc.II"],
			  ["Qc.III", "Qc.III"],
			  ["Qc.III.b", "Qc.III.Qc.III.b"],
			  ["Qc.III.c", "Qc.III.Qc.III.c"],
			  ["R.I", "R.I"],
			  ["R.II", "R.II"],
			  ["R.III", "R.III"],
			  ["R.IV", "R.IV"],
			  ["R.Reg", "R.Reg"],
			  ["SNAP.b", "SNAPb"],
			  ["SNAP.c", "SNAPc"]]

new_domtblout = "results/hmmscan_domtblout_new.txt"
old_domtblout = "results/hmmscan_domtblout_old.txt"

with open('results/hmmscan_results_summary.txt', 'w') as fo_summary:
	with open('results/hmmscan_results_all.txt', 'w') as fo_all:

		# HEADERS IN OUTPUTS
		fo_all.write("\t".join(["seq_id", "hmm_name", "model", "pval"]) + "\n")
		fo_summary.write("\t".join(["group", "sequences", "matches_in_new", "matches_in_old",
							"mean_diff_better", "len_better",
							"mean_diff_worse", "len_worse"]) + "\n")

		for group in hmm_groups:

			results = {}

			hmm_new = 'hmms/' + group[0] + ".hmm"
			hmm_old = 'traceyHmms/' + group[1] + ".hmm"
			in_fasta = 'initialSetSelection/iterations/' + group[0] + ".fasta"

			searchHMMs(in_fasta, hmm_new, new_domtblout)
			searchHMMs(in_fasta, hmm_old, old_domtblout)

			new_evals = retriveEvalues(new_domtblout, [group[0], group[1]])
			old_evals = retriveEvalues(old_domtblout, [group[0], group[1]])

			# Fix names for SNAP.b, SNAP.c, Qb.III and Qc.III
			if group[1] == "SNAPb":
				group[1] = "SNAP.b"
			if group[1] == "SNAPc":
				group[1] = "SNAP.c"
			if group[1] == "Qa.III.Qa.III.a":
				group[1] = "Qa.III.a"
			if group[1] == "Qa.III.Qa.III.b":
				group[1] = "Qa.III.b"
			if group[1] == "Qb.III.Qb.III.b":
				group[1] = "Qb.III.b"
			if group[1] == "Qb.III.Qb.III.d":
				group[1] = "Qb.III.d"
			if group[1] == "Qc.III.Qc.III.b":
				group[1] = "Qc.III.b"
			if group[1] == "Qc.III.Qc.III.c":
				group[1] = "Qc.III.c"

			keys = set(new_evals.keys() | old_evals.keys())
			for k in keys:
				results[k] = {"new": 10e10, "old": 10e10}
				if k in new_evals:
					results[k]["new"] = new_evals[k][group[0]]
				if k in old_evals:
					results[k]["old"] = old_evals[k][group[1]]

			# DETAILED RESULTS
			for k in results:
				if results[k]["new"] < 10e10:
					fo_all.write("\t".join( [k, group[0], "new", str(results[k]["new"])] ) + "\n")
				if results[k]["old"] < 10e10:
					fo_all.write("\t".join( [k, group[0], "old", str(results[k]["old"])] ) + "\n")

			# SUMMARY
			with open(in_fasta, 'r') as f:
				sequences = []
				for line in f:
					if line.startswith('>'):
						sequences.append(line[1:].strip())

			diff_better = []
			diff_worse = []
			for k in results:
				if -np.log(results[k]["new"]) > -np.log(results[k]["old"]):
					diff_better.append(-np.log(results[k]["new"]) + np.log(results[k]["old"]))
				else:
					diff_worse.append(-np.log(results[k]["old"]) + np.log(results[k]["new"]))

			fo_summary.write("\t".join([group[0],
								str(len(sequences)), str(len(new_evals)), str(len(old_evals)),
								str(np.mean(diff_better)), str(len(diff_better)),
								str(np.mean(diff_worse)), str(len(diff_worse))]) + "\n")
