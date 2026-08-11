new_hmms_file = "results/new_nseqs_per_hmm.txt"
old_hmms_file = "results/old_nseqs_per_hmm.txt"
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
			  ["Qa.III.a", "Qa.III.a"],
			  ["Qa.III.b", "Qa.III.b"],
			  ["Qa.IV", "Qa.IV"],
			  ["Qb.I", "Qb.I"],
			  ["Qb.II", "Qb.II"],
			  ["Qb.III", "Qb.III"],
			  ["Qb.III.b", "Qb.III.b"],
			  ["Qb.III.d", "Qb.III.d"],
			  ["Qc.I", "Qc.I"],
			  ["Qc.II", "Qc.II"],
			  ["Qc.III", "Qc.III"],
			  ["Qc.III.b", "Qc.III.b"],
			  ["Qc.III.c", "Qc.III.c"],
			  ["R.I", "R.I"],
			  ["R.II", "R.II"],
			  ["R.III", "R.III"],
			  ["R.IV", "R.IV"],
			  ["R.Reg", "R.Reg"],
			  ["SNAP.b", "SNAP.b"],
			  ["SNAP.c", "SNAP.c"]]

nseqs, nseqs_new, nseqs_old = {}, {}, {}
with open(new_hmms_file, "r") as f:
	for line in f.readlines():
		parts = line.strip().split()
		key = parts[0]
		val = parts[1]
		if key == "NAME":
			name = val
			nseqs_new[name] = 0
		else:
			nseqs_new[name] = val

with open(old_hmms_file, "r") as f:
	for line in f.readlines():
		parts = line.strip().split()
		key = parts[0]
		val = parts[1]
		if key == "NAME":
			name = val
			nseqs_old[name] = {"new": 0, "old": 0}
		else:
			nseqs_old[name] = val

for group in hmm_groups:
	hmm_new = group[0]
	hmm_old = group[1]
	nseqs[hmm_new] = {'new': nseqs_new[hmm_new], 'old': nseqs_old[hmm_old]}



with open("results/seqs_per_hmm_table.txt", "w") as out:
	out.write("HMM\tNew_Seqs\tOld_Seqs\n")
	for x in [[hmm, nseqs[hmm]['new'], nseqs[hmm]['old']] for hmm in nseqs if int(nseqs[hmm]['new']) > 0 and int(nseqs[hmm]['old']) > 0]:
		out.write( f'%s\t%s\t%s\n' % (x[0], str(x[1]), str(x[2])) )

# New	Old				(In common 40; Unique new 92)
# Ha
# Ha.I	Ha.I
# Ha.I.Syx18
# Ha.I.Syx18.Metazoa
# Ha.I.Syx18.Viridiplantae
# Ha.I.Ufe1
# Ha.II	Ha.II
# Ha.III	Ha.III
# Ha.III.a	Ha.III.a
# Ha.III.b	Ha.III.b
# Ha.III.b.Fungi
# Ha.III.b.Metazoa
# Ha.III.b.Syx17
# Ha.III.b.Viridiplantae
# Ha.IV	Ha.IV
# Ha.IV.Metazoa
# Ha.IV.Viridiplantae
# Ha.IV.Syx1x
# Ha.IV.Sso	Ha.IV.Sso
# Ha.IV.Syx	Ha.IV.Syx
# Hb
# Hb.I	Hb.I
# Hb.I.Metazoa
# Hb.I.Fungi
# Hb.I.Viridiplantae
# Hb.II	Hb.II
# Hb.II.Gos	Hb.II.a
# Hb.II.Gos.Fungi
# Hb.II.Gos.Zoophyta
# Hb.II.Membos	Hb.II.b
# Hb.II.Membos.Bos1
# Hb.II.Membos.Membrin
# Hb.III	Hb.III
# Hb.III.b
# Hb.III.b.Vti1a
# Hb.III.b.Vti1b
# Hb.III.b.Vtifungi
# Hb.III.b.Vtiviridiplantae
# Hc
# Hc.I	Hc.I
# Hc.I.Fungi
# Hc.I.Metazoa
# Hc.I.Viridiplantae
# Hc.III	Hc.III
# Hc.III.b	Hc.III.b
# Hc.III.b.Fungi
# Hc.III.b.Metazoa
# Hc.III.c	Hc.III.c
# Hc.III.c.Metazoa
# Hc.III.c.Fungi
# Hc.III.c.Vam7
# Qa
# Qa.I	Qa.I
# Qa.I.Syx18
# Qa.I.Syx18.Metazoa
# Qa.I.Syx18.SAR
# Qa.I.Syx18.Viridiplantae
# Qa.I.Ufe1
# Qa.II	Qa.II
# Qa.II.Sed5
# Qa.II.Syx5
# Qa.II.Syx5.Metazoa
# Qa.II.Syx5.SAR
# Qa.II.Syx5.Viridiplantae
# Qa.III	Qa.III
# Qa.III.a	Qa.III.a
# Qa.III.b	Qa.III.b
# Qa.IV	Qa.IV
# Qa.IV.Sso
# Qa.IV.Syx
# Qb
# Qb.I	Qb.I
# Qb.I.Metazoa
# Qb.I.Viridiplantae
# Qb.II	Qb.II
# Qb.II.Gos
# Qb.II.Gos.Fungi
# Qb.II.Gos.Viridiplantae
# Qb.II.Gos.Gos2
# Qb.II.Membos
# Qb.II.Membos.Bos1
# Qb.II.Membos.Membrin
# Qb.III	Qb.III
# Qb.III.b	Qb.III.b
# Qb.III.b.Fungi
# Qb.III.b.Metazoa
# Qb.III.b.SAR
# Qb.III.b.Viridiplantae
# Qb.III.d	Qb.III.d
# Qc
# Qc.I	Qc.I
# Qc.I.Metazoa
# Qc.I.Fungi
# Qc.I.Viridiplantae
# Qc.II	Qc.II
# Qc.II.Bet1
# Qc.II.Sft1
# Qc.III	Qc.III
# Qc.III.b	Qc.III.b
# Qc.III.b.Metazoa
# Qc.III.b.Syo7
# Qc.III.c	Qc.III.c
# Qc.III.c.Syp5
# Qc.III.c.Syp8
# Qc.III.c.Vam7
# R
# R.I	R.I
# R.I.Metazoa
# R.I.Fungi
# R.I.Viridiplantae
# R.II	R.II
# R.III	R.III
# R.III.Endob
# R.III.Nyv1
# R.III.Vamp4
# R.III.Vamp7
# R.IV	R.IV
# R.IV.Metazoa
# R.IV.Fungi
# R.Reg	R.Reg
# R.Reg.Metazoa
# R.Reg.Fungi
# R.Reg.Viridiplantae
# SNAP
# SNAPb	SNAPb
# SNAPb.Sec9
# SNAPb.SN25
# SNAPb.SN29
# SNAPc	SNAPc
# SNAPc.Sec9
# SNAPc.SN25
# SNAPc.SN29










