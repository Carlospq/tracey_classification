# tracey_classification

Classification pipeline for SNARE and Habc protein motifs. Sequences are aligned, used to
build phylogenetic trees and HMM profiles, and then classified into a hierarchy of
domains/groups/subgroups (e.g. `Qa` -> `Qa.I` -> `Qa.I.Syx18` -> ...) using logistic-regression
models trained on HMM e-values.

This repository accompanies a manuscript currently under review. It contains the analysis code;
the underlying sequence/taxonomy data is not published here (see
[Input data](#input-data-not-included-in-this-repository) below).

## Pipeline overview

1. **Sequences** (`fastas/*.fasta`) are cleaned and aligned with [MUSCLE](https://drive5.com/muscle5/).
2. Alignments (`alignments/`) are trimmed (gap-column/row filtering) and used to build
   maximum-likelihood trees with [IQ-TREE](http://www.iqtree.org/) (`trees/`).
3. Trees are used to manually curate sequence subgroups (`idsForSubgroups/`,
   `initialSetSelection/`), which are aligned and turned into
   [HMMER](http://hmmer.org/) profile HMMs (`hmms/`).
4. All sequences are searched against the HMM database (`hmmsearch`) and the resulting e-values
   are used to train a hierarchy of scikit-learn `LogisticRegression` classifiers
   (`regressionModel/lm_*.sav`).
5. New sequences are classified by running them through the same HMM search and walking the
   trained classifier hierarchy (`scripts/domainPredict.py`).

## Repository layout

| Path | Contents |
|---|---|
| `scripts/` | Core pipeline code (alignment, tree-building, HMM training, prediction) and R plotting scripts |
| `scripts/pipeline_utils.py` | Shared helper functions (FASTA I/O, alignment, HMM search, e-value parsing) used across the pipeline scripts |
| `scripts/motifNames.py` | The SNARE/Habc domain hierarchy (`SNARETree`, `HabcTree`) used to drive both training and prediction |
| `initialSetSelection/` | Iterative procedure for selecting a representative "core" sequence set per subgroup |
| `results/` | Scripts comparing this pipeline's HMMs against a legacy HMM database |
| `utils/` | Taxonomy reference tables and R/Shiny tree-visualization apps |
| `fastas/`, `alignments/`, `hmms/`, `trees/`, `idsForSubgroups/`, `regressionModel/`, `initialSetSelection/initialSets/`, `initialSetSelection/iterations/` | Pipeline inputs/outputs. Empty in this repository (see below) but kept as placeholders (`.gitkeep`) so the expected layout is visible |

## Installation

**Python** (developed with Python 3):
```
pip install -r requirements.txt
```

**R packages**: no lockfile is used in this repository; install what each script needs, e.g.:
```r
install.packages(c("shiny", "bslib", "ggplot2", "ggtree", "dplyr", "stringr", "tidytree", "DT",
                    "seqinr", "reshape2", "ggdendro", "gridExtra", "cowplot", "optparse",
                    "tidyr", "tibble", "ggseqlogo", "ggrepel", "patchwork", "ggpubr",
                    "tidyverse", "ape"))
```
(`ggtree` is on Bioconductor: `BiocManager::install("ggtree")`.)

**External binaries** (must be on `PATH`):
- [`muscle`](https://drive5.com/muscle5/) — sequence alignment
- [HMMER](http://hmmer.org/) (`hmmbuild`, `hmmsearch`, `hmmpress`) — profile HMM building/search
- [IQ-TREE](http://www.iqtree.org/) (`iqtree` / `iqtree2`) — phylogenetic tree inference
- `Rscript` — for the R plotting steps invoked from Python

## Running the pipeline

`scripts/generateStatisticalHmmModels.py` is the authoritative training pipeline (clean fastas ->
align -> trim -> build trees -> plot trees -> manually curate initial sequence sets via
`initialSetSelection/setSelection.py` -> copy the resulting fastas/alignments/HMMs back into the
top-level folders -> build per-motif and combined HMMs -> `hmmsearch` against all sequences and a
negative (`notSNARE`) set -> train and evaluate the logistic-regression models).

**It is not meant to be run start-to-finish in one go.** The alignment and tree-building steps
are flagged in the script's own comments as done on an HPC cluster — they're
computationally expensive (MUSCLE on large sequence sets, IQ-TREE with 1000 bootstraps) and are
normally run separately on a cluster rather than interactively. Every step in the script checks
whether its output file already exists and skips it if so, so the intended workflow is:

1. Run (or submit to a cluster) the alignment/tree-building commands for a given fasta ahead of
   time, so `alignments/*.trimmed.aln` and `trees/*.treefile` are already populated.
2. Then run the script from the repository root — it will skip the steps whose outputs already
   exist and continue from wherever the pipeline actually left off:
   ```
   python -m scripts.generateStatisticalHmmModels
   ```

`scripts/func.py` is an earlier, exploratory version of the same stages, kept for reference; it
is not the one used to produce the final models.

Note the tree-plotting step shells out to `Rscript scripts/PlotTree.R`; if that fails to run
non-interactively in your environment, open and run `scripts/PlotTree.R` manually in RStudio
(working directory: repository root) instead.

Most R scripts in `scripts/`, `results/` and `utils/` assume a specific working directory rather
than taking it as an argument — each file states which at the top:
- Repo root: `scripts/PlotTree.R`, `scripts/MutualInformation.R`, `scripts/newVSoldPvalues.R`,
  `scripts/plotEvaluesDistributionsHMMs.R`, `utils/shinyTrees.R`
- `initialSetSelection/`: `initialSetSelection/distanceMatrix.R`, `utils/plotTrees.R`
- `results/`: `results/plot_pvalues_hmmscan.R`

## Prediction / inference

`scripts/domainPredict.py` classifies either a FASTA file or a single sequence, using the trained
models in `regressionModel/`:
```python
from scripts.domainPredict import predictFromFasta, predictFromSeq
from scripts.motifNames import SNARETree

predictFromFasta("my_sequences.fasta", all_motifs, SNARETree)
predictFromSeq("MSTNL...", all_motifs, SNARETree)
```
Its `__main__` block runs a demo prediction against `test.fasta` (not included in this
repository) using `hmms/SNAREDb.hmm`; supply your own file to try it.

## Input data (not included in this repository)

The data below is required to actually run the pipeline but is not published here (see
`.gitignore`). The directories exist as empty placeholders so the expected structure is visible.
Inputs are split from outputs below: only the first table needs to be supplied from outside this
repository — everything in the second table is produced by running the pipeline itself.

### Required to start the pipeline

| Path | What it is | Status here |
|---|---|---|
| `fastas/<MainGroup>.fasta` (e.g. `Qa.fasta`), `fastas/notSNARE.fasta` | Raw per-main-group sequences and the negative training set, manually exported from the private TraceyDB (`generateStatisticalHmmModels.py` comment: *"step done manually"*) | Not included |
| `initialSetSelection/initialSets/*.initialSet.txt` | Manually curated seed sequence-ID lists per subgroup, the starting point for the iterative set-selection procedure | Not included |
| `utils/taxonomyTable.csv`, `utils/sequenceIDtaxonomyID.csv` | Taxonomy lookup tables used by the R tree-plotting scripts | **Included** (static exports; originally produced by an internal script querying a private database, which is not part of this repository) |

Optional, only needed for specific side scripts (not the core train/predict pipeline):

| Path | Needed for | Status here |
|---|---|---|
| `regressionModel/lm_*.sav` | Running `domainPredict.py` for inference **without** first retraining (see next table — these are training outputs too) | Not included |
| `traceyHmms/traceyHmmDb.hmm` | Legacy HMM database comparison (`scripts/HMMcomparison.py`, `results/hmmscan_to_table.py`) | Not included; provenance predates this repository |
| `results/new_nseqs_per_hmm.txt`, `results/old_nseqs_per_hmm.txt` | `results/seqs_to_table.py` reporting | Not included; no generator script exists in this repository |
| `test.fasta` (repo root) | `domainPredict.py`'s demo `__main__` block | Not included |

### Generated by the pipeline (do not copy — reproducible by running the code)

| Path | Produced by | From |
|---|---|---|
| `fastas/<Group>.<Subgroup>.fasta`, `idsForSubgroups/*.txt` | The "Copying files" step in `generateStatisticalHmmModels.py` | Copied/derived from `initialSetSelection/iterations/*.fasta` |
| `initialSetSelection/iterations/*.fasta/.aln/.hmm/.results/.data.csv` | `initialSetSelection/setSelection.py` | `initialSetSelection/initialSets/*.initialSet.txt` |
| `alignments/*.aln`, `*.trimmed.aln` | `alignFasta`/`trimAlignment` (`scripts/pipeline_utils.py`) — typically run separately/on a cluster, see [Running the pipeline](#running-the-pipeline) | `fastas/*.fasta` |
| `trees/*.treefile`, `*.log` | IQ-TREE — typically run separately/on a cluster, see [Running the pipeline](#running-the-pipeline) | `alignments/*.trimmed.aln` |
| `hmms/*.hmm`, `hmms/SNAREDb.hmm` (+ `.h3*` press files) | `hmmbuild`/`hmmpress` | `alignments/*.trimmed.aln` |
| `regressionModel/*.domtblout` | `hmmsearch` | `hmms/SNAREDb.hmm` + `fastas/*.fasta` |
| `regressionModel/lm_*.sav`, `report_*.txt`, `confusionMatrix_*.png` | Logistic-regression training in `generateStatisticalHmmModels.py` | `regressionModel/*.domtblout` |

## Known limitations

- Only SNARE classification models are currently trained; Habc classification is not yet wired
  up (`domainGroups["Habc"]` is an empty placeholder in `generateStatisticalHmmModels.py`, and
  `HabcTree` in `scripts/motifNames.py` is unused until that training path exists).
- `scripts/generateStatisticalHmmModels.py` computes e-values for the training data slightly
  differently (`pipeline_utils.retriveEvalues`) than `scripts/domainPredict.py` does at inference
  time (`pipeline_utils.retriveEvaluesLog`); this divergence predates this cleanup pass and has
  not been independently re-verified as intentional.
- There is no automated test suite or CI for this repository.

## License

[MIT](LICENSE)

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon
publication).
