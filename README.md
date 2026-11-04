# dpi_domain

Thesis experiments on drug-protein interaction (DPI) prediction with a focus on protein domains.
The idea: a compound that binds a domain on a single-domain protein does not necessarily bind the same
domain when it sits in a multi-domain protein ("context-sensitive" interactions). This repo builds
test sets from BindingDB that isolate such cases, and evaluates baselines and published predictors on
them.

## Layout

All code is under `experiments/`:

| Path | What it does |
| --- | --- |
| `building_datasets.ipynb` | Downloads BindingDB, binarizes affinities, builds the test/validation/training sets and writes them to `datasets/` |
| `baselines.ipynb`, `baselines.py` | Two sampling baselines (`baseline_P`: random affinity of the same protein, `baseline_C`: random affinity of the same compound), 100 repetitions, reported for binary and affinity targets |
| `evaluate_predictors.ipynb`, `report_preds.py` | Scores external predictors (MolTrans, DeepDTA) from prediction files in `predictions/` against the test sets |
| `evaluate.py` | Metric helpers: AUC, AUPR, MCC, F1, sensitivity, specificity (binary) and CI, MAE, MSE, Pearson, Spearman (affinity) |
| `save_dataset.py` | Writes a pair list to CSV with sequences and SMILES (rows are ordered by an MD5 of the pair, so the shuffle is reproducible) |
| `data_interfaces/access_bindingdb.py` | `bindingdb_records`: downloads/parses the BindingDB TSV, keeps single-chain targets with a UniProt ID and PubChem CID, binarizes Ki/Kd/IC50/EC50, caches to pickle |
| `data_interfaces/access_pfam.py` | `pfam_records`: Pfam domain annotations per protein, single-domain test, segment-string redundancy reduction |
| `data_interfaces/access_pubchem.py` | Canonical SMILES from the PubChem PUG REST API, batched and cached |
| `data_interfaces/access_uniprot.py` | Protein sequences from UniProt, cached |
| `data_interfaces/download_file.py` | `urlretrieve` with a tqdm progress bar |

## Datasets produced

`building_datasets.ipynb` writes these CSVs to `experiments/datasets/`, all with the columns
`protein_uniprot_id, compound_pubchem_cid, protein_sequence, compound_SMILES, binary_interaction, affinity[-log10(Ki/10e9)]`.
Every export calls `save_dataset_csv` with its default `only_with_ki=True`, so a pair is written only when
BindingDB has a Ki value for it; pairs whose binary label comes from Kd/IC50/EC50 alone are dropped rather
than written with a `none` affinity. The pair lists are selected and the positive/negative ratios computed
before this filter, so the ratios in the files can differ from those described below:

- `testset1.csv`: the difficult set. Domain-compound pairs where the compound binds a single-domain protein
  carrying the domain but does not bind some other protein carrying the same domain, after merging proteins
  with the same domain/non-domain segment layout.
- `testset2.csv`: all BindingDB pairs involving a protein or compound from testset1, resampled to testset1's
  positive/negative ratio.
- `testset3.csv`: a random 10% of all negative pairs plus positives at testset1's ratio.
- `comb123.csv`: union of the three test sets (for exclusion from training).
- `valset.csv`: a validation sample of the remaining pairs.
- `trainset_allrest.csv` and `trainset_same_ratio.csv`: the remaining pairs, all of them or resampled to the
  same ratio.

Binarization: an affinity below 1000 nM (1 uM) is an interaction, above 30000 nM (30 uM) a
non-interaction, and a pair is kept only when all its Ki/Kd/IC50/EC50 records agree
(`use_measure="all_consistent"`). Affinity labels are pKi from the mean Ki.

## Running

The notebooks and scripts use paths relative to `experiments/`, so run everything from there:

```
cd experiments
jupyter notebook
```

Dependencies are listed in `requirements.txt` (`pip install -r requirements.txt`). `enviornment.yml` is
the exported conda environment the work was done in (Python 3.9); its `name`/`prefix` point at the
original machine, so pass `-n <name>` if you create it with `conda env create`.

Order: `building_datasets.ipynb` first (it downloads the full BindingDB TSV, a multi-GB file, into
`data_sources/bindingdb/` and caches the parsed records under `data_interfaces/cached_data/`; the first
pass takes a long time), then `baselines.ipynb` and `evaluate_predictors.ipynb`.

## Not in the repo

`datasets/`, `data_sources/`, `predictions/` and `data_interfaces/cached_data/` are git-ignored. The
prediction files for `evaluate_predictors.ipynb` (`predictions/moltrans_custom.csv`,
`predictions/deepdta_custom.csv`; lines of `<uniprot>_<cid>,<score>`) come from running those models
elsewhere and are not generated here.

Two of the web services the data interfaces were written against have changed since: the legacy Pfam
site (`pfam.xfam.org`, now folded into InterPro) and UniProt's old `uploadlists` ID-mapping endpoint. The
PubChem and BindingDB downloads should still work.
