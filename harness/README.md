# NERVE comparison harness

Tools for comparing two runs of the pipeline without duplicating any pipeline
code. They exist so that a change can be shown to preserve results, rather than
asserted to.

- **Per-stage snapshots.** Running with `--keep_snapshots True` pickles the
  protein list after every stage into `<working_dir>/snapshots/NN_name.pkl`.
  `compare_states.py` diffs two snapshot sets attribute by attribute, matched by
  accession, reporting which stage, which attribute, which proteins and by how
  much they differ. This localises the *first* stage at which two runs diverge.
- **Final candidates.** `compare_candidates.py` diffs two
  `vaccine_candidates.csv` as biology: set overlap and Jaccard, score deltas,
  Pearson/Spearman correlation, and top-N overlap.

Both exit `0` when equal and `1` otherwise, so they can gate a regression script.
Duplicate or missing protein identifiers fail comparison rather than being
silently collapsed.

## Reduced same-source validation

`harness/validation/run_all.sh` checks image provenance, focused tests, an
interrupted and resumed small run against a fresh run, and selected *C. jejuni*
BLAST-thread and epitope-worker variants. It requires a built optimized image
and a running PSORTb service on `nerve-network`; set `IMAGE` to the verified
image ID and optionally set `RUN_ROOT` to a fresh directory. The checkpoint
case uses 100 proteins so the process is still running after the first saved
stage. Its snapshot comparison is exact except for equal-score ordering among
secondary PSORTb localization entries; the selected first localization and all
scores remain exact. The preflight pins both *C. jejuni* input hashes. It does
not establish equivalence with the pre-Zeamays
baseline or measure large-proteome scaling; those are separate validation
stages.

For the Section D gate, set `SECTION_D_REFERENCE` to the accepted Section C
`cjejuni/reference-t1-ew1` directory when running `run_all.sh`. The extra check
requires successful run manifests, consistent selection reasons and output
hashes, an interrupted first checkpoint attempt, production logs, and
byte-identical candidate and discarded CSVs against that reference. It writes
`section_d_gate.json` under `RUN_ROOT` only after every check passes.

## Tolerance

**The default is `--tol 0`, exact equality, and that is the intended setting.**

`code/NERVE.py` pins `OMP_NUM_THREADS` and the other numeric thread pools to 1
at import, before any of them is created. That matters because the adhesin and
virulence models are float32: a multi-threaded reduction reorders the sums and
moves `p_ad` / `p_vir` by around 1e-7, which is enough to reorder candidates of
similar score between two runs of *identical* code. Pinning removes that
variation rather than bounding it, so there is no epsilon to budget for.

Other compared values are required to match exactly. In particular
the peptide sums are integer arithmetic (`score += 1` or `score += len(...)`,
divided once by the protein length at the end), so they carry no float
accumulation error at all.

A difference of ~1e-7 in `p_ad` or `p_vir` therefore means the thread pinning did
not take effect in one of the runs -- not that the tolerance is too tight. Larger
differences, or any difference in a discrete attribute, indicate a real change.
Use `--tol` to characterise a difference you have already found, not to make one
go away.

## Ordering

Two attributes are order-sensitive in a way that is easy to miss.

The peptide sums fold `list_of_shared_*_peps` pairwise, comparing each peptide
with its predecessor and adding either 1 or the full peptide length. Reordering
that list therefore changes the sum by whole residues, not by an epsilon, and
`sapiens_peptides_sum > .15` is a hard cut in `select()`. Since blastp makes no
promise about the ordering of equally-scoring alignments and that ordering can
vary with `-num_threads`, `Immunity.ordered_hsps` imposes a total order on HSPs
before they are collected. Without it, `--threads` would be silently
result-affecting.

For the same reason `compare_states.py` compares lists element-wise, in order,
rather than as sets.

## Workflow

Run the reference and the variant with snapshots retained:

```bash
./NERVE.sh -p1 tests/input_data/proteome1.fasta -g n -wd runs/baseline -ks True
```

Then compare, stage by stage, to find the first divergence:

```bash
python3 harness/compare_states.py --dir runs/baseline/snapshots runs/variant/snapshots
```

Narrow to one stage and one attribute once you know where to look:

```bash
python3 harness/compare_states.py runs/baseline/snapshots/03_adhesin.pkl \
        runs/variant/snapshots/03_adhesin.pkl --attrs p_ad
```

And check whether the difference reaches the output:

```bash
python3 harness/compare_candidates.py runs/baseline/vaccine_candidates.csv \
                                      runs/variant/vaccine_candidates.csv
```

## Stage ordinals

| # | stage | writes |
|---|-------|--------|
| 1 | quality_control | id, original_id, accession, sequence, length |
| 2 | localization | localization |
| 3 | adhesin | model_raw_data, p_ad |
| 4 | topology | tmhmm_seq, transmembrane_doms |
| 5 | razor | sequence, original_sequence_if_razor, sequence_out |
| 6 | autoimmunity | list_of_shared_human_peps, sapiens_peptides_sum, ...mhcpep_sapiens |
| 7 | mouse | list_of_shared_mouse_peps, mouse_peptides_sum, ...mhcpep_mouse |
| 8 | conservation | list_of_shared_conserv_proteome_peps, conservation_score |
| 9 | virulence | p_vir |
| 10 | annotation | annotations |
| 11 | select | (the selected subset) |
| 12 | epitope | MHC1_binders, MHC2_binders, MHC1_pb_binders, MHC2_pb_binders |

Optional stages are absent from a run that skipped them, and
`compare_states.py` reports a stage present on only one side rather than
silently ignoring it.

## Unit-level equivalence

`tests/test_equivalence.py` and `tests/test_harness.py` assert the parts that
can be checked without running the pipeline, notably that the k-mer index in
`Immunity` reproduces the exhaustive epitope scan exactly, in order, against
the real mhcpep database. They need the pipeline's Python dependencies. Run
them in the built container with the source and harness mounted, or in a local
environment with those dependencies installed:

```bash
pytest tests/test_equivalence.py tests/test_harness.py -p no:debugging
```

(`-p no:debugging` is needed because the repository's package is named `code`,
which shadows the standard-library module of that name that pytest's debugging
plugin imports.)
