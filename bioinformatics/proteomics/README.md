# proteomics

A `uv` project — run `uv sync` to build the environment from `pyproject.toml` / `uv.lock`.

## Sub-projects
- `koina/` — Koina spectral-prediction models (fragment intensities, retention time, detectability).
- `claire/` — Claire's peptide / protein-group occupancy helper (`code.py` + notebook; the notebook does `from code import *`, so they stay together).
- `Zhuoxin/` — Zhuoxin's psmsToOccupancies project.
- `TestRunShortreedMSLReader/` — reader for Shortreed-format MSL spectral libraries, with tests.

## Top-level notebooks
- `MinSetTrypticPeptideMostHumanProteomeCoverage.ipynb` — minimal tryptic-peptide set covering the most human proteome (set-cover).
- `detectability_test.ipynb` — peptide detectability checks.
- `peptide_search_mod_stoich.ipynb` / `-v2.ipynb` — peptide search with modification stoichiometry.
- `JadeOccupancySampleDataset.ipynb` — Jade's occupancy sample dataset.
- `10DOT1126SLASHscisignalDOT2000475.ipynb` — SILAC phospho-occupancy derivation (sympy), after Science Signaling scisignal.2000475.
- `RandomMSPlay.ipynb`, `pyteomics_play.ipynb` — scratch and library exploration.

## Data
Local peptide / proteome lists (`unique_uniprot_human_*.txt`, `peptides_filtered_*.txt`, `synthetic_peptides.tsv`) are untracked via `.gitignore` (`*.txt` / `*.tsv`) — point scripts at them rather than committing them.
