# bioinformatics

- `proteomics/` — self-contained uv project (`pyproject.toml` + `uv.lock`) and the home for all proteomics work. Top-level notebooks cover occupancy sampling, tryptic-peptide proteome coverage, detectability, mod-stoichiometry search, pyteomics play, and a Science Signaling (scisignal.2000475) occupancy derivation. Sub-projects:
  - `koina/` — Koina spectral-prediction work: model benchmarking, input generation, API test.
  - `claire/` — Claire's peptide / protein-group occupancy helper (`code.py` + its notebook; the notebook does `from code import *`, so they stay together).
  - `Zhuoxin/` — Zhuoxin's psmsToOccupancies project (kept separate from claire's).
  - `TestRunShortreedMSLReader/` — MSL spectral-library reader with tests.
- `notebooks/` — non-proteomics bioinformatics: affine-gap alignment viz and AAindex play.
- `data/` — GENCODE v48 transcripts / translations, AAindex, and index lists. Kept on disk but untracked; point scripts here rather than baking in absolute paths.
