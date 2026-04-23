

__all__ = [
    'get_mods',
    'get_base_sequence',
    'nrmse',
    'SpectraFileInfo',
    'Peptide',
    'BasePeptide',
    'ProteinGroup',
    'first',
]

import regex as re
import toolz
import numpy as np
from dataclasses import dataclass, asdict
from Bio.SeqIO.UniprotIO import UniprotIterator


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def get_mods(fullSeq: str) -> dict[int, str]:
    """Parse modification annotations from a full peptide sequence string.

    Returns a dict mapping residue position (1-based) to modification name.
    Modifications are encoded as bracketed strings, e.g. 'PEPTM[Oxidation]IDE'.
    Leading '-' indicates a terminus modification.
    """
    mod_pattern = r"-?\[.+?\](?<!\[I+\])"
    mods = re.findall(mod_pattern, fullSeq)
    splits = re.split(mod_pattern, fullSeq)

    posmods = {}
    running_length = 0
    for i, mod in enumerate(mods):
        running_length += len(splits[i])
        if mod.startswith('-'):
            running_length += 1
        posmods[running_length] = re.sub(r'-?\[|\]', '', mod)

    return posmods


def get_base_sequence(full_seq: str) -> str:
    """Strip all modification annotations from a peptide sequence string."""
    mod_pattern = r"-?\[.+?\](?<!\[I+\])"
    return re.sub(mod_pattern, '', full_seq)


def nrmse(ideal: np.ndarray, observed: np.ndarray) -> float:
    """Normalised root-mean-square error between ideal and observed arrays."""
    return np.sum(np.square(ideal - observed)) / np.sum(np.square(observed))


def first(d: dict):
    """Return the first (key, value) pair from a dict."""
    k = next(iter(d))
    return (k, d[k])


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SpectraFileInfo:
    """Metadata associated with a single spectra (raw) file."""
    filename: str
    condition: int
    biorep: int
    techrep: int
    fraction: int
    cond: str
    digest: str


# ---------------------------------------------------------------------------
# Core classes
# ---------------------------------------------------------------------------

class Peptide:
    """A modified peptide sequence with per-file intensity and detection data."""

    def __init__(self,
                 sequence: str,
                 base_sequence: str,
                 protein_groups: str,
                 gene_names: str = '',
                 organism: str = '',
                 spectrafiles: dict[str, SpectraFileInfo] = dict(),
                 intensities_per_file: dict[str, float] = dict(),
                 detectiontypes_per_file: dict[str, str] = dict()):

        self.sequence = sequence
        self.base_sequence = base_sequence
        self.protein_groups = protein_groups
        self.gene_names = gene_names
        self.organism = organism
        self.spectrafiles = spectrafiles
        self.intensities_per_file = intensities_per_file
        self.detectiontypes_per_file = detectiontypes_per_file

    def group_spectrafiles(self, group_keys: list) -> None:
        """Aggregate per-file intensities by the given grouping keys.

        Results are stored in self.reduced_intensities as a dict mapping
        each group key tuple to the mean intensity across files in that group.
        """
        self.group_keys_labels = group_keys
        files = [asdict(f) for f in self.spectrafiles.values()]
        reduced_group_sums = toolz.reduceby(
            group_keys,
            lambda acc, x: [acc[0] + 1, acc[1] + self.intensities_per_file[x['filename']]],
            files,
            [0, 0],
        )
        # Divide summed intensity by count; guard against divide-by-zero
        self.reduced_intensities = toolz.valmap(
            lambda x: x[1] / x[0] if x[0] != 0 else 0,
            reduced_group_sums,
        )

    def normalize_file_intensities(self, norms: dict[str, float]) -> None:
        """Scale per-file intensities by the provided per-file normalisation factors."""
        self.normalized_intensities_per_file = {
            f: self.intensities_per_file[f] * norms[f]
            for f in self.spectrafiles
        }

    def __repr__(self) -> str:
        return "PeptideObj: " + self.sequence


class BasePeptide:
    """Groups all modified forms of a single base (unmodified) peptide sequence.

    Aggregates modification intensities across all modified peptides and
    computes per-residue modification occupancies relative to the total
    peptide intensity.

    self.mods structure:
        {position: {mod_name: cumulative_intensity}}
    """

    # Tab-separated column headers for write_output()
    output_header = '\t'.join([
        "Base Sequence",
        "Full Sequences",
        "Protein Groups",
        "Start Indices in Protein Group",
        "Occupancy",
        "Combined Intensity",
        "Condition",
        "Biorep",
        "Techrep",
        "Combined Spectra Files",
    ])

    def __init__(self,
                 sequence: str,
                 peptides: list[Peptide],
                 protein_groups: str,
                 prot_dict: dict,
                 file_group_keys: list,
                 file_group: list,
                 files: dict):

        self.sequence = sequence
        self.modified_peptides = peptides
        self.protein_groups = protein_groups
        self.file_group_keys = file_group_keys
        self.file_group = file_group
        self.grouped_spectrafiles = files

        # Map each protein accession to the 1-based start position of this peptide
        self.offset_in_proteins: dict[str, int] = {}
        self.protein_seqs: dict[str, str] = {}
        self._get_prot_pos_offset_from_prot_dict(prot_dict)

        self.mods: dict[int, dict[str, float]] = {}
        self.total_intensity: float = 0.0
        self._set_mods_from_peptides(file_group)

        # Map each modified position to a human-readable protein position string
        self.ambiguous_modpos: dict[int, str] = {}
        self._set_ambiguous_modpos()

    def _get_prot_pos_offset_from_prot_dict(self, prot_dict: dict) -> None:
        """Populate offset_in_proteins and protein_seqs from the protein dict."""
        for prot in re.split(r';|\|', self.protein_groups):
            self.offset_in_proteins[prot] = prot_dict[prot]['seq'].index(self.sequence) + 1
            self.protein_seqs[prot] = prot_dict[prot]['seq']

    def _set_mods_from_peptides(self, grouped_spectrafiles) -> None:
        """Accumulate modification intensities from all modified peptide forms."""
        self.file_group = grouped_spectrafiles

        for peptide in self.modified_peptides:
            peptide_mods = get_mods(peptide.sequence)
            peptide_intensity = sum(
                peptide.intensities_per_file[f] for f in self.grouped_spectrafiles
            )
            self.total_intensity += peptide_intensity

            for modpos, mod in peptide_mods.items():
                self.mods.setdefault(modpos, {}).setdefault(mod, 0)
                self.mods[modpos][mod] += peptide_intensity

    def _set_ambiguous_modpos(self) -> None:
        """Build human-readable position strings for each modified residue.

        Terminus modifications are labelled 'N-term' or 'C-term'; internal
        residues are labelled with the amino-acid letter and 1-based protein
        position (e.g. 'S123').
        """
        for modpos in self.mods:
            pos_string = self.protein_groups

            for prot in self.offset_in_proteins:
                if modpos == 0 or modpos == len(self.sequence) + 1:
                    # Terminus modification
                    if modpos == modpos + self.offset_in_proteins[prot] == 0:
                        pos_string = pos_string.replace(prot, 'N-term')
                    elif modpos + self.offset_in_proteins[prot] + 1 == len(self.protein_seqs[prot]) + 1:
                        pos_string = pos_string.replace(prot, 'C-term')
                    else:
                        break
                else:
                    protein_modpos = modpos + self.offset_in_proteins[prot] - 1
                    aa_label = self.protein_seqs[prot][protein_modpos - 1] + str(protein_modpos)
                    pos_string = pos_string.replace(prot, aa_label)

            if pos_string != self.protein_groups:
                self.ambiguous_modpos[modpos] = pos_string

    def get_occupancy_string(self) -> str:
        """Return a formatted string of modification occupancies.

        Format: 'ProtPos[ModName(occupancy)]; ...'
        Occupancy is the fraction of total peptide intensity carrying that mod.
        """
        aa_info = []
        for modpos in sorted(self.ambiguous_modpos):
            mod_occupancies = [
                f"{mod}({self.mods[modpos][mod] / self.total_intensity:.4f})"
                for mod in self.mods[modpos]
                if self.mods[modpos][mod] > 0
            ]
            if mod_occupancies:
                aa_info.append(self.ambiguous_modpos[modpos] + '[' + ', '.join(mod_occupancies) + ']')

        return '; '.join(aa_info)

    def write_output(self) -> str:
        """Return a tab-separated output row for this base peptide."""
        offsets = self.protein_groups
        for prot in self.offset_in_proteins:
            offsets = offsets.replace(prot, str(self.offset_in_proteins[prot]))

        return '\t'.join([
            self.sequence,
            ', '.join(p.sequence for p in self.modified_peptides),
            self.protein_groups,
            offsets,
            self.get_occupancy_string(),
            str(self.total_intensity),
            *map(str, self.file_group),
            ', '.join(self.grouped_spectrafiles),
        ])


class ProteinGroup:
    """A group of peptides belonging to the same protein (group).

    NOTE: normalize_intensities is currently broken — the chi-squared
    normalisation logic still needs to be implemented.
    """

    def __init__(self, name: str, peptides: list[Peptide]):
        self.proteingroup_name = name
        self.peptides = peptides
        self.spectrafiles = peptides[0].spectrafiles

    def normalize_intensities(self, good_thresh: float = 0.1) -> None:
        """Normalise per-file peptide intensities using median-deviation filtering.

        Peptides whose per-file deviation from the median is within good_thresh
        are considered 'good' and used to compute per-file normalisation factors.

        TODO: implement chi-squared normalisation.
        """
        total = len(self.peptides)

        # Mean intensity per file across all peptides
        mean_per_file = {
            f: sum(p.intensities_per_file[f] for p in self.peptides) / total
            for f in self.spectrafiles
        }

        # Absolute deviation of each peptide from the per-file mean
        mds = {
            p: {f: mean_per_file[f] - p.intensities_per_file[f] for f in self.spectrafiles}
            for p in self.peptides
        }

        # Median deviation across all peptides for each file
        median_deviations_per_file: dict[str, float] = {
            f: np.median([mds[p][f] for p in self.peptides])
            for f in self.spectrafiles
        }

        # Classify peptides as 'good' or 'bad' per file based on deviation threshold
        self.good_bad_per_file: dict[str, dict[str, list]] = {
            f: {'good': [], 'bad': []} for f in self.spectrafiles
        }
        for f in self.spectrafiles:
            for p in self.peptides:
                rel_dev = abs(mds[p][f] - median_deviations_per_file[f]) / median_deviations_per_file[f]
                if rel_dev < good_thresh:
                    self.good_bad_per_file[f]['good'].append(p)
                else:
                    self.good_bad_per_file[f]['bad'].append(p)

        # Compute per-file normalisation factors from 'good' peptides
        self.norms_per_file: dict[str, float] = {}
        for f in self.spectrafiles:
            good_peptides = self.good_bad_per_file[f]['good']
            total_good_intensity = (
                sum(p.intensities_per_file[f] for p in good_peptides) / len(good_peptides)
                if good_peptides else 0
            )
            self.norms_per_file[f] = 1 / total_good_intensity if total_good_intensity else 0

        for p in self.peptides:
            p.normalize_file_intensities(self.norms_per_file)

    def __repr__(self) -> str:
        return f'ProteinGroupObj {self.proteingroup_name} w/ {len(self.peptides)} peptides.'
