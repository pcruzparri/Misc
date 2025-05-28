

__all__ = ['get_mods', 
           'get_base_sequence', 
           'nrmse', 
           'SpectraFileInfo', 
           'Peptide', 
           'BasePeptide', 
           'ProteinGroup', 
           'first']

from collections import defaultdict
import regex as re
import toolz
from dataclasses import dataclass, asdict
from functools import reduce
import random
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import typing
import numpy as np
from Bio.SeqIO.UniprotIO import UniprotIterator
import json

def get_mods(fullSeq: str):
    mod_pattern = r"-?\[.+?\](?<!\[I+\])"
    mods = re.findall(mod_pattern, fullSeq)
    splits = re.split(mod_pattern, fullSeq)
    #print(mods, splits)
    
    m = {}
    running_length = 0
    for i, mod in enumerate(mods):
        running_length += len(splits[i])
        if mod.startswith('-'): 
            running_length += 1
        m[running_length] = re.sub(r'-?\[|\]', '', mod)
        
    return m

def get_base_sequence(full_seq):
    mod_pattern = r"-?\[.+?\](?<!\[I+\])"
    return re.sub(mod_pattern, '', full_seq)

def nrmse(ideal: np.array, observed: np.array):
    return np.sum( np.square(ideal-observed) ) / np.sum( np.square(observed) )

@dataclass
class SpectraFileInfo:
    filename: str
    condition: int
    biorep: int
    techrep: int
    fraction: int
    cond: str
    digest: str


class Peptide:
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
        self.group_keys_labels = group_keys
        files = [asdict(f) for f in spectrafiles.values()]
        reduced_group_sums = toolz.reduceby(group_keys, 
                                            lambda acc, x: [acc[0]+1, 
                                                            acc[1]+self.intensities_per_file[x['filename']]],
                                            files, [0, 0])
        self.reduced_intensities = toolz.valmap(lambda x: x[1]/x[0] if x[0]!=0 else 0, reduced_group_sums) #prevent div by 0 errors
        

    def normalize_file_intensities(self, norms: dict[str, float]):
        self.normalized_intensities_per_file = dict()
        for f in self.spectrafiles:
            self.normalized_intensities_per_file[f] = self.intensities_per_file[f]*norms[f]

    def __repr__(self) -> str:
        return "PeptideObj: " + self.sequence


class BasePeptide: 
    def __init__(self, sequence:str, peptides: list[Peptide], protein_groups: str, prot_dict: dict, file_group_keys: list, file_group: list, files: dict):
        self.modified_peptides = peptides
        self.protein_groups = protein_groups
        self.sequence = sequence
        self.file_group_keys = file_group_keys
        self.file_group = file_group
        self.grouped_spectrafiles = files        
        self.offset_in_proteins = {} # prot_name: pos
        self.protein_seqs = {}
        self._get_prot_pos_offset_from_prot_dict(prot_dict)
        
        """
        structure of self.mods
        {ambiguous_pos: 
            {mod_name: intensity}
        }
        """
        self.mods = {}
        self.total_intensity = 0.0
        self._set_mods_from_peptides(file_group)
        self.ambiguous_modpos = {}
        self._set_ambiguous_modpos()
        self.output_header = '\t'.join(["Base Sequence",
                                        "Full Sequences", 
                                        "Protein Groups", 
                                        "Start Indices in Protein Group",
                                        "Occupancy", 
                                        "Combined Intensity",
                                        "Condition",
                                        "Biorep",
                                        "Techrep", 
                                        "Combined Spectra Files"])

    def _get_prot_pos_offset_from_prot_dict(self, prot_dict):
        for prot in re.split(r';|\|', self.protein_groups):
            self.offset_in_proteins[prot] = prot_dict[prot]['seq'].index(self.sequence) + 1
            self.protein_seqs[prot] = prot_dict[prot]['seq']
        
    def _set_mods_from_peptides(self, grouped_spectrafiles):
        self.file_group = grouped_spectrafiles
        
        for peptide in self.modified_peptides:
            peptide_mods = get_mods(peptide.sequence)
            peptide_intensity = sum([peptide.intensities_per_file[f] for f in self.grouped_spectrafiles])
            self.total_intensity += peptide_intensity
            
            for modpos, mod in peptide_mods.items():
                if modpos not in self.mods.keys():
                    self.mods[modpos] = {}
                if mod not in self.mods[modpos].keys():
                    self.mods[modpos][mod] = 0
                self.mods[modpos][mod]+=peptide_intensity

    def _set_ambiguous_modpos(self):
        for modpos in self.mods:
            ambiguous_pos_string = self.protein_groups    
            
            for prot in self.offset_in_proteins.keys():
                if modpos==0 or modpos==len(self.sequence)+1: #if terminus mod
                    if modpos == modpos + self.offset_in_proteins[prot] == 0:
                        ambiguous_pos_string = ambiguous_pos_string.replace(prot, 'N-term')
                    elif modpos + self.offset_in_proteins[prot] + 1 == len(self.protein_seqs[prot]) + 1:
                        ambiguous_pos_string = ambiguous_pos_string.replace(prot, 'C-term')
                    else:
                        ambiguous_pos_string = ambiguous_pos_string.replace(prot, 'X') # terminal in peptide but not in protein.
                else:
                    protein_modpos = modpos + self.offset_in_proteins[prot] - 1
                    repl = self.protein_seqs[prot][ protein_modpos - 1 ] + str(protein_modpos)
                    ambiguous_pos_string = ambiguous_pos_string.replace(prot, repl)
            self.ambiguous_modpos[modpos] = ambiguous_pos_string    
            

    def get_occupancy_string(self):
        #self.occupancy: dict[] = {} # {modpos: {modname: {occupancy: #.##, total_intensity: ####}}}
        aa_info = []
        for modpos in sorted(self.mods):
            modpos_string = self.ambiguous_modpos[modpos]+'['
            mod_occupancies = []
            for mod in self.mods[modpos]:
                if self.mods[modpos][mod]>0:
                    mod_occupancies.append(f"{mod}({self.mods[modpos][mod]/self.total_intensity:.4f})")
            modpos_string+=', '.join(mod_occupancies)
            modpos_string+=']'
            aa_info.append(modpos_string)
        return '; '.join(aa_info)
            
    def write_output(self):
        offsets = self.protein_groups
        for prot in self.offset_in_proteins:
            offsets = offsets.replace(prot, str(self.offset_in_proteins[prot]))
        
        return '\t'.join([self.sequence,
                          ', '.join([p.sequence for p in self.modified_peptides]), 
                          self.protein_groups, 
                          offsets,
                          self.get_occupancy_string(), 
                          str(self.total_intensity),
                          *map(str, self.file_group),
                          ', '.join(self.grouped_spectrafiles)])
        


class ProteinGroup:
    def __init__(self, name: str, peptides: list[Peptide]):
        self.proteingroup_name = name
        self.peptides = peptides
        self.spectrafiles = peptides[0].spectrafiles

    # STILL BROKEN! NEED TO FIGURE OUT THE X^2 NORMALIZATION 
    def normalize_intensities(self, good_thresh=0.1):

        # get mean and relative deviations from the mean per file/experiment
        mds = dict() #mean deviations (deviations from the mean) per f per peptide
        total = len(self.peptides)
        mean_per_files = {f: sum([p.intensities_per_file[f] for p in peptides])/total for f in self.spectrafiles}
        for p in self.peptides:
            mds[p] = {f: mean_per_files[f] - p.intensities_per_file[f] for f in self.spectrafiles}

        # get relative median deviation for the means
        median_deviations_per_file: dict[str, float] = dict()
        for f in spectrafiles:
            median_deviations_per_file[f] = np.median([mds[p][f] for p in self.peptides])

        self.good_bad_peptides = {'good': [], 'bad':[]}
        for f in self.spectrafiles:
            for p in self.peptides:
                if abs(mds[p][f]-median_deviations_per_file[f])/median_deviations_per_file[f] < good_thresh:
                    self.good_bad_per_file[f]['good'].append(p)
                    #print("Good: ", p.intensities_per_file[f])
                else:
                    self.good_bad_per_file[f]['bad'].append(p)
                    #print("Bad: ", p.intensities_per_file[f])
        
        #return good_bad_per_file
        self.norms_per_file = dict()
        for f in self.spectrafiles:
            total_good_intensity = sum([p.intensities_per_file[f] for p in self.good_bad_per_file[f]['good']])/len(self.good_bad_per_file[f]['good'])
            self.norms_per_file[f] = 1/total_good_intensity if total_good_intensity else 0
        print([len(self.good_bad_per_file[f]['good']) for f in self.spectrafiles])
        print(self.norms_per_file)
        for p in self.peptides: 
            p.normalize_file_intensities(self.norms_per_file)

    
    def __repr__(self):
        return 'ProteinGroupObj ' + self.proteingroup_name + f' w/ {len(self.peptides)} peptides.'


def first(d: dict):
    k = next(iter(d))
    return (k, d[k])