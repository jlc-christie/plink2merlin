import os
import sys
import copy
import enum
import time
import uuid
import argparse
import subprocess

from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple

class PlinkWrapper:
    class InputType(enum.Enum):
        BED = 1
        PED = 2
        VCF = 3

        def get_plink_flag(self):
            if self.value == 1:
                return '--bfile'
            elif self.value == 2:
                return '--file'
            elif self.value == 3:
                return '--vcf'

    def __init__(self, args, uuid=str(uuid.uuid4())):
        self.args = args
        self.uuid = uuid # Not currently used, can be used for temp files
        self.input_set = False
        self._validate_input_path()

    def _validate_input_path(self):
        def check_input_init():
            if self.input_set:
                print("Maximum of 1 input allowed, found input:")
                print("{}: {}".format(self.input_type, self.input_str))
                print("exiting...")
                sys.exit()

        if self.args.bed is not None:
            check_input_init()
            self.input_type = self.InputType.BED
            self.input_str = self.args.bed
            self.input_set = True
        elif self.args.ped is not None:
            check_input_init()
            self.input_type = self.InputType.PED
            self.input_str = self.args.ped
            self.input_set = True
        elif self.args.vcf is not None:
            check_input_init()
            self.input_type = self.InputType.VCF
            self.input_str = self.args.vcf
            self.input_set = True
        else:
            if not self.input_set:
                print("No valid input file path has been given")
                print("exiting...")
                sys.exit()

    def run(self, options):
        plink_options = [
            self.args.plink_binary,
            self.input_type.get_plink_flag(), self.input_str,
        ]
        plink_options += options
        subprocess.run(plink_options, capture_output=True)


def inout(f):
    def in_out(*args, **kw):
        start = time.time()
        a = '...'
        print(f"Entering {f.__name__}{a:25}", end='')
        res = f(*args, **kw)
        end = time.time()
        print(f"Exiting. (Finished in {end-start:2.4} seconds)")
        
        return res
    return in_out

class Sex(enum.Enum):
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2

@dataclass
class MerlinRecord:
    fid: str
    iid: int
    mid: int
    pid: int
    sex: Sex
    genotypes: List[str]

    def as_string(self):
        line = f"{self.fid}\t{self.iid}\t{self.mid}\t{self.pid}\t{self.sex.value}\t"
        for i in range(0, len(self.genotypes) - 2, 2):
            line += f"{self.genotypes[i]}/{self.genotypes[i+1]}\t"
        line += f"{self.genotypes[-2]}/{self.genotypes[-1]}\n"
        
        return line

@dataclass
class Individual:
    iid: int
    fid: str
    mid: int
    pid: int
    sex: Sex

    def __eq__(self, obj):
        iid_eq = self.iid == obj.iid
        fid_eq = self.fid == obj.fid
        return iid_eq and fid_eq

    def __hash__(self):
        return hash((self.iid, self.fid))

@dataclass
class Variant:
    chrom: int
    rsid: str
    pos: float
    bp: int

    def __eq__(self, obj):
        return self.rsid == obj.rsid # sketchy

    def __hash__(self):
        return hash(self.rsid)

def get_individuals_from_fam(filename: str) -> Set[Individual]:
    individual_list = []
    with open(filename, 'r') as f:
        for line in f:
            fid, iid, pid, mid, sex, _ = line.split()
            indiv = Individual(int(iid), fid, int(mid), int(pid), Sex(int(sex)))
            individual_list.append(indiv)
    
    return set(individual_list)

def add_missing_individuals(indivs: Set[Individual]):
    updated_indivs = []
    for indiv in indivs:
        if indiv.mid == 0 and indiv.pid == 0:
            updated_indivs.append(indiv)
            continue
        tmp_mother = Individual(indiv.mid, indiv.fid, 0, 0, Sex(2))
        tmp_father = Individual(indiv.pid, indiv.fid, 0, 0, Sex(1))
        if tmp_mother not in indivs:
            updated_indivs.append(tmp_mother)
        if tmp_father not in indivs:
            updated_indivs.append(tmp_father)
        updated_indivs.append(indiv)

    return set(updated_indivs)

def generate_family_map(individuals: Set[Individual]) -> Dict[str, List[Individual]]:
    fam_map = {}
    for indiv in individuals:
        if indiv.fid not in fam_map.keys():
            fam_map[indiv.fid] = []
        fam_map[indiv.fid].append(indiv)

    return fam_map

def filter_useful_fams(fam_map: Dict[str, List[Individual]]) -> Dict[str, List[Individual]]:
    filtered = {}
    for fid, indivs in fam_map.items():
        founders = []
        non_founders = []
        for indiv in indivs:
            if indiv.pid == 0 and indiv.mid == 0:
                founders.append(indiv)
            else:
                non_founders.append(indiv)

        if len(non_founders) > 1:
            filtered[fid] = indivs

    return filtered

def find_disjoint_fams(fam_map: Dict[str, List[Individual]]) -> Dict[str, List[set]]:
    fid_disjoint_map = {}
    for fid, indivs in fam_map.items():
        family_sets = []
        for indiv in indivs:
            tmp_fs = set([indiv.iid])
            if indiv.mid != 0:
                tmp_fs.add(indiv.mid)
            if indiv.pid != 0:
                tmp_fs.add(indiv.pid)

            if len(family_sets) == 0:
                family_sets.append(tmp_fs)
                continue

            for fs in family_sets:
                if not fs.isdisjoint(tmp_fs):
                    tmp_fs |= fs
                    break

            family_sets.append(tmp_fs)

        did_merge = True
        while did_merge:
            did_merge = False
            for i in range(len(family_sets)):
                fs = family_sets[i]
                for j in range(i, len(family_sets)):
                    if i == j:
                        continue
                    tmp = family_sets[j]
                    if not fs.isdisjoint(tmp):
                        family_sets.remove(tmp)
                        family_sets.remove(fs)
                        fs |= tmp
                        family_sets.append(fs)
                        did_merge = True
                        break
                if did_merge:
                    break

        fid_disjoint_map[fid] = family_sets

    return fid_disjoint_map

# Rename, and change func sig
def get_fidpid_genotype_map(plink: PlinkWrapper, fm: Dict[str, List[set]]):
    '''
    This is one of the two big bottlenecks in the script, uses plink to generate
    the .ped/.map files from the .bed files. When the genotypes could be
    directly read from the binary, but from experience doing this, it's very
    easy to make mistakes and not very easy to know IF you've made a mistake.
    Because of this, this just reads the plaintext genotypes from the .ped
    file. 
    '''
    files_to_delete = []
    indiv_tuples = [(indiv.fid, indiv.iid) for indivs in fm.values() for indiv in indivs]
    with open('keep.txt', 'w+') as f:
        for fid, pid in indiv_tuples:
            f.write(f"{fid} {pid}\n")
    files_to_delete.append('keep.txt')

    plink.run([
        '--keep', 'keep.txt',
        '--maf', '0.2',
        '--indep-pairwise', '50', '5', '0.05',
        '--out', 'plink',
    ])
    files_to_delete += ['plink.prune.in', 'plink.prune.out', 'plink.log']
    plink.run([
        '--extract', 'plink.prune.in',
        '--recode',
        '--out', 'pruned'
    ])
    files_to_delete += ['pruned.ped', 'pruned.map', 'pruned.log']
    
    cm_rsid_set = set()
    variants = []
    indices = []
    index = 0
    with open('pruned.map', 'r') as f:
        for line in f:
            chrom, rsid, pos_cm, pos_bp = line.split()
            if pos_cm in cm_rsid_set:
                indices.append(index)
                continue
            cm_rsid_set.add(pos_cm)
            variants.append(Variant(int(chrom), rsid, float(pos_cm), int(pos_bp)))
            index += 1

    fp_geno_map = {}
    with open('pruned.ped', 'r') as f:
        for line in f:
            data = line.split()
            fid, iid, pid, mid, sex = data[:5]
            genotypes = data[6:]
            for i in reversed(indices):
                del genotypes[i:i+2]
            fp_geno_map[fid+" "+iid] = genotypes

    [os.remove(fn) for fn in files_to_delete]
    return fp_geno_map, variants

def create_merlin_records(fam_map, disjoint_fams, fp_geno_map):
    def swap_fid(disjoint_fams, fid, pid, iid):
        tmp_pid = pid
        if pid is 0:
            tmp_pid = iid # Sets founders pid to itself, to get correct fid

        for i in range(len(disjoint_fams[fid])):
            disjoint_fam = disjoint_fams[fid][i]
            if tmp_pid in disjoint_fam:
                return f"{fid}_{i+1}"
    
    records = []
    genos_len = len(list(fp_geno_map.values())[0])
    for fid, indivs in fam_map.items():
        for indiv in indivs:
            new_fid = swap_fid(disjoint_fams, fid, indiv.pid, indiv.iid)
            genotypes = None
            try:
                genotypes = fp_geno_map[f"{fid} {indiv.iid}"]
            except KeyError:
                genotypes = ['0' for _ in range(genos_len)]

            r = MerlinRecord(new_fid, indiv.iid, indiv.mid, indiv.pid,
                             indiv.sex, genotypes)
            records.append(r)

    return records

def write_merlin_ped(records, out_filename):
    '''
    This is the biggest bottleneck, for obvious reasons. Maybe optimising the
    as_string method of the MerlinRecord class would speed it up, not sure how
    much though. 
    '''
    with open(out_filename, 'w+') as f:
        for record in records:
            f.write(record.as_string()) 

def write_merlin_dat(variants, out_filename):
    with open(out_filename, 'w+') as f:
        for variant in variants:
              f.write(f"M {variant.rsid}\n")

def write_merlin_map(variants: List[Variant], outfile: str):
    with open(outfile, 'w+') as f:
        f.write("CHROMOSOME\tMARKER\tPOSITION\n")
        for variant in variants:
            line = f"{variant.chrom}\t{variant.rsid}\t{variant.pos}\n"
            f.write(line)

@inout
def process_chrom(chr_n: int, plink: PlinkWrapper,
                  indir='split_by_chromosome', outdir='merlin_input_files'):
    ''' 
    This function is not ideal, ideally I would have written the
    PlinkWrapper object with a better constructor so that I could easily make a
    new PlinkWrapper object without using output from argparse. Instead I've chosen
    to just deepcopy the object and then manually change the input string.
    
    This function is also very redundant, the family information never changes,
    only the genotypes...
    '''
    pc = copy.deepcopy(plink)
    pc.input_str = f"{indir}/{chr_n}"
    indivs = get_individuals_from_fam(f"{indir}/{chr_n}.fam")
    indivs = add_missing_individuals(indivs)
    fam_map = generate_family_map(indivs)
    # Uncomment line below to only include families with >1 offspring (i.e.
    #   families with sib-pairs)
    #fam_map = filter_useful_fams(fam_map)
    disjoint_fams = find_disjoint_fams(fam_map)
    fp_geno_map, variants = get_fidpid_genotype_map(pc, fam_map)
    records = create_merlin_records(fam_map, disjoint_fams, fp_geno_map)
    fp_geno_map = None # Makes sure we're not holding excessive memory
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    write_merlin_ped(records, f"{outdir}/{chr_n}.ped")
    write_merlin_map(variants, f"{outdir}/{chr_n}.map")
    write_merlin_dat(variants, f"{outdir}/{chr_n}.dat")

def plink_split_by_chrom(plink: PlinkWrapper, outdir='split_by_chromosome'):
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    for chrom in range(1,23):
        plink.run([
            '--chr', str(chrom),
            '--make-bed',
            '--out', f"{outdir}/{chrom}",
        ])

if __name__ == '__main__':
    desc = '''
           File format converter to take a plink binary file or vcf and convert it
           to a format which MERLIN accepts, which consists of a .ped and a .dat
           file. It is important to note that the .ped file here is not the same as
           the plink .ped file.
           '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--bed', type=str, help='PLINK binary input file path stub')
    parser.add_argument('--ped', type=str, help='PLINK input file path stub')
    parser.add_argument('--vcf', type=str, help='VCF file path')
    parser.add_argument('--plink-binary', type=str, default='plink', help='''Path
                        to PLINK binary, useful if plink isn't globally accessible
                        or you want to use a specific version of plink''')
    args = parser.parse_args()
    plink = PlinkWrapper(args)

    plink_split_by_chrom(plink)
    for i in range(1, 23):
        print(f"Processing chromosome {i}...")
        process_chrom(i, plink)
