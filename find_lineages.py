from collections import defaultdict
from os import listdir, mkdir, system
from os.path import isdir, isfile, join

import json

import numpy as np


def find_mutants_in_grom(grom_path, aa_mutations, filter_primers=True):
    with open(grom_path, 'r') as f:
        mut_lines = [line.split(',') for line in f.readlines()[1:]]
    mut_results = defaultdict(float)
    for line in mut_lines:
        pos = int(line[0])
        label = line[1]
        mut = line[2]
        mut_name = '{}/{}@{}'.format(mut, label, pos)
        aa_mutations.add(mut_name)
        freq = float(line[3])
        # if freq < 0.05:
        #     freq = 0
        mut_results[mut_name] = freq

    return mut_results


def mut_idx(mut):
    # Sort by genomic index of mutations
    return int(mut.split('@')[1])


def make_mut_data(virus):
    data_paths = ['./data/{}/runs/{}'.format(virus, f) for f in listdir('./data/{}/runs'.format(virus)) if isdir('./data/{}/runs/{}'.format(virus, f))]
    aa_mutations = set()
    mut_results = []

    for data_path in data_paths:
        print('Getting data in {}'.format(data_path))
        map_csvs = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.mapped.csv')]
        map_csvs.sort()
        for map_csv in map_csvs:
            grom_path = '{}/{}'.format(data_path, map_csv)
            mr = find_mutants_in_grom(grom_path, aa_mutations)
            mut_results.append(mr)

    sorted_muts = sorted(aa_mutations, key=mut_idx)

    mut_data = [[mr[mut] for mut in sorted_muts] for mr in mut_results]

    with open('data/{}/mut_data.npy'.format(virus), 'wb') as f:
        np.save(f, np.array(sorted_muts))
        np.save(f, np.array(mut_data))


def load_muts(virus):
    with open('data/{}/mut_data.npy'.format(virus), 'rb') as f:
        muts = np.load(f)
    return muts


def load_mut_data(virus):
    with open('data/{}/mut_data.npy'.format(virus), 'rb') as f:
        muts = np.load(f)
        mut_data = np.load(f)
    return muts, mut_data


def do_nmf(virus, mut_data, n_components=5):
    from sklearn.decomposition import NMF
    # nmf = NMF(n_components=n_components, init='nndsvd')
    nmf = NMF(n_components=n_components, init='nndsvdar')
    nmf.fit(mut_data)
    with open('data/{}/mut_nmf.npy'.format(virus), 'wb') as f:
        np.save(f, nmf.components_)


def load_nmf(virus):
    with open('data/{}/mut_nmf.npy'.format(virus), 'rb') as f:
        nmf = np.load(f, allow_pickle=True)
    return nmf


def print_big_muts(muts, comp, cutoff=0.25, just_muts=True):
    max_val = max(comp)
    for i in range(len(muts)):
        if comp[i] > max_val * cutoff:
            mut, label, pos = muts[i].replace('/', '@').split('@')
            if just_muts:
                print(mut)
            else:
                print('{}: {}'.format(muts[i], comp[i]/max_val))


def save_seqs(virus, muts, comp, comp_name, cutoff=0.25):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    ref_seq = next(SeqIO.parse('data/{}/reference.fasta'.format(virus), 'fasta'))
    new_seq = [nt for nt in str(ref_seq.seq)]
    max_val = max(comp)
    labels = []
    for i in range(len(muts)):
        if comp[i] > max_val * cutoff:
            mut, label, pos = muts[i].replace('/', '@').split('@')
            labels.append(label)
    for label in labels:
        if label[0] == '~':
            pos = int(label[1:-1])
            nt = label[-1]
            new_seq[pos-1] = nt
        elif label[0] == '-':
            pos, del_l = label[1:].split('.')
            pos = int(pos)
            del_l = int(del_l)
            for p in range(pos, pos+del_l):
                new_seq[p] = ''
    new_seq = ''.join(new_seq)
    new_record = SeqRecord(
        Seq(new_seq),
        id=comp_name,
    )
    f_path = "data/{}/{}.fasta".format(virus, comp_name)
    SeqIO.write([new_record], f_path, "fasta")
    print('Saved to {}'.format(f_path))


def only_orf(virus, orf, muts, mut_data=None):
    if orf == 'all':
        if mut_data is not None:
            return muts, mut_data
        return muts
    self.orfs = {
        'sars-cov-2': {
            'orf1a': (265, 13468),
            'orf1b': (13467, 21555),
            'S': (21562, 25384),
            'orf3a': (25392, 26220),
            'E': (26244, 26472),
            'M': (26522, 27191),
            'orf6': (27201, 27387),
            'orf7a': (27393, 27759),
            'orf7b': (27755, 27887),
            'orf8': (27893, 28259),
            'N': (28273, 29533),
            'orf10': (29557, 29674),
        },
        'tobrfv': {
            'RdRp': (76, 4924),
            'MP': (4910, 5711),
            'CP': (5713, 6193),
        },
    }
    min_idx = False
    max_idx = False
    start, end = orfs[orf]
    for i in range(len(muts)):
        mut = muts[i]
        pos = int(mut.split('@')[1])
        if not min_idx and pos > start:
            min_idx = i
        if pos < end:
            max_idx = i
    if mut_data is None:
        return muts[min_idx:max_idx]
    return muts[min_idx:max_idx], mut_data[:, min_idx:max_idx]


virus = 'sars-cov-2' # sars-cov-2 or tobrfv
orf = 'all' # See orfs in only_orf function
num_lineages = 3

make_mut_data(virus) # Can be commented out after this is run the first time
muts, mut_data = load_mut_data(virus)
muts, mut_data = only_orf(virus, orf, muts, mut_data)
do_nmf(virus, mut_data, n_components=num_lineages) # Can be commented out if no changes are made
muts = load_muts(virus)
muts = only_orf(virus, orf, muts)
nmf = load_nmf(virus)
idx = 0
for comp in nmf:
    idx += 1
    print('Lineage #{}'.format(idx))
    print_big_muts(muts, comp, just_muts=True)
    save_seqs(virus, muts, comp, 'lineage{}'.format(idx))
    print()

