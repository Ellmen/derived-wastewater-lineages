from collections import defaultdict
from os import listdir, mkdir, system
from os.path import isdir, isfile, join

import json

import numpy as np


def find_mutants_in_grom(grom_path, aa_mutations):
    with open(grom_path, 'r') as f:
        mut_lines = [line.split(',') for line in f.readlines()[1:]]
    mut_results = defaultdict(float)
    for line in mut_lines:
        pos = int(line[0])
        label = line[1]
        mut = line[2]
        freq = float(line[3])
        cov = int(line[4])
        mut_name = '{}/{}@{}'.format(mut, label, pos)
        if freq < 0.005: # These are almost always sequencing errors which add dimensions
            continue
        if cov < 20:
            freq = np.nan
        aa_mutations.add(mut_name)
        mut_results[mut_name] = freq

    return mut_results


def get_sample_cov(cov_path):
    with open(cov_path, 'r') as f:
        cov_lines = [int(line.split(',')[1]) for line in f.readlines()[1:]]
    return cov_lines


def mut_idx(mut):
    # Sort by genomic index of mutations
    return int(mut.split('@')[1])


def make_mut_data(virus, only_runs=None):
    from sklearn.impute import KNNImputer
    data_paths = ['./data/{}/runs/{}'.format(virus, f) for f in listdir('./data/{}/runs'.format(virus)) if isdir('./data/{}/runs/{}'.format(virus, f))]
    if only_runs is not None:
        data_paths = [path for path in data_paths if any(run in path for run in only_runs)]
    aa_mutations = set()
    mut_results = []
    sample_covs = []

    for data_path in data_paths:
        print('Getting data in {}'.format(data_path))
        map_csvs = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.mapped.csv')]
        cov_csvs = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.coverage.csv')]
        map_csvs.sort()
        cov_csvs.sort()
        for sample_idx in range(len(map_csvs)):
            map_csv = map_csvs[sample_idx]
            cov_csv = cov_csvs[sample_idx]
            grom_path = '{}/{}'.format(data_path, map_csv)
            cov_path = '{}/{}'.format(data_path, cov_csv)
            mr = find_mutants_in_grom(grom_path, aa_mutations)
            mut_results.append(mr)
            sample_covs.append(get_sample_cov(cov_path))

    sorted_muts = np.array(sorted(aa_mutations, key=mut_idx))

    mut_data = []
    for sample_idx in range(len(mut_results)):
        sample_results = []
        mr = mut_results[sample_idx]
        covs = sample_covs[sample_idx]
        for mut in sorted_muts:
            if mut in mr:
                sample_results.append(mr[mut])
            else:
                pos = int(mut.split('@')[1])
                if covs[pos] < 20:
                    sample_results.append(np.nan)
                else:
                    sample_results.append(0)
        mut_data.append(sample_results)
    mut_data = np.array(mut_data)
    # mut_data = np.array([[mr[mut] for mut in sorted_muts] for mr in mut_results])
    # imputer = KNNImputer(n_neighbors=2)
    imputer = KNNImputer(n_neighbors=5)
    # imputer = IterativeImputer(min_value=0, max_value=1)
    mut_data = imputer.fit_transform(mut_data)
    sorted_muts = imputer.get_feature_names_out(input_features=sorted_muts)

    with open('data/{}/mut_data.npy'.format(virus), 'wb') as f:
        np.save(f, sorted_muts)
        np.save(f, mut_data)


def load_muts(virus):
    with open('data/{}/mut_data.npy'.format(virus), 'rb') as f:
        # muts = np.load(f)
        muts = np.load(f, allow_pickle=True)
    return muts


def load_mut_data(virus):
    with open('data/{}/mut_data.npy'.format(virus), 'rb') as f:
        # muts = np.load(f)
        # mut_data = np.load(f)
        muts = np.load(f, allow_pickle=True)
        mut_data = np.load(f, allow_pickle=True)
    return muts, mut_data


def do_nmf(virus, mut_data, n_components=5, save_model=True):
    from sklearn.decomposition import NMF, DictionaryLearning
    # nmf = NMF(n_components=n_components, init='nndsvd')
    nmf = NMF(n_components=n_components, init='nndsvdar', max_iter=1000)
    nmf.fit(mut_data)
    if save_model:
        with open('data/{}/mut_nmf.npy'.format(virus), 'wb') as f:
            np.save(f, nmf.components_)
    return nmf.reconstruction_err_


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


def save_seqs(virus, muts, comps, fasta_name, cutoff=0.25):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    ref_seq = next(SeqIO.parse('data/{}/reference.fasta'.format(virus), 'fasta'))
    new_seqs = []
    idx = 0
    for comp in comps:
        idx += 1
        comp_name = 'lineage{}'.format(idx)
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
        new_seqs.append(new_record)
    f_path = "data/{}/{}.fasta".format(virus, fasta_name)
    SeqIO.write(new_seqs, f_path, "fasta")
    print('Saved to {}'.format(f_path))


def only_orf(virus, orf, muts, mut_data=None):
    if orf == 'all':
        if mut_data is not None:
            return muts, mut_data
        return muts
    orfs = {
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
    start, end = orfs[virus][orf]
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


def only_samples(mut_data, samples):
    data_paths = ['./data/{}/runs/{}'.format(virus, f) for f in listdir('./data/{}/runs'.format(virus)) if isdir('./data/{}/runs/{}'.format(virus, f))]
    # sample_paths = []
    data_mask = []
    for data_path in data_paths:
        map_csvs = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.mapped.csv')]
        map_csvs.sort()
        for map_csv in map_csvs:
            data_mask.append(any(run in join(data_path, map_csv) for run in samples))
    # data_mask = [any(run in data_path for run in samples) for data_path in data_paths]
    return mut_data[data_mask]


def plot_lineage_mutations(muts, comps, cutoff=0.25, lin_names=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set_theme()
    # Normalize the component values
    for i in range(len(comps)):
        max_val = max(comps[i])
        for j in range(len(comps[i])):
            comps[i][j] = comps[i][j]/max_val
    # muts, comps = only_orf('sars-cov-2', 'S', muts, comps)
    # muts, comps = only_orf('sars-cov-2', 'N', muts, comps)
    mask = (np.max(comps, axis=0) > cutoff) # only show mutations in at least one lineage
    sig_comps = comps.T[mask]
    sig_muts = muts[mask]
    mut_names = [mut.split('/')[0] for mut in sig_muts]
    if lin_names is None:
        lin_names = ['Lineage {}'.format(i+1) for i in range(len(comps))]
    sns.heatmap(
        sig_comps,
        cmap=sns.cm.rocket_r,
        vmin=0,
        vmax=1,
        xticklabels=lin_names,
        yticklabels=mut_names,
    )
    plt.xticks(rotation=0)
    plt.xlabel('Predicted lineage')
    plt.ylabel('Mutation')
    plt.show()


def plot_lineages_with_real(muts, comps, cutoff=0.25):
    from get_outbreak_lineages import get_muts
    real_muts, real_comps = get_muts('N')
    real_muts = list(real_muts)
    comp_cols = []
    # muts, comps = only_orf('sars-cov-2', 'S', muts, comps)
    muts, comps = only_orf('sars-cov-2', 'N', muts, comps)
    # Normalize the component values
    for i in range(len(comps)):
        max_val = max(comps[i])
        for j in range(len(comps[i])):
            comps[i][j] = comps[i][j]/max_val
    mask = (np.max(comps, axis=0) > cutoff) # only show mutations in at least one lineage
    comps = comps[:, mask]
    muts = muts[mask]
    parsed_muts = [mut.split('/')[0] for mut in muts]
    for mut in real_muts:
        if mut not in parsed_muts:
            print(mut)
    for real_comp in real_comps:
        comp_col = [real_comp[real_muts.index(mut)] if mut in real_muts else 0 for mut in parsed_muts]
        comp_cols.append(comp_col)
    new_comps = []
    for i in range(len(comps)):
        new_comps.append(comps[i])
        new_comps.append(comp_cols[i])
    lin_names = ['Lineage 1', 'BA.2', 'Lineage 2', 'BA.1.1', 'Lineage 3', 'B.1.617.2']
    # new_comps = np.concatenate([comps, np.array(comp_cols)])
    new_comps = np.array(new_comps)
    plot_lineage_mutations(muts, new_comps, lin_names=lin_names)


def test_multiple_n(virus, mut_data, max_n, label="synthetic"):
    from matplotlib import pyplot as plt
    errs = []
    for num_lineages in range(max_n):
        err = do_nmf(virus, mut_data, n_components=num_lineages+1, save_model=False)
        errs.append(err)
    plt.plot([i+1 for i in range(max_n)], errs)
    plt.xlabel('Number of components')
    plt.ylabel('Reconstruction error')
    plt.savefig(f'figures/num_lineages_{label}.pdf')


virus = 'sars-cov-2' # sars-cov-2 or tobrfv
# virus = 'tobrfv' # sars-cov-2 or tobrfv
orf = 'all' # See orfs in only_orf function
num_lineages = 5
fasta_name = 'synthetic_lineages'
only_runs = [
    'synthetic',
]

make_mut_data(virus, only_runs) # Can be commented out after this is run the first time
muts, mut_data = load_mut_data(virus)
muts, mut_data = only_orf(virus, orf, muts, mut_data)
print(mut_data.shape)
err = do_nmf(virus, mut_data, n_components=num_lineages+1) # Run for single number of lineages
# test_multiple_n(virus, mut_data, 10)
muts = load_muts(virus)
muts = only_orf(virus, orf, muts)
nmf = load_nmf(virus)
print()
idx = 0
for comp in nmf:
    idx += 1
    print('Lineage #{}'.format(idx))
    print_big_muts(muts, comp, just_muts=False)
    print()

save_seqs(virus, muts, nmf, fasta_name)

plot_lineage_mutations(muts, nmf)
# plot_lineages_with_real(muts, nmf)

