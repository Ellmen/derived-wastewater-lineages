import json
import json5

import numpy as np


data_dir = 'outbreak_lineages'
lins = ['BA.2', 'BA.1.1', 'B.1.617.2']
# lins = ['B.1.617.2', 'BA.2', 'BA.1.1']

def get_mutation_jsons():
    lins_data = []

    for lin in lins:
        print(lin)
        with open('{}/{}.json'.format(data_dir, lin), 'r') as f:
            js = f.read()
            lins_data.append(json5.loads(js))

    for i in range(len(lins)):
        lin = lins[i]
        lin_data = lins_data[i]
        with open('{}/{}_muts.json'.format(data_dir, lin), 'w') as f:
            f.write(json.dumps({r['mutation']: r['prevalence'] for r in lin_data['results'][lin]}))

rename = {
    # "aa:S:DEL25/27": 'del:21633:9/-21633.9@21632',
    # "aa:S:DEL69/70": 'del:21765:6/-21765.6@21764',
    # "aa:S:DEL143/145": 'del:21987:9/-21987.9@21986',
    # # "aa:S:DEL157/158": 'del:22194:3/-22194.3@22193',
    # "aa:S:DEL212/212": 'del:22194:3/-22194.3@22193',
    "aa:S:DEL25/27": 'del:21633:9',
    "aa:S:DEL69/70": 'del:21765:6',
    "aa:S:DEL143/145": 'del:21987:9',
    # "aa:S:DEL157/158": 'del:22194:3/-22194.3@22193',
    "aa:S:DEL212/212": 'del:22194:3',
    "aa:N:DEL31/33": "del:28362:9",
}

def rename_mut(mut):
    new_mut = 'aa:{}'.format(mut.upper())
    if new_mut in rename:
        new_mut = rename[new_mut]
    return new_mut


def get_muts(gene='S'):
    all_muts = set()
    lins_muts = []
    for lin in lins:
        with open('{}/{}_muts.json'.format(data_dir, lin), 'r') as f:
            old_muts = json.loads(f.read())
            muts = {rename_mut(mut): old_muts[mut] for mut in old_muts}
            lins_muts.append(muts)
            for mut in muts:
                if muts[mut] > 0.25 and (mut.startswith('aa:{}'.format(gene)) or mut in rename.values()):
                    all_muts.add(mut)
    all_muts = sorted(list(all_muts))
    # print(all_muts)
    mut_data = [[lin_muts[mut] if mut in lin_muts else 0 for mut in all_muts] for lin_muts in lins_muts]
    return np.array(all_muts), np.array(mut_data)

