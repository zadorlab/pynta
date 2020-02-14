from itertools import product

import numpy as np

from networkx import is_isomorphic

from .graph_utils import node_test

def find_equivalent_atoms(graph):
    """Given a CatKit Gratoms object, find all groups of equivalent atoms."""
    equiv = dict()
    for idx in graph:
        testgraph1 = graph.copy()
        testgraph1.nodes[idx]['number'] *= -1
        for testgraph2, subgroup in equiv.items():
            if is_isomorphic(testgraph1, testgraph2, node_test):
                equiv[testgraph2].append(idx)
                break
        else:
            equiv[testgraph1] = [idx]
    return list(equiv.values())
    print(graph)
    # return print('aa')


def find_equivalent_permutations(gratoms, indices):
    """Given a CatKit Gratoms object and a set of indices corresopnding
    to a reaction, find all equivalent sets of indices based on equivalency
    of the nodes"""

    natoms = len(gratoms)
    # print(natoms)
    groups = find_equivalent_atoms(gratoms.graph)
    groupids = -np.ones(natoms, dtype=int)
    for i, group in enumerate(groups):
        for idx in group:
            groupids[idx] = i

    assert not np.any(groupids == -1)

    basegraph = gratoms.graph.copy()

    for i, idxi in enumerate(indices):
        basegraph.nodes[idxi]['number'] = -i

    permgroups = [groups[groupids[idx]] for idx in indices]
    valid_permutations = []
    for perm_indices in product(*permgroups):
        print(perm_indices)
        testgraph = gratoms.graph.copy()
        for i, idxi in enumerate(perm_indices):
            testgraph.nodes[idxi]['number'] = -i
        if is_isomorphic(basegraph, testgraph, node_test):
            valid_permutations.append(perm_indices)
    return valid_permutations
