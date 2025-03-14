from pysidt.sidt import SubgraphIsomorphicDecisionTree, Rule, Datum
from ase.data import covalent_radii, atomic_numbers
import logging 
import numpy as np

class SurfaceDatum(Datum):
    """
    Datum class that contains metal and facet information
    """
    def __init__(self, mol, value, metal, facet, comment="", weight=1.0, uncertainty=0.0) -> None:
        super().__init__(mol, value, weight=weight, uncertainty=uncertainty)
        self.metal = metal
        self.facet = facet
        self.comment = comment

class SurfaceBondLengthSIDT(SubgraphIsomorphicDecisionTree):
    def fit_tree(self, data=None):
        """
        fit rule for each node
        """
        if data:
            self.clear_data()
            self.root.items = data[:]
            self.descend_training_from_top(only_specific_match=False)

        for node in self.nodes.values():
            if not node.items:
                logging.warning(f"Node: {node.name} was empty")
                node.rule = None 
                continue

            n = len(node.items)
            
            if n == 1:
                node.rule = Rule(value=[0.0,node.items[0].value], uncertainty=None, num_data=1)
            else:
                metals = [d.metal for d in node.items]
                crad = [covalent_radii[atomic_numbers[x]] for x in metals]
                bdlens = [d.value for d in node.items]
                coefs, residuals, rank, singular_values, rcond = np.polyfit(crad,bdlens,deg=1,full=True)
                
                if coefs[0] < 0.0: #bond length can't decrease with larger metal atoms
                    coefs[0] = 0.0
                    coefs[1] = np.mean(bdlens)
                
                data_var = np.var([bdlens[i] -  np.poly1d(coefs)(crad[i]) for i in range(len(bdlens))])
                
                
                node.rule = Rule(value=coefs.tolist(), uncertainty=data_var, num_data=n)
        
        for node in self.nodes.values():
            n = node
            while n.rule is None:
                n = n.parent
            node.rule = n.rule
            if node.rule.uncertainty is None:
                node.rule.uncertainty = node.parent.rule.uncertainty

    def evaluate(self, mol, metal, trace=False, estimate_uncertainty=False):
        """
        Evaluate tree for a given possibly labeled mol
        """
        children = self.root.children
        node = self.root
        crad = covalent_radii[atomic_numbers[metal]]

        while children:
            for child in children:
                if mol.is_subgraph_isomorphic(
                    child.group, generate_initial_map=True, save_order=True
                ):
                    children = child.children
                    node = child
                    break
            else:
                if trace and estimate_uncertainty:
                    return np.poly1d(node.rule.value)(crad), np.sqrt(node.rule.uncertainty), node.name
                elif estimate_uncertainty:
                    return np.poly1d(node.rule.value)(crad), np.sqrt(node.rule.uncertainty)
                elif trace:
                    return np.poly1d(node.rule.value)(crad), node.name
                else:
                    return np.poly1d(node.rule.value)(crad)

        if trace and estimate_uncertainty:
            return np.poly1d(node.rule.value)(crad), np.sqrt(node.rule.uncertainty), node.name
        elif estimate_uncertainty:
            return np.poly1d(node.rule.value)(crad), np.sqrt(node.rule.uncertainty)
        elif trace:
            return np.poly1d(node.rule.value)(crad), node.name
        else:
            return np.poly1d(node.rule.value)(crad)