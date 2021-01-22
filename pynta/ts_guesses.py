from pynta.excatkit.molecule import Molecule
from pynta.excatkit.gratoms import Gratoms
from pynta.io import IO
from typing import Dict, List, Tuple


class GeneralTSGuessesGenerator():
    def __init__(self,
                 ts_est: str,
                 rxn: Dict[str, str],
                 rxn_name: str,
                 easier_to_build: str,
                 scfactor: float) -> None:
        ''' Initialize :meth:GeneralTSGuessesGenerator method

        Parameters
        ----------
        ts_est : str
            string representations of the species which is use as a TS guess
            skeleton - the one that is used to construct TS guess,
            e.g. ``'OH'`` for ``'OH_O+H'``
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file
        rxn_name : str
            a reaction name
            e.g ``'OH_O+H'``
        easier_to_build : str
            a keys to ts_guess_skeleton; ``'product'`` or ``'reatant'``
        scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction

        '''
        self.ts_est = ts_est
        self.rxn = rxn
        self.rxn_name = rxn_name
        self.easier_to_build = easier_to_build
        self.scfactor = scfactor
        self.ts_guess_skeleton = self.rxn[self.easier_to_build]
        self.reacting_species_connectivity = self.ts_guess_skeleton.split(
            '\n')

    def decide(
            self,
            reacting_idxs: List[int]) -> Tuple[Gratoms, int]:
        ''' Main method to decide to which family of reactions given
        reaction belongs

        Parameters
        ----------
        reacting_idxs : List[int]
            list with indicies of reacting atoms, in order as they appear in
            :literal:`*.yaml` file but not necessary with the same indicies

        Returns
        -------
        ts_guess_image, s_bonded_idx : Tuple[Gratoms, int]
            ts_geuss_el is Gratoms representation of TS guess, whereas
            s_bonded_idx is the index of adsorbate atom that bonds to surface

        '''

        ts_guess_image, s_bonded_idx = None, None
        how_many_atoms_react = len(reacting_idxs)

        if how_many_atoms_react == 2:
            ts_guess_image, s_bonded_idx = Diatomic(
                self.ts_est,
                self.rxn,
                self.rxn_name,
                self.easier_to_build,
                self.scfactor).get_ts_guess_and_bonded_idx(reacting_idxs)

        elif how_many_atoms_react == 3:
            ts_guess_image, s_bonded_idx = Triatomic(
                self.ts_est,
                self.rxn,
                self.rxn_name,
                self.easier_to_build,
                self.scfactor).get_ts_guess_and_bonded_idx(reacting_idxs)

        return ts_guess_image, s_bonded_idx

    def build_ts_guess(self) -> Gratoms:
        ''' Convert ``'ts_est'`` string into a list of Gratoms objects.

        Numner of elements in the list depends is equal to number of
        distinct topologies available for given ``'ts_est'``.

        For diatomics there will always be one element in the list.
        For other types, more than one structure is possible, e.g.

        >>> ts_est = 'COH'
        >>> print(ts_guess_list)
        ['COH' (sp), 'CHO' (sp2)]

        Returns
        -------
        ts_guess_list : List[Gratoms]
            a list of Gratoms objects with all distinct topologies for a given
            ts_est

        '''
        ts_guess_list = Molecule().molecule(self.ts_est)
        return ts_guess_list

    def get_s_bonded_idx(self) -> int:
        ''' Get index of the atoms that connects adsorbate to the surface

        Returns
        -------
        s_bonded_idx: int
            an index of an adsorbate atom that bonds it to surface

        Raises
        ------
        NotImplementedError
            if more then 1 atoms connects adsorbate to the surface - curently,
            only monodentate adsorbtion is supported.

        '''
        surface_indicies = []
        s_bonded_idxs = []
        for line in self.reacting_species_connectivity:
            if 'X' in line:
                index = line.split()[0]
                surface_indicies.append(index)
        for index in surface_indicies:
            keyphrase = '{' + '{}'.format(index)
            for line in self.reacting_species_connectivity:
                if keyphrase in line:
                    s_bonded_idxs.append(line.split()[0])

        # surface reaction, at least one atom connected to the surface
        if s_bonded_idxs:
            if len(s_bonded_idxs) > 1:
                raise NotImplementedError('Only monodendate type of adsorbtion'
                                          ' is currently supported.')
            else:
                s_bonded_idx = int(s_bonded_idxs[0]) - 1

        # gas phase reaction
        else:
            if self.ts_est == 'CH2O':
                s_bonded_idx = 0
            # if self.ts_est in edge_cases_bonded_dict.keys():
            #     s_bonded_idx = [edge_cases_bonded_dict[self.ts_est]]
            else:
                s_bonded_idx = self.get_the_most_connected_atom()
        return s_bonded_idx

    def get_the_most_connected_atom(self) -> int:
        ''' Return index of an atom that has the highest numner of connections
        to the other atoms.

        Returns
        -------
        max_conected_atom_idx : int
            an index of an atom that has the most connections to other atoms

        '''
        number_of_connections = {}
        n_surf_at_befor_ads = 0

        for line in (self.reacting_species_connectivity):
            if 'X' in line:
                n_surf_at_befor_ads += 1
            else:
                break

        for num, line in enumerate(self.reacting_species_connectivity):
            if 'X' not in line:
                number_of_connections[num] = line.count('{')

        for tmp_idx, n_connect in number_of_connections.items():
            if n_connect == max(number_of_connections.values()):
                max_conected_atom_idx = (tmp_idx - n_surf_at_befor_ads - 1)
                return max_conected_atom_idx

    def get_connectivity(
            self,
            yaml_ref_idx1: int,
            yaml_ref_idx2: int,
            visited_atoms: List[int],
            n_surf_at_befor_ads: int,
            multi_line: int,
            first_line: int) -> List[int]:
        ''' Recursvely retrive information about connectivity of the atom_2

        Parameters
        ----------
        yaml_ref_idx1 : int
            a line in the :literal:`*.yaml` file containing information about
            connectivity of the ``'atom_1'``
        yaml_ref_idx2 : int
            a line in the :literal:`*.yaml` file containing information about
            connectivity of the ``'atom_2'``
        visited_atoms : List[int]
            a list with all atomic indicies (tags - as they appear in the
            :literal:`*.yaml`
            file, starts from 0, ingoring surface atoms and multiplicity line)
            that have been alread visited.
        n_surf_at_befor_ads : int
            number of surface atoms 'X' that are present before any adsorbate
            in the :literal:`*.yaml` file
        multi_line : int
            ``1`` if multiplicity line is present in the :literal:`*.yaml`
            file, ``0`` otherwise
        first_line : int
            ``0`` if multiplicity line is present in the :literal:`*.yaml`
            file, ``1`` otherwise

        Returns
        -------
        visited_atoms: List[int]
            a list with all atomic indicies (tags - as they appear in the
            :literal:`*.yaml` file, starts from 0, ingoring surface atoms and
            multiplicity line) that have been alread visited.

        '''
        ref_line_idx = yaml_ref_idx2 - first_line
        reference_line = self.reacting_species_connectivity[ref_line_idx]

        # retrive connectivity information
        yaml_connected_atom_idxs = [int(item.split(',')[0][1:])
                                    for item in reference_line.split()
                                    if '{' in item]

        visited_atoms.append(yaml_ref_idx1)

        for tag in yaml_connected_atom_idxs:
            tag = tag - n_surf_at_befor_ads - multi_line - first_line
            if tag not in visited_atoms and tag != yaml_ref_idx1:
                self.get_connectivity(
                    tag, yaml_ref_idx2,
                    visited_atoms,
                    n_surf_at_befor_ads,
                    multi_line,
                    first_line)

        return visited_atoms

    def get_connected_atoms_tag_idxs(
            self,
            tag_atom_idx1: int,
            tag_atom_idx2: int) -> List[int]:
        ''' Get a list with all atomic indicies of the atoms and subatoms
        and so that are connected to the atom_2, excluding atom_1.

        Parameters
        ----------
        tag_atom_idx1 : int
            an index of the atom_1 as it appears in yaml file, start from 0,
            ignore surface atoms multiplicity line.
            The same idx is used as tag in ts_guess_image.
        tag_atom_idx2 : int
            an index of the atom_2 as it appears in yaml file, start from 0,
            ignore surface atoms multiplicity line.
            The same idx is used as tag in ts_guess_image.

        Returns
        -------
        tag_indicies : List[int]
            a list with all atomic indicies of atoms and subatoms that are
            connected to the atom_2, excluding atom_1. Indicies are as they
            appear in the yaml file, tart from 0, ignore surface atoms and
            multiplicity line.

        '''
        visited_atoms = []
        n_surf_at_befor_ads = 0
        multi_line = 0

        first_line = 1
        if 'multiplicity' in self.reacting_species_connectivity[0]:
            first_line = 0

        for line in self.reacting_species_connectivity:
            if 'multiplicity' in line:
                multi_line += 1
                continue
            if 'X' in line:
                n_surf_at_befor_ads += 1
            else:
                break

        yaml_ref_idx1 = tag_atom_idx1 + n_surf_at_befor_ads \
            + multi_line + first_line
        yaml_ref_idx2 = tag_atom_idx2 + n_surf_at_befor_ads \
            + multi_line + first_line

        tag_indicies = self.get_connectivity(yaml_ref_idx1,
                                             yaml_ref_idx2,
                                             visited_atoms,
                                             n_surf_at_befor_ads,
                                             multi_line, first_line)

        tag_indicies.remove(tag_atom_idx1)
        return tag_indicies

    def convert_tag_to_correct_idx(
            self,
            ts_guess_image: Gratoms,
            tag_atom_idx1: int,
            tag_atom_idx2: int) -> List[int]:
        ''' Convert tag indicies (:literal:`*.yaml` file is a source) to
        indicies as they appear in the ts_guess_image Gratoms object

        Parameters
        ----------
        ts_guess_image : Gratom
            Gratoms representation of the TS guess
        tag_atom_idx1 : int
            an index of the atom_1 as it appears in yaml file, start from 0,
            ignore surface atoms multiplicity line.
            The same idx is used as tag in ts_guess_image.
        tag_atom_idx2 : int
            an index of the atom_2 as it appears in yaml file, start from 0,
            ignore surface atoms multiplicity line.
            The same idx is used as tag in ts_guess_image.

        Returns
        -------
        connected_atoms_idx : List[int]
            a list with all atomic indicies of the atoms and subatoms
            connected to the atom_2, excluding atom_1, as they appear in the
            ts_guess_image, Gratom object.

        '''
        tag_connected_atoms_idx = self.get_connected_atoms_tag_idxs(
            tag_atom_idx1, tag_atom_idx2)
        connected_atoms_idx = []
        for atom in ts_guess_image:
            if atom.tag in tag_connected_atoms_idx:
                connected_atoms_idx.append(atom.index)
        return connected_atoms_idx


class Diatomic(GeneralTSGuessesGenerator):
    def get_ts_guess_and_bonded_idx(
            self,
            reacting_idxs: List[int]) -> Tuple[Gratoms, int]:
        ''' Get ts_guess (Gratom) and index of atom that connects ts_guess
        to the surface. Currently, only monodentate adsobrtion is supported

        Parameters
        ----------
        reacting_idxs : List[int]
            list with indicies of reacting atoms, in order as they appear in
            .yaml file but not necessary with the same indicies

        Returns
        -------
        ts_guess_image, s_bonded_idx : Tuple[Gratoms, int]
            ts_geuss_el is Gratoms representation of the TS guess, whereas
            s_bonded_idx is the index of adsorbate atom that bonds to surface

        '''
        # get TS guess image
        ts_guess_image = IO.get_TS_guess_image(self.rxn, self.easier_to_build)

        # get index of surface bonded atom
        s_bonded_idx = self.get_s_bonded_idx()

        # get reacting atoms indices, as they appear in the .yaml file
        # surface atoms and multiplicity line are ignored
        tag_react_atom_idx_1, tag_react_atom_idx_2 = reacting_idxs

        # convert tag indices to indices as they appear in the Gratoms object
        react_atom_idx_1 = [
            atom.index for atom in ts_guess_image
            if atom.tag == tag_react_atom_idx_1][0]

        react_atom_idx_2 = [
            atom.index for atom in ts_guess_image
            if atom.tag == tag_react_atom_idx_2][0]

        # get atomic indices of all atoms connected to atom2
        connected_atoms = self.convert_tag_to_correct_idx(ts_guess_image,
                                                          tag_react_atom_idx_1,
                                                          tag_react_atom_idx_2)

        # add atom2 idx to the list of connected_atoms
        if react_atom_idx_2 not in connected_atoms:
            connected_atoms.append(react_atom_idx_2)

        # get lenght of the bond between reacting atoms
        bondlen = ts_guess_image.get_distance(
            react_atom_idx_1, react_atom_idx_2)

        n_total_ads_atoms = len(ts_guess_image)

        # create a better TS_guess for a couple of common edge cases
        if n_total_ads_atoms == 2:
            # proper diatomic
            ts_guess_image.rotate(90, 'y')

        elif n_total_ads_atoms == 3:
            # proper triatomic
            remaining_atom_idx = n_total_ads_atoms - \
                (react_atom_idx_1 + react_atom_idx_2)

            if react_atom_idx_2 != 2:
                # Symetrically not important, but structure looks
                # better visually
                ts_guess_image.rotate(-90, 'z')
            else:
                ts_guess_image.rotate(90, 'z')

            # set angle that puts atom_2 closer to the surface
            ts_guess_image.set_angle(
                remaining_atom_idx,
                react_atom_idx_1,
                react_atom_idx_2,
                30,
                indices=[remaining_atom_idx,
                         react_atom_idx_1,
                         react_atom_idx_2],
                add=True)

        else:
            # continue with defaults
            # TODO for now it should work bu should be improved later
            # pass
            ts_guess_image.rotate(90, 'z')

            # scale the bond distance between reacting atoms
        ts_guess_image.set_distance(react_atom_idx_1, react_atom_idx_2,
                                    bondlen * self.scfactor, fix=0,
                                    indices=connected_atoms)

        return ts_guess_image, s_bonded_idx


class Triatomic(GeneralTSGuessesGenerator):
    def get_ts_guess_and_bonded_idx(self, reacting_idxs):
        # will be done later
        raise NotImplementedError('Only diatomic reactions at this moment')
