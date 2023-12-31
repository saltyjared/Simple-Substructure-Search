import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from provided import parse_sdf


class Molecule:
    """Molecule object constructed from a .sdf file"""

    def __init__(self, filepath: str):
        """
        Constructs a Molecule object from a .sdf file, complete
        with its own graphical representation in NetworkX

        Parameters
        ----------
        filepath : str
        Filepath of .sdf file containing a molecule of interest

        Returns
        -------
        Molecule object with callable `.graph` attribute
        """
        self.sdf_contents = parse_sdf(filepath, include_hydrogen=False)
        self.element_dict = self.sdf_contents[0]
        self.bond_tuples = self.sdf_contents[1]
        self.graph = nx.Graph()
        self.graph.add_nodes_from((list(self.element_dict.keys())))
        self.graph.add_weighted_edges_from(self.bond_tuples)
        self.edge_labels = nx.get_edge_attributes(self.graph, "weight")

    def draw_graph(self):
        """
        Draws graph of molecule based on its .sdf contents

        Parameters
        ----------
        self : Molecule

        Returns
        -------
        Matplotlib visualization, accessed in _ipython_display_()
        """
        plt.figure(figsize=(10, 10))
        pos = nx.kamada_kawai_layout(self.graph)
        nx.draw_networkx(G=self.graph, pos=pos, edgecolors="black", with_labels=False)
        nx.draw_networkx_labels(G=self.graph, pos=pos, labels=self.element_dict)
        nx.draw_networkx_edge_labels(
            G=self.graph, pos=pos, edge_labels=self.edge_labels
        )

    def _ipython_display_(self):
        """
        Method required to access a graphical representation of a Molecule
        by running an instance as the last line of a cell
        """
        self.draw_graph()

    @property
    def fingerprint(self):
        """Property for accessing the molecular fingerprint of a Molecule"""
        return self.get_fingerprint()

    def get_fingerprint(self):
        """
        Returns a 1024-bit molecular fingerprint of a molecule by converting
        all bond paths up to length of 7 to its SMILES equivalent, hashing
        each path, and assigning 2 bits/path in the fingerprint as true (1)

        Parameters
        ----------
        self : Molecule

        Returns
        -------
        fingerprint : numpy.ndarray
        Array representation of a molecular fingerprint containing 1024
        bits (elements) marked as either 0 or 1
        """
        fingerprint = np.zeros(1024, dtype=int)
        nodes = list(self.graph.nodes)
        smiles_paths = {}

        for start in nodes:
            for end in nodes:
                for path in nx.all_simple_paths(
                    self.graph, source=start, target=end, cutoff=6
                ):
                    bond_path = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
                    smiles_path = ""
                    for node1, node2 in bond_path:
                        bond_weight = int(self.graph[node1][node2]["weight"])
                        bond_symbol = {1: "", 2: "=", 3: "#"}[bond_weight]
                        smiles_path += self.element_dict[node1] + bond_symbol

                    smiles_path += self.element_dict[end]
                    smiles_paths[(start, end)] = smiles_path
        unique_smiles = set(smiles_paths.values())

        for smile in unique_smiles:
            random.seed(hash(str(smile)))
            rand_ints = set(random.sample(range(1024), 2))
            for rand_int in rand_ints:
                fingerprint[rand_int] = 1

        return fingerprint

    def __eq__(self, other: "Molecule"):
        """
        Equivalence operator (=) to compare two Molecule objects based
        on their fingerprints

        Parameters
        ----------
        self : Molecule

        other : Molecule
        Another instance of a Molecule to be compared to, based on fingerprint

        Returns
        -------
        bool
        Returns whether two Molecules are equivalent based on their fingerprint
        (whether their array representations have the same elements in order)
        """
        if not isinstance(other, Molecule):
            raise TypeError("Right-hand side must be a Molecule instance.")

        return (self.fingerprint == other.fingerprint).all()

    def check_substructure_sdf(self, other: "Molecule"):
        """
        Checks if the current Molecule object contains a substructure
        represented by a Molecule object created from a .sdf file

        Parameters
        ----------
        self : Molecule

        other : Molecule
        Molecule object of a potential substructure constructed from a
        .sdf file

        Returns
        -------
        bool
        Returns whether the element-wise difference of the two fingerprints
        results in an array with all elements that are greater than or
        equal to 0
        """
        if not isinstance(other, Molecule):
            raise TypeError("Input must be a Molecule instance.")

        match = self.fingerprint - other.fingerprint
        return np.all(match >= 0)

    def check_substructure_smiles(self, other: str):
        """
        Checks if the current Molecule object contains a substructure
        represented by a SMILES path

        Parameters
        ----------
        self : Molecule

        other : str
        SMILES string path of a potential substructure

        Returns
        -------
        bool
        Returns whether the element-wise difference of the two fingerprints
        results in an array with all elements that are greater than or
        equal to 0
        """
        if not isinstance(other, str):
            raise TypeError("Input must be a SMILES string.")

        other_fp = np.zeros(1024, dtype=int)
        random.seed(hash(str(other)))
        rand_ints = set(random.sample(range(1024), 2))
        for rand_int in rand_ints:
            other_fp[rand_int] = 1

        match = self.fingerprint - other_fp
        return np.all(match >= 0)
