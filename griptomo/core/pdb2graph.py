# August George, 2022, PNNL

import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import (
    PDBParser,
    MMCIFParser,
)  # BioPython PDB doc: https://biopython.org/docs/1.75/api/Bio.PDB.html
from scipy.spatial.distance import pdist, squareform
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

# for convenience
import warnings

warnings.filterwarnings("ignore")


def get_hydrophobicity(name, warn=False):
    """
    Gets hydrophobicity based on amino acid name.

    Parameters
    ----------
    name : str
        Amino acid name.
    warn : bool, optional
        If True, print a warning when hydrophobicity can't be determined. Default is False.

    Returns
    -------
    float
        Hydrophobicity value.

    Notes
    -----
    Returns NaN if the amino acid name is not recognized.

    Ref: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
    """
    acid_names = [
        "Ile",
        "Val",
        "Leu",
        "Phe",
        "Cys",
        "Met",
        "Ala",
        "Gly",
        "Thr",
        "Ser",
        "Trp",
        "Tyr",
        "Pro",
        "His",
        "Glu",
        "Gln",
        "Asp",
        "Asn",
        "Lys",
        "Arg",
    ]
    acid_kdH = [
        4.5,
        4.2,
        3.8,
        2.8,
        2.5,
        1.9,
        1.8,
        -0.4,
        -0.7,
        -0.8,
        -0.9,
        -1.3,
        -1.6,
        -3.2,
        -3.5,
        -3.5,
        -3.5,
        -3.5,
        -3.9,
        -4.5,
    ]
    hydro_dict = dict(zip(acid_names, acid_kdH))
    try:
        hydrophobicity = hydro_dict[
            f"{name.title()}"
        ]  # make sure string is in correct format
    except:
        hydrophobicity = np.nan
        if warn == True:
            print(
                f"warning! could not assign hydrophobicity for {name.title()}. setting hydrophobicity=NaN"
            )
    return hydrophobicity


def PDB_to_df(pdb_code, fname, pdbx, offset, CA_only=1):
    """
    Loads a PDB (or PDBx) file and stores the atom coordinates and residue name and number into a dataframe.

    Note: if the PDB file has more than one model, the first model is chosen.

    Parameters
    ----------
    pdb_code : str
        PDB ID / label for the protein of interest.
    fname : str
        Filename for the protein of interest. Can be PDB or PDBx format.
    pdbx : int
        Set to 1 if using the newer PDBx file format.
    offset : int
        Index offset in case the first residue ID in PDB file is not the first physical residue (e.g. PDB starts at 5th residue).
    CA_only : int, optional
        Set to 1 [default] if using only alpha carbons, else all atoms are used.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing every atom's x, y, z coordinates and serial number.
    """
    # pick a file reader based on the filetype. Newer/bigger structures will use PDBx
    if pdbx == 1:
        parser = MMCIFParser()  # for PDBx files (.CIF)
    else:
        parser = PDBParser()  # for PDB files (.PDB)
    structure = parser.get_structure(
        pdb_code, fname
    )  # parse pdb file and store as a PDB structure object

    # makes lists storing the alpha carbon (CA) coordinates and serial number
    atom_coord_list = []
    res_id_list = []
    hydro_list = []

    # loop through PDB structure object to get alpha carbon (or all atom) coordinates and serial number
    for i, model in enumerate(
        structure
    ):  # iterate through each model in the structure, exit loop after first model
        if (
            i > 0
        ):  # NMR files have multiple models, but only one model is used for the graph
            break
        for chain in model:  # iterate through each chain in the model
            for residue in chain:  # iterate through each residue in the chain
                for atom in residue:  # iterate through each atom in residue
                    if CA_only == 1:
                        if (
                            atom.get_name() == "CA"
                        ):  # if atom is an alpha carbon then append its coordinates and serial number
                            atom_coord_list.append(
                                atom.get_vector()[:]
                            )  # convert vector object into list of numbers
                            atom_sn = atom.get_serial_number()

                            # biopython PDB module residue.__repr__() returns a string representing the residue:
                            # e.g. '<Residue THR het=  resseq=1 icode= >'
                            # 'resseq' is the residue serial number
                            # --> split the string until the number after 'resseq=' is all that remains
                            res_repr = residue.__repr__()
                            res_repr_tmp = res_repr.split("resseq=", 1)[1]
                            res_sn = (
                                int(res_repr_tmp.split(" ", 1)[0]) + offset
                            )  # get residue serial number (+ offest)
                            res_name = (
                                residue.get_resname()
                            )  # get first letters of residue name
                            hydro_list.append(
                                get_hydrophobicity(res_name)
                            )  # get hydrophobicity based on residue name

                            # store residue ID: <first letter of residue name><residue serial number + offset>
                            res_id_list.append(f"{res_name}{res_sn}-{atom_sn}")
                    else:  # get all atoms
                        if residue.get_resname() != "HOH":
                            atom_coord_list.append(
                                atom.get_vector()[:]
                            )  # convert vector object into list of numbers
                            atom_sn = atom.get_serial_number()

                            # biopython PDB module residue.__repr__() returns a string representing the residue:
                            # e.g. '<Residue THR het=  resseq=1 icode= >'
                            # 'resseq' is the residue serial number
                            # --> split the string until the number after 'resseq=' is all that remains
                            res_repr = residue.__repr__()
                            res_repr_tmp = res_repr.split("resseq=", 1)[1]
                            res_sn = (
                                int(res_repr_tmp.split(" ", 1)[0]) + offset
                            )  # get residue serial number (+ offest)
                            res_name = (
                                residue.get_resname()
                            )  # get first letters of residue name
                            hydro_list.append(
                                get_hydrophobicity(res_name)
                            )  # get hydrophobicity based on residue name

                            # store residue ID: <first letter of residue name><residue serial number + offset>
                            res_id_list.append(f"{res_name}{res_sn}-{atom_sn}")

    # create dataframes for atom coordinates, serial numbers, and then combine them.
    # there's probably a better way to do this in a single line
    atom_df = pd.DataFrame(atom_coord_list, columns=["x", "y", "z"])
    res_id_df = pd.DataFrame(res_id_list, columns=["residue id"])
    hydro_df = pd.DataFrame(hydro_list, columns=["hydrophobicity"])
    df = pd.concat([atom_df, res_id_df, hydro_df], axis=1)
    return df


def PDB_df_to_G(PDB_df, d_cut=8.0):
    """
    Converts a dataframe containing alpha carbon / atom coordinates (in Angstroms) into a graph, G(V,E).

    Each vertex, V, is an alpha carbon / atom. Two alpha carbons / atoms with a distance (in Angstroms) less than a cutoff, d_cut, are connected by an edge, E.

    Parameters
    ----------
    PDB_df : pandas.DataFrame
        A dataframe containing alpha carbon / atom coordinate columns labeled: 'x', 'y', and 'z'.
    d_cut : float, optional
        Threshold for two alpha carbons / atoms to be connected (in Angstroms) by an edge. Defaults to 8.0.

    Returns
    -------
    networkx.Graph
        Protein structure network graph, G(V,E).
    """
    # create distance matrix where element i,j is the Euclidean distance between vertex i and vertex j
    d_matrix = squareform(pdist(PDB_df[["x", "y", "z"]], "euclid"))
    d_matrix_thresh = np.where(
        d_matrix > d_cut, 0, d_matrix
    )  # if the distance is > d_cut, replace it with 0 (i.e. remove edge)
    G = nx.convert_matrix.from_numpy_array(
        d_matrix_thresh
    )  # convert distance matrix to networkx graph object

    # create mapping to relabel nodes from index (i.e. 0...n-1) to residue id (i.e. str(<resname><res s/n + offset>) for each residue)
    label_mapping = {}
    hydro_mapping = {}
    for node in G.nodes():  # here the node value is an index (0-n-1)
        label_mapping[node] = PDB_df["residue id"][node]
        hydro_mapping[node] = PDB_df["hydrophobicity"][node]

    betweenness = nx.betweenness_centrality(G)
    nx.set_node_attributes(G, betweenness, "betweenness")
    nx.set_node_attributes(G, hydro_mapping, "hydrophobicity")
    H = nx.relabel_nodes(
        G, label_mapping
    )  # update the nodes and store as new graph object
    return H


def save_data(df, G, df_name, G_name):
    """
    Convenience function that stores dataframe as .csv and graph as .gexf file.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to save.
    G : networkx.Graph
        Graph to save.
    df_name : str
        Output filename for dataframe .csv.
    G_name : str
        Output filename for graph .gexf.
    """
    df.to_csv(f"{df_name}.csv")  # write dataframe to .csv
    nx.write_gexf(G, f"{G_name}.gexf")  # write graph to .gexf (graph network XML file)


def save_data_at_this_folder(data_path, df, G, df_name, G_name):
    """
    Convenience function that stores dataframe as .csv and graph as .gexf file.

    Parameters
    ----------
    data_path : str or Path
        Output directory path.
    df : pandas.DataFrame
        Dataframe to save.
    G : networkx.Graph
        Graph to save.
    df_name : str
        Output filename for dataframe .csv.
    G_name : str
        Output filename for graph .gexf.
    """
    df_name = Path(data_path, df_name)
    df.to_csv(f"{df_name}.csv")  # write dataframe to .csv
    G_name = Path(data_path, G_name)
    nx.write_gexf(G, f"{G_name}.gexf")  # write graph to .gexf (graph network XML file)


def plot_coordinates(xyz_data, figsize=5):
    """
    Creates a 3D scatter plot containing the xyz data.

    Parameters
    ----------
    xyz_data : np.ndarray
        A[0] = [x0, y0, z0].
    figsize : int, optional
        Size of figure (figsize x figsize). Defaults to 5.

    Returns
    -------
    matplotlib.pyplot.Figure
        3D scatter plot figure.
    """
    fig = plt.figure(figsize=(figsize, figsize))
    plt.title("x y z coordinates")
    ax = fig.add_subplot(projection="3d")
    ax.scatter(
        xyz_data[:, 0], xyz_data[:, 1], xyz_data[:, 2], c="purple", s=2, alpha=0.3
    )
    return fig


def plot_FA_and_CA_coordinates(FA_xyz, CA_xyz, figsize=5):
    """
    Creates a 3D scatter plot containing CA and FA atom coordinates.

    Parameters
    ----------
    FA_xyz : np.ndarray
        A[0] = [x0, y0, z0] for all atom coordinate data.
    CA_xyz : np.ndarray
        A[0] = [x0, y0, z0] for alpha carbon only coordinate data.
    figsize : int, optional
        Size of figure (figsize x figsize). Defaults to 5.

    Returns
    -------
    matplotlib.pyplot.Figure
        3D scatter plot figure.
    """
    fig = plt.figure(figsize=(figsize, figsize))

    plt.title("x y z coordinates - FA and CA")
    ax = fig.add_subplot(projection="3d")
    ax.scatter(
        FA_xyz[:, 0], FA_xyz[:, 1], FA_xyz[:, 2], c="purple", s=5, alpha=0.3, label="FA"
    )
    ax.scatter(
        CA_xyz[:, 0], CA_xyz[:, 1], CA_xyz[:, 2], c="black", s=15, alpha=0.5, label="CA"
    )
    ax.legend()
    return fig


def main(args):
    """
    Takes a .pdb(x) file, converts it into a graph, and saves the atom coordinates to .csv and graph as .gexf.

    Parameters
    ----------
    args : argparse.Namespace
        - pdb_code : str
            PDB id / protein name.
        - fname : str
            PDB/PDBx filename.
        - d_cut : float
            Alpha Carbon / atom pairwise contact distance cutoff (in Angstroms).
        - o : int
            PDB residue index offset integer. Default is 0.
        - pdbx : int
            Set=1 to use pdbx file parser.
        - CA_only : int
            Set=1 to use only alpha carbons (0 for all atoms
    """
    # get arguments
    pdb_code = args.pdb_code
    fname = args.fname
    d_cut = args.d_cut
    o = args.o
    pdbx = args.pdbx
    CA_only = args.CA_only

    df = PDB_to_df(
        pdb_code, fname, pdbx, o, CA_only
    )  # convert pdb file into dataframe of atom/alpha carbon coordinates and s/n
    G = PDB_df_to_G(df, d_cut)  # convert coordinate dataframe into network graph
    save_data(
        df, G, pdb_code, pdb_code
    )  # save dataframe as .csv and graph as .gexf using the same name as the pdb code


if __name__ == "__main__":
    # example: >> python pdb2graph.py PDBid pdb_file.pdb 8 0 0 1
    # get input arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_code", help="PDB id / protein name", type=str)
    parser.add_argument("fname", help="PDB/PDBx filename", type=str)
    parser.add_argument(
        "d_cut",
        help="Alpha Carbon / atom pairwise contact distance cutoff (in Angstroms)",
        type=float,
    )
    parser.add_argument(
        "o", help="PDB residue index offset integer. Default is 0.", type=int
    )
    parser.add_argument("pdbx", help="set=1 to use pdbx file parser", type=int)
    parser.add_argument(
        "CA_only", help="set=1 to use only alpha carbons (0 for all atoms)", type=int
    )
    args = parser.parse_args()
    main(args)
