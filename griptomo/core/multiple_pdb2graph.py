import argparse
import fnmatch
import os, sys, time
import multiprocessing
import pandas as pd
import griptomo.core.pdb2graph as p2g


# to check cpu count use multiprocessing.cpu_count()

# to check the platform look at the link below
# https://docs.python.org/3/library/sys.html#sys.platform

"""
AIX : 'aix'
Linux : 'linux'
Windows : 'win32'
Windows/Cygwin : 'cygwin'
macOS : 'darwin'
"""

_args = sys.argv[0:]
_py_file = _args[0]


_code_location = os.path.dirname(os.path.abspath(_py_file))

try:
    index_of_latest = _code_location.index("griptomo")
    scripts_path = os.path.join(_code_location[:index_of_latest], "griptomoml", "core")
    sys.path.insert(0, scripts_path)
except ValueError:
    print(f"Error: 'griptomo' not found in the directory path {_code_location}")


_folder_location = _args[1]


def show_time(process, time_start, time_end):
    """
    Calculate and format the time taken for a process.

    Parameters
    ----------
    process : str
        Name of the process.
    time_start : float
        Start time of the process.
    time_end : float
        End time of the process.

    Returns
    -------
    str
        Formatted string showing the time taken for the process.
    """
    time_took = "\n" + str(process) + " finished in "
    if round((time_end - time_start) / 60, 1) < 1:
        time_took = time_took + str(round((time_end - time_start), 1)) + " seconds "
    elif round((time_end - time_start) / 60 / 60, 1) < 1:
        time_took = (
            time_took + str(round((time_end - time_start) / 60, 1)) + " minutes "
        )
    else:
        time_took = (
            time_took + str(round((time_end - time_start) / 60 / 60, 1)) + " hours "
        )
    time_took = time_took + "(wall clock)."
    return time_took


# define the function to call in main
def generate_graph(pdb_code, fname, t, o, pdbx, CA_only):
    """
    Generate a graph from a PDB file.

    Parameters
    ----------
    pdb_code : str
        PDB ID / label for the protein of interest.
    fname : str
        Filename for the protein of interest. Can be PDB or PDBx format.
    t : float
        Alpha Carbon / atom pairwise contact distance cutoff (in Angstroms).
    o : int
        Index offset in case the first residue ID in PDB file is not the first physical residue.
    pdbx : int
        Set to 1 if using the newer PDBx file format.
    CA_only : int
        Set to 1 if using only alpha carbons, else all atoms are used.

    Returns
    -------
    str
        PDB code.
    """
    df = p2g.PDB_to_df(pdb_code, fname, pdbx, o, CA_only)
    G = p2g.PDB_df_to_G(df, t)
    p2g.save_data(df, G, pdb_code, pdb_code)
    return pdb_code


def pdb2graph_list(t, o, pdbx, CA_only=1):
    """
    Convert multiple PDBs to graphs from a folder of PDBs.

    Parameters
    ----------
    t : float
        Alpha Carbon / atom pairwise contact distance cutoff (in Angstroms).
    o : int
        Index offset in case the first residue ID in PDB file is not the first physical residue.
    pdbx : int
        Set to 1 if using the newer PDBx file format.
    CA_only : int, optional
        Set to 1 if using only alpha carbons, else all atoms are used. Default is 1.

    Returns
    -------
    list
        List of results from the graph generation process.
    """
    # List for storing file names and the arguments
    file_list = []
    for file_name in os.listdir(_folder_location):
        if fnmatch.fnmatch(file_name, "*.pdb"):
            pdb_code = file_name[0 : file_name.find(".pdb")] + "_" + str(int(t))
            fname = os.path.join(_folder_location, file_name)
            tuple = (pdb_code, fname, t, o, pdbx, CA_only)
            file_list.append(tuple)

    # process pool for passing pdb2graph_list to multiple processes
    pool = multiprocessing.Pool()
    result = pool.starmap_async(generate_graph, file_list)
    pool.close()
    print(result.get())
    pool.join()
    return result.get()


def main(args):
    """
    Main function to process multiple PDB files and convert them to graphs.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments.
    """
    pdb_code = ""
    fname = ""
    t = args.t
    o = args.o
    pdbx = args.pdbx
    CA_only = args.CA_only

    pdb2graph_list(t, o, pdbx, CA_only)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("f", help="Folder name with the pdb files", type=str)
    parser.add_argument(
        "t", help="Alpha Carbon contact distance threshold (in Angstroms)", type=float
    )
    parser.add_argument(
        "o", help="PDB residue index offset integer. Default is 0.", type=int
    )
    parser.add_argument("pdbx", help="set=1 to use pdbx file parser", type=int)
    parser.add_argument("CA_only", help="set=1 to use only alpha carbons", type=int)
    args = parser.parse_args()

    start_time = time.time()
    main(args)
    print(show_time("pdb to graph", start_time, time.time()))

    # example running: python /griptomo/core/multiple_pdb2graph.py 8 0 0 1

    # parameters used for the 1st paper
    # t = 8 # pairwise distance cutoff for assigning edges, in Angstroms
    # o = 0  # residue indexing offest (default = 0)
    # pdbx = 0  # using .pdb (0) or .pdbx (1) file format
    # CA_only = 1 # using alpha carbons only

    # ref. With 8 cores at mac,         it took 20 minutes (wall clock) to generate graphs for 64 pdb files
    # ref. With 10 M1 max cores at mac, it took 11 minutes (wall clock) to generate graphs for 100 apoferritin-sized pdb files
