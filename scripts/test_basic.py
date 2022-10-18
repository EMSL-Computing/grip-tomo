# August George, 2022, PNNL

import unittest
import os
from pathlib import Path
import density2graph as d2g
import graph2class as g2c
import pdb2graph as p2g
import networkx as nx
import pandas as pd
import numpy as np


class TestGripTomo(unittest.TestCase):
    """ A set of basic unit/integration tests for GRIP-tomo.
    Usage: run python test_basic.py from ./scripts

    Args:
        unittest ([TestCase]): unittest TestCase
    """

    def setUp(self):
       self.data_path =  Path(Path(os.getcwd()).parent,'example_data')  # datafile path
 

    #test if parallelized functions work
    def test_package_versions(self):
        '''checks the package versions of numpy, pandas, and networkx are correct'''
        print(f"networkx version: {nx.__version__}")
        print(f"numpy version: {np.__version__}")
        print(f"pandas version: {pd.__version__}")
        assert(nx.__version__ >= "2.8")
        assert(np.__version__ >= "1.21")
        assert(pd.__version__ >= "1.4")
     

    def test_pdb2graph_CA_only(self):
        ''' creates an example graph from a pdbfile and compares it to the expected output - using alpha carbons only'''

        G_true = nx.read_gexf(Path(self.data_path,'3vjf_pdb2graph.gexf'))

        pdb_code = '3vjf' # name of protein (for output)
        fname = Path(self.data_path,'3vjf.pdb')  # file to turn into graph
        d_cut = 8 # pairwise distance cutoff for assigning edges, in Angstroms
        o = 0  # residue indexing offest (default = 0)
        pdbx = 0  # using .pdb (0) or .pdbx (1) file format
        CA_only = 1  # using only alpha carbons

        # run conversion script
        df = p2g.PDB_to_df(pdb_code, fname, pdbx, o, CA_only)  # convert pdb file into dataframe of atom (only alpha carbon) coordinates and s/n
        df2 = p2g.PDB_to_df(pdb_code, fname, pdbx, o)  # convert pdb file into dataframe of atom (only alpha carbon) coordinates and s/n (check if default works)
        G = p2g.PDB_df_to_G(df, d_cut)  # convert coordinate dataframe into network graph
        G2 = p2g.PDB_df_to_G(df2, d_cut)  # convert coordinate dataframe into network graph

        # test 1 - dataframe is not empty
        self.assertTrue(not df.empty)
        self.assertTrue(not df2.empty)

        # test 2 - graph is not empty
        self.assertTrue(not nx.is_empty(G))
        self.assertTrue(not nx.is_empty(G2))

        # test 3 - graph has the same number of nodes as expected
        self.assertEqual(G.number_of_nodes(), G_true.number_of_nodes())
        self.assertEqual(G2.number_of_nodes(), G_true.number_of_nodes())

        # test 4 - graph has the same number of edges as expected
        self.assertEqual(G.number_of_edges(), G_true.number_of_edges())
        self.assertEqual(G2.number_of_edges(), G_true.number_of_edges())


    def test_density2graph_pdb2density(self):
        ''' creates an example graph from an ideal density (pdb to mrc density) file and compares it to the expected output'''
        G_true = nx.read_gexf(Path(self.data_path,'3vjf_density2graph.gexf'))

        # set parameters
        fname =  Path(self.data_path,'3vjf_density.mrc')
        t = 0.425  # pixel intensity threshold - unnormalized
        DBSCAN_epsilon = 1  # DBSCAN  epsilon
        DBSCAN_min_samples = 4  # DBSCAN min samples
        d_cut = 8  # pairwise cutoff distance in pixels (edge is assigned to two nodes if Euclidean distance <= cutoff distance)
        out_fname = fname.with_suffix('')  # output file name

        # conversion pipeline
        mrc = d2g.load_density_file(fname)  # load density file
        xyz_data = d2g.normalize_and_threshold_data(mrc,t)  # normalize data and threshold and then apply threshold
        model = d2g.cluster_data(xyz_data,DBSCAN_epsilon,DBSCAN_min_samples)  # cluster thresholded data using DBSCAN
        coarse_model = d2g.get_cluster_centroids(xyz_data,model)  # coarse grain model by getting cluster centroids
        G = d2g.create_and_save_graph(coarse_model,d_cut,out_fname, save=False)  # create graph where nodes = centroids, edges assigned by pairwise cutoff
      
        # test 1 - xyz_data is not empty
        self.assertTrue(xyz_data.size > 0)

        # test 2 - coarse_model is not empty
        self.assertTrue(coarse_model.size > 0)

        # test 3 - graph has the same number of nodes as expected
        self.assertEqual(G.number_of_nodes(), G_true.number_of_nodes())

        # test 4 - graph has the same number of edges as expected
        self.assertEqual(G.number_of_edges(), G_true.number_of_edges())


    def test_graph2class(self):
            ''' classifies graphs and compares it to the expected output'''

            class_files = ['3vjf_pdb2graph.gexf']

            sample_files = [
                    '3vjf_density2graph.gexf', 
                    ]

            similarity_features_list = [
                        'n nodes',
                        'n edges',
                        'density',
                        'diameter',
                        'avg path length',
                        'avg clustering',
                        'max closeness centrality',
                        'max eigenvector centrality',
                        'max betweenness centrality',
                        'degree assortativity',
                        'max clique number',
                        'n communities',
                ]
            
            # get graph network features for each control/reference graph (classes) and save to a .csv file 
            class_list = [Path(self.data_path, fname) for fname in class_files]
           
            # get graph network features for each non-control/non-reference graph (sampled) and save to a .csv file
            sample_list = [Path(self.data_path, fname) for fname in sample_files]

            # classify based on network feature similarities and save to a .csv file
            classify_df = g2c.classify_graphs(class_list, sample_list, similarity_features_list)

            # test 1 - check if similarty values match expected value 
            true_class_dict = {
                'name': ['3vjf_density2graph'],
                '3vjf_pdb2graph': [0.895314],
            } 

            true_class_df = pd.DataFrame.from_dict(true_class_dict)
            for col in true_class_df.columns:
                if col != 'name':
                    self.assertTrue(np.allclose(true_class_df[col],classify_df[col]))
           

if __name__ == '__main__':
    unittest.main()
