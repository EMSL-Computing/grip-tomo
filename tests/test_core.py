# August George, 2022, PNNL

import unittest
from pathlib import Path
from griptomo.core import density2graph as d2g
from griptomo.core import graph2class as g2c
from griptomo.core import ig_extract_features as igf
from griptomo.core import pdb2graph as p2g
import networkx as nx
import pandas as pd
import numpy as np
import mrcfile
import igraph as ig
from scipy.spatial.distance import pdist, squareform

if d2g.IS_GPU_AVAILABLE:
    import cupy as cp
    from griptomo.core import gpu_extract_features as gef
    import cugraph as cg
    import cudf


class TestGripTomo(unittest.TestCase):
    """
    A set of basic unit tests for GRIP-tomo.

    Usage:
        Run via pytest within the repository root.

    Attributes
    ----------
    data_path : Path
        Path to the example data directory.
    """

    def setUp(self):
        """Load shared example data paths and detect GPU availability once per test."""
        repo_root = Path(__file__).resolve().parents[1]
        self.data_path = repo_root / "tests" / "fixtures" / "example_data"
        self.is_gpu_available = d2g.IS_GPU_AVAILABLE

    # test if parallelized functions work
    def test_package_versions(self):
        """
        Checks the package versions of numpy, pandas, and networkx are correct.
        """
        print(f"networkx version: {nx.__version__}")
        print(f"numpy version: {np.__version__}")
        print(f"pandas version: {pd.__version__}")
        assert nx.__version__ >= "2.8"
        assert np.__version__ >= "1.21"
        assert pd.__version__ >= "1.4"

    def test_pdb2graph_CA_only(self):
        """
        Creates an example graph from a pdbfile and compares it to the expected output using alpha carbons only.
        """
        G_true = nx.read_gexf(Path(self.data_path, "3vjf.gexf"))

        pdb_code = "3vjf"  # name of protein (for output)
        fname = Path(self.data_path, "3vjf.pdb")  # file to turn into graph
        d_cut = 8  # pairwise distance cutoff for assigning edges, in Angstroms
        o = 0  # residue indexing offest (default = 0)
        pdbx = 0  # using .pdb (0) or .pdbx (1) file format
        CA_only = 1  # using only alpha carbons

        # run conversion script
        df = p2g.PDB_to_df(
            pdb_code, fname, pdbx, o, CA_only
        )  # convert pdb file into dataframe of atom (only alpha carbon) coordinates and s/n
        df2 = p2g.PDB_to_df(
            pdb_code, fname, pdbx, o
        )  # convert pdb file into dataframe of atom (only alpha carbon) coordinates and s/n (check if default works)
        G = p2g.PDB_df_to_G(
            df, d_cut
        )  # convert coordinate dataframe into network graph
        G2 = p2g.PDB_df_to_G(
            df2, d_cut
        )  # convert coordinate dataframe into network graph

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

    def test_pdb2graph_all_atom(self):
        """
        Creates an example graph from a pdbfile and compares it to the expected output using all atoms.
        """
        # Note: the 3vjf_all_atom.gexf file has no water molecules.
        G_true = nx.read_gexf(Path(self.data_path, "3vjf_all_atom.gexf"))

        pdb_code = "3vjf_all_atom"  # name of protein (for output)
        fname = Path(self.data_path, "3vjf.pdb")  # file to turn into graph
        d_cut = 8  # pairwise distance cutoff for assigning edges, in Angstroms
        o = 0  # residue indexing offest (default = 0)
        pdbx = 0  # using .pdb (0) or .pdbx (1) file format
        CA_only = 0  # using all atoms

        # run conversion script
        df = p2g.PDB_to_df(
            pdb_code, fname, pdbx, o, CA_only
        )  # convert pdb file into dataframe of atom coordinates and s/n
        G = p2g.PDB_df_to_G(
            df, d_cut
        )  # convert coordinate dataframe into network graph

        # test 1 - dataframe is not empty
        self.assertTrue(not df.empty)

        # test 2 - graph is not empty
        self.assertTrue(not nx.is_empty(G))

        # test 3 - graph has the same number of nodes as expected
        self.assertEqual(G.number_of_nodes(), G_true.number_of_nodes())

        # test 4 - graph has the same number of edges as expected
        self.assertEqual(G.number_of_edges(), G_true.number_of_edges())

    def test_density2graph_pdb2density(self):
        """
        Creates an example graph from an ideal density (pdb to mrc density) file and compares it to the expected output.
        """
        G_true = nx.read_gexf(
            Path(self.data_path, "3vjf_CA_apix_1_box_dimension_400_res_3.gexf")
        )

        # set parameters
        fname = Path(self.data_path, "3vjf_CA_apix_1_box_dimension_150_res_3.mrc")
        t = 0.425  # pixel intensity threshold - unnormalized
        DBSCAN_epsilon = 1  # DBSCAN  epsilon
        DBSCAN_min_samples = 4  # DBSCAN min samples
        d_cut = 8  # pairwise cutoff distance in pixels (edge is assigned to two nodes if Euclidean distance <= cutoff distance)
        out_fname = fname.with_suffix("")  # output file name

        # conversion pipeline
        mrc = d2g.load_density_file(fname)  # load density file
        xyz_data = d2g.normalize_and_threshold_data(
            mrc, t
        )  # normalize data and threshold and then apply threshold
        model = d2g.cluster_data(
            xyz_data, DBSCAN_epsilon, DBSCAN_min_samples
        )  # cluster thresholded data using DBSCAN
        coarse_model = d2g.get_cluster_centroids(
            xyz_data, model
        )  # coarse grain model by getting cluster centroids
        G = d2g.create_and_save_graph(
            coarse_model, d_cut, out_fname, save=False
        )  # create graph where nodes = centroids, edges assigned by pairwise cutoff

        # test 1 - xyz_data is not empty
        self.assertTrue(xyz_data.size > 0)

        # test 2 - coarse_model is not empty
        self.assertTrue(coarse_model.size > 0)

        # test 3 - graph has the same number of nodes as expected
        self.assertEqual(G.number_of_nodes(), G_true.number_of_nodes())

        # test 4 - graph has the same number of edges as expected
        self.assertEqual(G.number_of_edges(), G_true.number_of_edges())

    def test_density2graph_tomo_density(self):
        """
        Creates an example graph from a reconstructed density (pdb to synthetic tomogram) file and compares it to the expected output.
        """
        G_true = nx.read_gexf(
            Path(
                self.data_path,
                "3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed.gexf",
            )
        )

        # set parameters
        fname = Path(
            self.data_path, "3vjf_CA_apix_1_box_dimension_150_res_3_reconstructed.mrc"
        )
        t = 26985  # pixel intensity threshold - unnormalized (found empirically)
        DBSCAN_epsilon = 1  # DBSCAN  epsilon
        DBSCAN_min_samples = 4  # DBSCAN min samples
        d_cut = 8  # pairwise cutoff distance in pixels (edge is assigned to two nodes if Euclidean distance <= cutoff distance)
        out_fname = fname.with_suffix("")  # output file name

        # conversion pipeline
        mrc = d2g.load_density_file(fname)  # load density file
        xyz_data = d2g.normalize_and_threshold_data(
            mrc, t
        )  # normalize data and threshold and then apply threshold
        model = d2g.cluster_data(
            xyz_data, DBSCAN_epsilon, DBSCAN_min_samples
        )  # cluster thresholded data using DBSCAN
        coarse_model = d2g.get_cluster_centroids(
            xyz_data, model
        )  # coarse grain model by getting cluster centroids
        G = d2g.create_and_save_graph(
            coarse_model, d_cut, out_fname, save=False
        )  # create graph where nodes = centroids, edges assigned by pairwise cutoff

        # test 1 - xyz_data is not empty
        self.assertTrue(xyz_data.size > 0)

        # test 2 - coarse_model is not empty
        self.assertTrue(coarse_model.size > 0)

        # test 3 - graph has the same number of nodes as expected
        self.assertEqual(G.number_of_nodes(), G_true.number_of_nodes())

        # test 4 - graph has the same number of edges as expected
        self.assertEqual(G.number_of_edges(), G_true.number_of_edges())

    def test_graph2class(self):
        """
        Classifies graphs and compares it to the expected output.
        """
        class_files = ["3vjf.gexf", "4rlc.gexf", "4rly.gexf"]

        sample_files = [
            "3vjf_CA_apix_1_box_dimension_400_res_3.gexf",
            "3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed.gexf",
        ]

        similarity_features_list = [
            "n nodes",
            "n edges",
            "density",
            "diameter",
            "avg path length",
            "avg clustering",
            "max closeness centrality",
            "max eigenvector centrality",
            "max betweenness centrality",
            "degree assortativity",
            "max clique number",
            "n communities",
        ]

        # get graph network features for each control/reference graph (classes) and save to a .csv file
        class_list = [Path(self.data_path, fname) for fname in class_files]

        # get graph network features for each non-control/non-reference graph (sampled) and save to a .csv file
        sample_list = [Path(self.data_path, fname) for fname in sample_files]

        # classify based on network feature similarities and save to a .csv file
        classify_df = g2c.classify_graphs(
            class_list, sample_list, similarity_features_list, 1
        )
        # print(classify_df)

        # test 1 - check if similarity values match expected value using networkx
        # true_class_dict = {
        #     'name': ['3vjf_CA_apix_1_box_dimension_400_res_3', '3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed'],
        #     '3vjf': [0.895314, 0.815112],
        #     '4rlc': [0.744599, 0.693560],
        #     '4rly': [0.749838, 0.759750]
        # }

        true_class_dict = {
            "name": [
                "3vjf_CA_apix_1_box_dimension_400_res_3",
                "3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed",
            ],
            "3vjf": [0.914831, 0.834077],
            "4rlc": [0.761140, 0.709476],
            "4rly": [0.760443, 0.769583],
        }

        true_class_df = pd.DataFrame.from_dict(true_class_dict)
        for col in true_class_df.columns:
            if col != "name":
                # print("TEST 1")
                # print(f"{col} true values:\n{true_class_df[col]}")  # debugging
                # print(f"{col} predicted values:\n{classify_df[col]}\n")  # debugging
                self.assertTrue(np.allclose(true_class_df[col], classify_df[col]))

        # test 2 - check if similarity values match expected value using igraph
        classify_df_igraph = g2c.classify_graphs(
            class_list, sample_list, similarity_features_list, 2
        )
        true_class_dict_igraph = {
            "name": [
                "3vjf_CA_apix_1_box_dimension_400_res_3",
                "3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed",
            ],
            "3vjf": [0.931504, 0.850751],
            "4rlc": [0.777810, 0.726148],
            "4rly": [0.774519, 0.783632],
        }
        true_class_df_igraph = pd.DataFrame.from_dict(true_class_dict_igraph)
        for col in true_class_df_igraph.columns:
            if col != "name":
                print("TEST 2")
                print(f"{col} true values:\n{true_class_df_igraph[col]}")  # debugging
                print(
                    f"{col} predicted values:\n{classify_df_igraph[col]}\n"
                )  # debuggingdifference = np.abs(true_class_df_igraph[col] - classify_df_igraph[col])
                self.assertTrue(
                    np.allclose(true_class_df_igraph[col], classify_df_igraph[col])
                )

    def test_parallel_graph_feature_calculation(self):
        """
        Checks that parallel calculations for avg path length and betweenness centrality are equal to the serial calculation.
        """
        file_list = ["3vjf_all_atom.gexf", "3vjf.gexf", "4rlc.gexf", "4rly.gexf"]
        graph_list = [Path(self.data_path, fname) for fname in file_list]
        true_values = {
            # old values with h2o'3vjf_all_atom' : [4.2964395564404665, 0.011613536965548748],
            "3vjf_all_atom": [4.2888853996097, 0.012322327209593288],
            "3vjf": [5.263030507711359, 0.09337110821076301],
            "4rlc": [3.812382531785517, 0.06094428517917971],
            "4rly": [5.983015355979525, 0.20395424151212188],
        }
        for i, j in zip(graph_list, true_values):
            calculated_features_nx = g2c.calc_graph_features(nx.read_gexf(i))
            self.assertEqual(
                calculated_features_nx["avg path length"], true_values[j][0]
            )
            self.assertEqual(
                calculated_features_nx["max betweenness centrality"], true_values[j][1]
            )
            calculated_features_ig = g2c.igraph_calc_graph_features(nx.read_gexf(i))
            self.assertEqual(
                calculated_features_ig["avg path length"], round(true_values[j][0], 16)
            )
            self.assertEqual(
                calculated_features_ig["max betweenness centrality"],
                round(true_values[j][1], 15),
            )
        # Compare the values from networkx to values from igraph (except for # of communities)
        for i in graph_list:
            calculated_features_nx = g2c.calc_graph_features(nx.read_gexf(i))
            calculated_features_ig = g2c.igraph_calc_graph_features(nx.read_gexf(i))
            for j in calculated_features_nx:
                if j == "n communities":
                    continue
                self.assertEqual(
                    round(calculated_features_nx[j], 2),
                    round(calculated_features_ig[j], 2),
                )

    @unittest.skipUnless(d2g.IS_GPU_AVAILABLE, "GPU is not available")
    def test_cugraph_calc_graph_features(self):
        """
        Tests the cugraph_calc_graph_features function.
        """
        G_ig = ig.Graph.Famous("Zachary")  # Example igraph graph
        G_ig.to_undirected()
        edgelist = G_ig.get_edgelist()
        df = pd.DataFrame(edgelist, columns=["source", "destination"])
        gdf = cudf.DataFrame.from_pandas(df)
        symmetrized_gdf = cg.symmetrize_df(gdf, "source", "destination")
        G_cg = cg.Graph(directed=False)
        G_cg.from_cudf_edgelist(
            symmetrized_gdf, source="source", destination="destination"
        )

        self.assertEqual(G_cg.number_of_nodes(), G_ig.vcount())
        self.assertEqual(G_cg.number_of_edges(), G_ig.ecount())

        features = gef.cugraph_calc_graph_features(G_ig, G_cg, skip_clique_num=True)
        self.assertIsInstance(features, dict)
        self.assertIn("total n nodes", features)
        self.assertIn("total n edges", features)
        self.assertIn("n nodes", features)
        self.assertIn("n edges", features)
        self.assertIn("density", features)
        self.assertIn("diameter", features)
        self.assertIn("avg clustering", features)
        self.assertIn("max closeness centrality", features)
        self.assertIn("max eigenvector centrality", features)
        self.assertIn("degree assortativity", features)
        self.assertIn("n communities", features)
        self.assertIn("avg path length", features)
        self.assertIn("max betweenness centrality", features)

    @unittest.skipUnless(d2g.IS_GPU_AVAILABLE, "GPU is not available")
    def test_cugraph_calc_graph_features_timed(self):
        """
        Tests the cugraph_calc_graph_features_timed function.
        """
        G_ig = ig.Graph.Famous("Zachary")  # Example igraph graph
        G_ig.to_undirected()
        edgelist = G_ig.get_edgelist()
        df = pd.DataFrame(edgelist, columns=["source", "destination"])
        gdf = cudf.DataFrame.from_pandas(df)
        symmetrized_gdf = cg.symmetrize_df(gdf, "source", "destination")
        G_cg = cg.Graph(directed=False)
        G_cg.from_cudf_edgelist(
            symmetrized_gdf, source="source", destination="destination"
        )

        features, times = gef.cugraph_calc_graph_features_timed(
            G_ig, G_cg, skip_clique_num=True
        )
        self.assertIsInstance(features, dict)
        self.assertIsInstance(times, dict)
        self.assertIn("total n nodes", features)
        self.assertIn("total n edges", features)
        self.assertIn("n nodes", features)
        self.assertIn("n edges", features)
        self.assertIn("density", features)
        self.assertIn("diameter", features)
        self.assertIn("avg clustering", features)
        self.assertIn("max closeness centrality", features)
        self.assertIn("max eigenvector centrality", features)
        self.assertIn("degree assortativity", features)
        self.assertIn("n communities", features)
        self.assertIn("avg path length", features)
        self.assertIn("max betweenness centrality", features)
        self.assertIn("total n nodes", times)
        self.assertIn("total n edges", times)
        self.assertIn("n nodes", times)
        self.assertIn("n edges", times)
        self.assertIn("density", times)
        self.assertIn("diameter", times)
        self.assertIn("avg clustering", times)
        self.assertIn("max closeness centrality", times)
        self.assertIn("max eigenvector centrality", times)
        self.assertIn("degree assortativity", times)
        self.assertIn("n communities", times)
        self.assertIn("avg path length", times)
        self.assertIn("max betweenness centrality", times)

    @unittest.skipUnless(d2g.IS_GPU_AVAILABLE, "GPU is not available")
    def test_compare_cpu_gpu_features(self):
        """
        Compares the graph features calculated on CPU and GPU.
        """
        G_ig = ig.Graph.Famous("Zachary")  # Example igraph graph
        G_ig.to_undirected()
        edgelist = G_ig.get_edgelist()
        df = pd.DataFrame(edgelist, columns=["source", "destination"])
        gdf = cudf.DataFrame.from_pandas(df)
        symmetrized_gdf = cg.symmetrize_df(gdf, "source", "destination")
        G_cg = cg.Graph(directed=False)
        G_cg.from_cudf_edgelist(
            symmetrized_gdf, source="source", destination="destination"
        )

        # self.assertTrue(G_cg.has_self_loops() == False)

        # Calculate features on CPU
        cpu_features, cpu_times = igf.igraph_calc_graph_features_timed(
            G_ig, skip_clique_num=True
        )

        # Calculate features on GPU
        gpu_features, gpu_times = gef.cugraph_calc_graph_features_timed(
            G_ig, G_cg, skip_clique_num=True
        )

        # Compare features
        for key in cpu_features:
            print(f"Testing feature: {key}")
            if key != "n communities":
                self.assertAlmostEqual(
                    cpu_features[key],
                    gpu_features[key],
                    places=2,
                    msg=f"Feature {key} does not match between CPU and GPU.",
                )

    def test_normalize_mrc_data(self):
        """
        Tests if a test mrc file is normalized correctly.
        """
        # generate test MRC data
        fname = "test.mrc"
        mrc = mrcfile.open(Path(self.data_path, fname), mode="r+")
        numpy_data = np.float16(
            np.random.uniform(low=-100, high=100, size=(100, 100, 100))
        )
        mrc.set_data(numpy_data)
        header_string = np.frombuffer(b"Test MRC data (random)".ljust(52), dtype="S52")
        mrc.set_extended_header(header_string)

        # test when input is an mrc file, using default lower and upper bounds [-1,1]
        normalized_mrc = d2g.normalize_mrc_data(mrc)
        assert normalized_mrc.min() == -1
        assert normalized_mrc.max() == 1

        # test when input is an mrc file, using [0,1] bounds
        normalized_mrc_0_to_1 = d2g.normalize_mrc_data(mrc, 0, 1)
        assert normalized_mrc_0_to_1.min() == 0
        assert normalized_mrc_0_to_1.max() == 1

        # test when input is an mrc file, using [1,10] bounds
        normalized_mrc_1_to_10 = d2g.normalize_mrc_data(mrc, 1, 10)
        assert normalized_mrc_1_to_10.min() == 1
        assert normalized_mrc_1_to_10.max() == 10
        mrc.close()

        # test when input is numpy array, using default lower and upper bounds [-1,1]
        normalized_data = d2g.normalize_mrc_data(numpy_data)
        assert normalized_data.min() == -1
        assert normalized_data.max() == 1

        # test when input is an mrc file, using [0,1] bounds
        normalized_data_0_to_1 = d2g.normalize_mrc_data(numpy_data, 0, 1)
        assert normalized_data_0_to_1.min() == 0
        assert normalized_data_0_to_1.max() == 1

        # test when input is an mrc file, using [1,10] bounds
        normalized_data_1_to_10 = d2g.normalize_mrc_data(numpy_data, 1, 10)
        assert normalized_data_1_to_10.min() == 1
        assert normalized_data_1_to_10.max() == 10

        if self.is_gpu_available:
            # test when input is a cupy array, using default lower and upper bounds [-1,1]
            normalized_data_gpu = d2g.normalize_mrc_data(
                cp.asarray(numpy_data), device="gpu"
            )
            assert cp.allclose(normalized_data_gpu.min(), -1)
            assert cp.allclose(normalized_data_gpu.max(), 1)

            # test when input is a cupy array, using [0,1] bounds
            normalized_data_0_to_1_gpu = d2g.normalize_mrc_data(
                cp.asarray(numpy_data), 0, 1
            )
            assert cp.allclose(normalized_data_0_to_1_gpu.min(), 0)
            assert cp.allclose(normalized_data_0_to_1_gpu.max(), 1)

            # test when input is a cupy array, using [1,10] bounds
            normalized_data_1_to_10_gpu = d2g.normalize_mrc_data(
                cp.asarray(numpy_data), 1, 10
            )
            assert cp.allclose(normalized_data_1_to_10_gpu.min(), 1)
            assert cp.allclose(normalized_data_1_to_10_gpu.max(), 10)

    def test_threshold_mrc_data(self):
        """
        Tests if data from MRC is correctly thresholded.
        """
        # generate test MRC data
        threshold = 50  # arbitrary threshold value for testing
        fname = "test.mrc"
        mrc = mrcfile.open(Path(self.data_path, fname), mode="r+")
        numpy_data = np.float16(
            np.random.uniform(low=-100, high=100, size=(100, 100, 100))
        )
        mrc.set_data(numpy_data)
        header_string = np.frombuffer(b"Test MRC data (random)".ljust(52), dtype="S52")
        mrc.set_extended_header(header_string)

        # test when input is an mrc file
        thresholded_mrc_data = d2g.threshold_mrc_data(mrc, threshold)
        assert np.all(
            (thresholded_mrc_data == np.min(numpy_data))
            | (thresholded_mrc_data > threshold)
        )
        mrc.close()

        # test when input is numpy array
        thresholded_data = d2g.threshold_mrc_data(numpy_data, threshold)
        assert np.all(
            (thresholded_data == np.min(numpy_data)) | (thresholded_data > threshold)
        )

        if self.is_gpu_available:
            # test when input is a cupy array
            thresholded_data_gpu = d2g.threshold_mrc_data(
                cp.asarray(numpy_data), threshold, device="gpu"
            )
            assert cp.all(
                (thresholded_data_gpu == cp.min(cp.asarray(numpy_data)))
                | (thresholded_data_gpu > threshold)
            )

    def test_generate_point_cloud_from_mrc_data(self):
        """
        Tests if data from MRC is correctly converted into point cloud (xyz coordinates).
        """
        # generate test MRC data
        threshold = 50  # arbitrary threshold value for testing
        fname = "test.mrc"
        mrc = mrcfile.open(Path(self.data_path, fname), mode="r+")
        numpy_data = np.float16(
            np.random.uniform(low=-100, high=100, size=(100, 100, 100))
        )
        mrc.set_data(numpy_data)
        header_string = np.frombuffer(b"Test MRC data (random)".ljust(52), dtype="S52")
        mrc.set_extended_header(header_string)

        # test when input is an mrc file
        threshold_mrc_xyz = d2g.generate_point_cloud_from_mrc_data(mrc, threshold)
        assert np.all(
            numpy_data[
                threshold_mrc_xyz[:, 0],
                threshold_mrc_xyz[:, 1],
                threshold_mrc_xyz[:, 2],
            ]
            > threshold
        )
        mrc.close()

        # test when input is numpy array
        threshold_data_xyz = d2g.generate_point_cloud_from_mrc_data(
            numpy_data, threshold
        )
        assert np.all(
            numpy_data[
                threshold_data_xyz[:, 0],
                threshold_data_xyz[:, 1],
                threshold_data_xyz[:, 2],
            ]
            > threshold
        )

        if self.is_gpu_available:
            # test when input is a cupy array
            threshold_data_xyz_gpu = d2g.generate_point_cloud_from_mrc_data(
                cp.asarray(numpy_data), threshold, device="gpu"
            )
            assert cp.all(
                cp.asarray(numpy_data)[
                    threshold_data_xyz_gpu[:, 0],
                    threshold_data_xyz_gpu[:, 1],
                    threshold_data_xyz_gpu[:, 2],
                ]
                > threshold
            )

    def test_augment_mrc_data(self):
        """
        Tests if data from MRC is augmented correctly.
        """
        # generate test MRC data
        offset_percent = 10  # arbitrary percent offset
        fname = "test.mrc"
        mrc = mrcfile.open(Path(self.data_path, fname), mode="r+")
        mrc_data = np.float16(
            np.random.uniform(low=-100, high=100, size=(100, 100, 100))
        )
        mrc.set_data(mrc_data)
        header_string = np.frombuffer(b"Test MRC data (random)".ljust(52), dtype="S52")
        mrc.set_extended_header(header_string)

        # test when input is an mrc file and n=1
        n = 1
        adjusted_mrc_data = d2g.augment_mrc_data(mrc, offset_percent, n)
        relative_difference = np.abs(adjusted_mrc_data - mrc_data) / mrc_data
        assert not np.all(np.equal(adjusted_mrc_data, mrc_data))
        assert np.all(relative_difference <= offset_percent / 100)

        # test when input is an mrc file and n=3
        n = 3
        adjusted_mrc_data_list = d2g.augment_mrc_data(mrc, offset_percent, n)
        for adjusted_mrc_data in adjusted_mrc_data_list:
            relative_difference = np.abs(adjusted_mrc_data - mrc_data) / mrc_data
            assert not np.all(np.equal(adjusted_mrc_data, mrc_data))
            assert np.all(relative_difference <= offset_percent / 100)
        mrc.close()

        # test when input is a numpy array and n=1
        n = 1
        adjusted_mrc_data = d2g.augment_mrc_data(mrc_data, offset_percent, n)
        relative_difference = np.abs(adjusted_mrc_data - mrc_data) / mrc_data
        assert not np.all(np.equal(adjusted_mrc_data, mrc_data))
        assert np.all(relative_difference <= offset_percent / 100)

        # test when input is a numpy array and n=3
        n = 3
        adjusted_mrc_data_list = d2g.augment_mrc_data(mrc_data, offset_percent, n)
        for adjusted_mrc_data in adjusted_mrc_data_list:
            relative_difference = np.abs(adjusted_mrc_data - mrc_data) / mrc_data
            assert not np.all(np.equal(adjusted_mrc_data, mrc_data))
            assert np.all(relative_difference <= offset_percent / 100)

    def test_calculate_point_cloud_mass(self):
        """
        Tests whether protein mass is calculated as expected.
        """
        # Create test data set
        threshold = 50
        expected_number_of_voxels_above_threshold = 1000
        mrc_data = np.full((100, 100, 100), threshold - 1)
        indices_above_threshold = np.unravel_index(
            np.random.choice(
                mrc_data.size,
                size=expected_number_of_voxels_above_threshold,
                replace=False,
            ),
            (100, 100, 100),
        )
        mrc_data[indices_above_threshold] = np.random.uniform(
            low=threshold + 1,
            high=threshold + 100,
            size=expected_number_of_voxels_above_threshold,
        )

        point_cloud_data = d2g.generate_point_cloud_from_mrc_data(mrc_data, threshold)
        assert point_cloud_data.shape[0] == expected_number_of_voxels_above_threshold

        voxel_dim_nm3 = 1.2
        density_kDa_nm3 = 0.82
        expected_mass = (
            expected_number_of_voxels_above_threshold * voxel_dim_nm3 * density_kDa_nm3
        )
        actual_mass = d2g.calculate_point_cloud_mass(
            point_cloud_data, voxel_dim_nm3, density_kDa_nm3
        )
        assert actual_mass == expected_mass

    def test_identify_threshold_ratio(self):
        """
        Tests whether threshold ratio is calculated as expected.

        Note: This test validates the identify_threshold_ratio function which
        performs a binary search to find the threshold value that produces a
        target ratio of voxels above threshold. The test is marked for review
        to ensure edge cases are properly covered (e.g., target ratio = 0.0 or 1.0,
        uniform density data, etc.).

        TODO: Add test cases for:
        - Edge case: target_ratio = 0.0 (should return max density)
        - Edge case: target_ratio = 1.0 (should return min density)
        - Uniform density array (all values identical)
        - Verify binary search convergence tolerance
        """
        ratios = [0.0, 0.1, 0.5, 0.7, 1]

        fname = Path(self.data_path, "3vjf_CA_apix_1_box_dimension_150_res_3.mrc")
        mrc = d2g.load_density_file(fname)  # load density file
        normalized_mrc = d2g.normalize_mrc_data(mrc)

        mrc_datas = [normalized_mrc]
        for mrc_data in mrc_datas:
            for ratio in ratios:
                print(f"Testing ratio {ratio} on mrc_data of shape {mrc_data.shape}")
                count_before = np.sum(mrc_data > -1)
                try:
                    threshold_value = d2g.identify_threshold_ratio(mrc_data, ratio)
                except ValueError:
                    if ratio == 1.0:
                        print("Good job! You caught the ValueError!")
                        continue
                    else:
                        raise ValueError("You should not have caught a ValueError!")
                thresholded_data = d2g.threshold_mrc_data(mrc_data, threshold_value)
                count_after = np.sum(thresholded_data > -1)
                if ratio == 0.0:
                    assert count_before == count_after
                else:
                    expected_count = min(np.size(mrc_data) * (1 - ratio), count_before)
                    actual_count = count_after
                    if expected_count != actual_count:
                        print(
                            f"Expected count: {expected_count}, Actual count: {actual_count}, ratio {ratio}"
                        )

    def test_standardize_mrc_data(self):
        """
        Tests if a test mrc file is standardized correctly.
        """
        # generate test MRC data
        fname = "test.mrc"
        mrc = mrcfile.open(Path(self.data_path, fname), mode="r+")
        numpy_data = np.float32(
            np.random.uniform(low=-100, high=100, size=(100, 100, 100))
        )
        mrc.set_data(numpy_data)
        header_string = np.frombuffer(b"Test MRC data (random)".ljust(52), dtype="S52")
        mrc.set_extended_header(header_string)

        # test when input is an mrc file
        standardized_mrc = d2g.standardize_mrc_data(mrc)
        if not np.isclose(standardized_mrc.mean(), 0, atol=1e-2):
            raise ValueError(f"Mean is not close to 0: {standardized_mrc.mean()}")
        if not np.isclose(standardized_mrc.std(), 1, atol=1e-2):
            raise ValueError(
                f"Standard deviation is not close to 1: {standardized_mrc.std()}"
            )

        # test when input is numpy array
        standardized_data = d2g.standardize_mrc_data(numpy_data)
        if not np.isclose(standardized_data.mean(), 0, atol=1e-2):
            raise ValueError(f"Mean is not close to 0: {standardized_data.mean()}")
        if not np.isclose(standardized_data.std(), 1, atol=1e-2):
            raise ValueError(
                f"Standard deviation is not close to 1: {standardized_data.std()}"
            )

        if self.is_gpu_available:
            # test when input is a cupy array
            standardized_data_gpu = d2g.standardize_mrc_data(
                cp.asarray(numpy_data), device="gpu"
            )
            if not cp.isclose(standardized_data_gpu.mean(), 0, atol=1e-2):
                raise ValueError(
                    f"Mean is not close to 0: {standardized_data_gpu.mean()}"
                )
            if not cp.isclose(standardized_data_gpu.std(), 1, atol=1e-2):
                raise ValueError(
                    f"Standard deviation is not close to 1: {standardized_data_gpu.std()}"
                )

        mrc.close()

    @unittest.skipUnless(d2g.IS_GPU_AVAILABLE, "GPU is not available")
    def test_create_cugraph_from_point_cloud(self):
        """
        Tests if a point cloud is correctly converted into a cuGraph object.
        """
        # generate test point cloud data
        num_nodes = 100
        point_cloud_data = cp.random.uniform(low=0, high=10, size=(num_nodes, 3))

        cutoff = 5

        # test when input is a numpy array
        G = d2g.create_cugraph_from_point_cloud(point_cloud_data, 5)

        dist_matrix = squareform(pdist(cp.asnumpy(point_cloud_data), "euclid"))
        # if the distance is > t, replace it with 0 (i.e. remove edge)
        dist_matrix_thresh = np.where(dist_matrix > cutoff, 0, dist_matrix)
        # create an unweighted, undirected graph from the distance matrix
        G_reference = ig.Graph.Adjacency(
            (dist_matrix_thresh > 0).tolist(), mode="undirected"
        )

        self.assertEqual(G.number_of_nodes(), G_reference.vcount())
        self.assertEqual(G.number_of_edges(), G_reference.ecount())

    def test_create_igraph_from_point_cloud(self):
        """
        Tests if a point cloud is correctly converted into an igraph object.
        """
        # generate test point cloud data
        point_cloud_data = np.float32(np.random.uniform(low=0, high=10, size=(100, 3)))

        cutoff = 5

        # test when input is a numpy array
        G = d2g.create_igraph_from_point_cloud(point_cloud_data, 5)

        dist_matrix = squareform(pdist(point_cloud_data, "euclid"))
        # if the distance is > t, replace it with 0 (i.e. remove edge)
        dist_matrix_thresh = np.where(dist_matrix > cutoff, 0, dist_matrix)
        # create an unweighted, undirected graph from the distance matrix
        G_reference = ig.Graph.Adjacency(
            (dist_matrix_thresh > 0).tolist(), mode="undirected"
        )

        self.assertTrue(G.isomorphic(G_reference))


if __name__ == "__main__":
    unittest.main()
