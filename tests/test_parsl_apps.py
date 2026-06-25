import unittest
import os
from pathlib import Path
import tempfile
import numpy as np
import networkx as nx
import tracemalloc

try:
    # Parsl is optional; grab it if present so the smoke test can run.
    import parsl  # type: ignore[import-not-found]
except ImportError:  # pragma: no cover - optional dependency
    parsl = None  # type: ignore[assignment]
    PARSL_AVAILABLE = False
    make_ci_config = None  # type: ignore[assignment]
    python_app = None  # type: ignore[assignment]
else:
    # Parsl is installed locally; reuse the lightweight CI config and python_app decorator.
    from griptomo.parsl.parsl_config_ci import make_ci_config
    from parsl import python_app  # type: ignore[import-not-found]

    PARSL_AVAILABLE = True

RUN_PARS_TESTS = (
    os.getenv("RUN_PARS_TESTS") == "1"
)  # heavier workflows opt-in via env flag

if PARSL_AVAILABLE and RUN_PARS_TESTS:
    from griptomo.parsl.apps import (
        convert_density_file_to_point_cloud_hdbscan_timed,
        extract_graph_features_from_centroids_timed,
    )
else:  # pragma: no cover - optional dependency
    convert_density_file_to_point_cloud_hdbscan_timed = None  # type: ignore[assignment]
    extract_graph_features_from_centroids_timed = None  # type: ignore[assignment]


@unittest.skipUnless(PARSL_AVAILABLE, "Parsl not installed")
class TestParslSmoke(unittest.TestCase):
    """Lightweight Parsl smoke test that runs by default when Parsl is present."""

    def setUp(self):
        """Start a clean Parsl DFK using the CI configuration."""
        parsl.clear()
        parsl.load(make_ci_config())

    def tearDown(self):
        """Ensure Parsl state is cleared after each smoke test."""
        parsl.clear()

    def test_parallel_betweenness_matches_networkx(self):
        """Verify calc_bc executed via Parsl matches the NetworkX baseline."""
        repo_root = Path(__file__).resolve().parents[1]
        graph_path = repo_root / "tests" / "fixtures" / "example_data" / "3vjf.gexf"
        graph = nx.read_gexf(graph_path)
        expected = max(nx.betweenness_centrality(graph).values())

        @python_app
        def betweenness_task(gexf_path):
            import networkx as _nx
            from multiprocessing import Manager
            from griptomo.core.graph2class import calc_bc

            G = _nx.read_gexf(gexf_path)
            with Manager() as manager:
                result_dict = manager.dict()
                calc_bc(G, result_dict)
                return result_dict[1]

        future1 = betweenness_task(str(graph_path))
        future2 = betweenness_task(str(graph_path))

        value1 = future1.result()
        value2 = future2.result()

        self.assertAlmostEqual(value1, expected, places=12)
        self.assertAlmostEqual(value2, expected, places=12)
        self.assertAlmostEqual(value1, value2, places=12)


@unittest.skipUnless(
    PARSL_AVAILABLE and RUN_PARS_TESTS,
    "Parsl tests require RUN_PARS_TESTS=1 and Parsl installed",
)
class TestParslApps(unittest.TestCase):
    """
    A set of unit tests for griptomo.parsl.apps workflows.

    Usage:
        Run via pytest within the repository root.

    Attributes
    ----------
    data_path : Path
        Path to the example data directory.
    output_dir : Path
        Path to the output directory.
    """

    def setUp(self):
        """Prime temporary directories and Parsl runtime for workflow tests."""
        repo_root = Path(__file__).resolve().parents[1]
        self.data_path = repo_root / "tests" / "fixtures" / "example_data"
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)  # output path

        parsl.load(make_ci_config())
        tracemalloc.start()

    def tearDown(self):
        """Dispose of Parsl resources and temporary directories after tests."""
        tracemalloc.stop()
        parsl.clear()
        self.temp_dir.cleanup()

    def log_memory_usage(self, test_name):
        """Print current and peak memory usage for long-running workflows."""
        current, peak = tracemalloc.get_traced_memory()
        print(
            f"{test_name} - Current memory usage: {current / 1024**2:.2f} MB; Peak: {peak / 1024**2:.2f} MB"
        )

    def test_convert_density_file_to_point_cloud_hdbscan_timed(self):
        """
        Tests the convert_density_file_to_point_cloud_hdbscan_timed workflow.
        """
        self.log_memory_usage(
            "test_convert_density_file_to_point_cloud_hdbscan_timed begin"
        )

        protein_name = "test_protein"
        mrc_file_path = Path(self.data_path, "test.mrc")
        HDBSCAN_min_cluster_size = 6
        HDBSCAN_min_samples = 6
        coarsening_dim = 3
        density_threshold = 0.98
        output_dir = self.output_dir

        result = convert_density_file_to_point_cloud_hdbscan_timed(
            protein_name,
            mrc_file_path,
            HDBSCAN_min_cluster_size,
            HDBSCAN_min_samples,
            coarsening_dim,
            density_threshold,
            output_dir,
            averaged=True,
            density_coarsening=True,
            precoarsened=True,
            standardize=True,
            threshold_ratio=True,
        ).result()

        centroids, centroid_dir = result
        self.assertTrue(os.path.exists(centroid_dir))
        self.assertTrue(centroids.size > 0)
        self.log_memory_usage(
            "test_convert_density_file_to_point_cloud_hdbscan_timed end"
        )

    def test_extract_graph_features_from_centroids_timed(self):
        """
        Tests the extract_graph_features_from_centroids_timed workflow.
        """
        self.log_memory_usage("test_extract_graph_features_from_centroids_timed begin")

        # Generate 100 random points in a 10x10x10 box
        np.random.seed(42)
        centroids = np.random.rand(100, 3) * 10
        cutoff = 5
        output_dir = self.output_dir

        result = extract_graph_features_from_centroids_timed(
            centroids,
            cutoff,
            output_dir,
            force_overwrite=False,
            skip_clique_num=True,
        ).result()

        self.assertTrue(os.path.exists(result))
        with open(result, "r") as f:
            print(f.read())
        self.log_memory_usage("test_extract_graph_features_from_centroids_timed end")

    def test_full_workflow(self):
        """
        Tests the full workflow of converting a density file to a point cloud, then extracting graph features.
        """
        self.log_memory_usage("test_full_workflow begin")

        protein_name = "test_protein"
        mrc_file_path = Path(self.data_path, "test.mrc")
        HDBSCAN_min_cluster_size = 6
        HDBSCAN_min_samples = 6
        coarsening_dim = 3
        density_threshold = 0.99
        output_dir = self.output_dir
        cutoff = 5

        # Convert density file to point cloud
        centroid_result = convert_density_file_to_point_cloud_hdbscan_timed(
            protein_name,
            mrc_file_path,
            HDBSCAN_min_cluster_size,
            HDBSCAN_min_samples,
            coarsening_dim,
            density_threshold,
            output_dir,
            averaged=True,
            density_coarsening=True,
            precoarsened=True,
            standardize=True,
            threshold_ratio=True,
        )
        centroids = centroid_result[0]
        centroid_dir = centroid_result[1]

        # Extract graph features
        result = extract_graph_features_from_centroids_timed(
            centroids,
            cutoff,
            centroid_dir,
            force_overwrite=False,
            skip_clique_num=True,
        ).result()

        self.assertTrue(os.path.exists(result))
        with open(result, "r") as f:
            print(f.read())
        self.log_memory_usage("test_full_workflow end")


if __name__ == "__main__":
    unittest.main()
