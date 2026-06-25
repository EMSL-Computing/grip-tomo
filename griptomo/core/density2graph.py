# August George, 2025, PNNL

import argparse
import mrcfile
import matplotlib.pyplot as plt

try:
    import cupy as cp
    from cuml.cluster import DBSCAN as DBSCAN_GPU, HDBSCAN as HDBSCAN_GPU
    from cuml.preprocessing import StandardScaler as StandardScaler_GPU
    import cugraph as cg
    import cudf

    if cp.cuda.is_available():
        IS_GPU_AVAILABLE = True
        print("GPU support for GRIP-Tomo clustering is available.")
    else:
        raise Exception("No GPU detected.")
except Exception as e:
    print(
        f"GPU support for GRIP-Tomo clustering is not available because of {e}. Using CPU."
    )
    IS_GPU_AVAILABLE = False
import numpy as np
from hdbscan import HDBSCAN
from sklearn.cluster import DBSCAN, OPTICS
from sklearn.preprocessing import StandardScaler
import networkx as nx
from pathlib import Path
from warnings import warn
from scipy.spatial import KDTree
import igraph as ig


def _ensure_array(data, device=None):
    """Ensures that the input data is a numpy array."""
    if isinstance(data, np.ndarray):
        if device == "gpu" and IS_GPU_AVAILABLE:
            return cp.asarray(data)
        return data
    elif isinstance(data, mrcfile.mrcfile.MrcFile):
        array_data = np.asarray(data.data)
        if device == "gpu" and IS_GPU_AVAILABLE:
            return cp.asarray(array_data)
        return array_data
    if IS_GPU_AVAILABLE:
        if isinstance(data, cp.ndarray):
            return data
    else:
        raise ValueError(f"Unsupported data type: {type(data)}")


def load_density_file(fname, print_warning=True):
    """
    Load a .mrc file using the mrcfile package.

    Parameters
    ----------
    fname : str
        Filename or filepath of the .mrc file.
    print_warning : bool, optional
        If True, prints a warning message about data inversion. Default is True.

    Returns
    -------
    mrcfile.mrcfile.MrcFile
        MRC file object containing the header and data properties.
    """
    # load .mrc tomogram file as a MRC object which has header and data properties.
    # see: https://mrcfile.readthedocs.io/en/latest/usage_guide.html
    mrc = mrcfile.mmap(fname, mode="r")  # memory mapped mode for large files
    # warn('WARNING: check the inversion of the data. White voxels should correspond to high density.')
    if print_warning:
        print(
            "WARNING: check the inversion of the data. White voxels should correspond to high density."
        )
    return mrc


def normalize_and_threshold_data(mrc, mrc_t, noise_stdev=0.0, norm_T=True):
    """
    Normalizes threshold value and densities then applies a cutoff threshold.

    Parameters
    ----------
    mrc : mrcfile.mrcfile.MrcFile
        MRC data.
    mrc_t : float
        Raw (unnormalized) pixel intensity cutoff threshold.
    noise_stdev : float, optional
        Standard deviation of Gaussian noise (mean=0) to add. Default is 0 (no noise added).
    norm_T : bool, optional
        If True, threshold value is normalized. Default is True.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates which are above the cutoff threshold. A[0] = [x0, y0, z0].
    """
    warn("This method is deprecated.", DeprecationWarning, stacklevel=2)
    # load and normalize data, normalize threshold value
    if noise_stdev == 0:
        D = mrc.data
    else:
        assert noise_stdev >= 0
        print("add_Gaussian_noise")
        D = add_Gaussian_noise(mrc, scale=noise_stdev)

    # get x,y,z coordinates above threshold
    if norm_T == False:
        print("Do not normalize")
        x, y, z = np.where(D > mrc_t)
    else:
        D_min = np.min(D)
        D_max = np.max(D)
        D_norm = (D - D_min) / (
            D_max - D_min
        )  # normalize to 0,1: (x_i-x_min) / (x_max - x_min)

        t_norm = (mrc_t - D_min) / (D_max - D_min)
        # Doo Nam confirmed that normalization of density ran well.
        # Doo Nam confirmed that as long as he puts reasonable mrc_t, t_norm is 0~1
        # However, since both density and t_norm are normalized, there was no practice difference for xyz_data (regardless of norm_T)

        x, y, z = np.where(D_norm > t_norm)
    xyz_data = np.transpose(np.stack([x, y, z]))
    return xyz_data


### end of def normalize_and_threshold_data(mrc, mrc_t, noise_stdev=0.0, norm_T=True):


def normalize_mrc_data(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray,
    lower_bound: float = -1,
    upper_bound: float = 1,
    device: str = "cpu",
) -> np.ndarray:
    """Normalizes the data of an MRC file.

    Normalizes the data from `mrc` to the range of
    [`lower_bound`, `upper_bound`].

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray or cupy.ndarray
        MRC file or MRC data (numpy array or cupy array).
    lower_bound : float, optional
        Lower bound of the normalized data. Default is -1.
    upper_bound : float, optional
        Upper bound of the normalized data. Default is 1.
    device : str, optional
        Device to use for computation. Default is cpu.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Array of normalized MRC data.

    Notes
    -----
    x_norm = lower_bound + (x - min(x))(upper_bound-lower_bound)/(max(x)-min(x))
    """
    D = _ensure_array(mrc_data, device)
    assert lower_bound < upper_bound
    D_min = np.min(D)
    D_max = np.max(D)
    D_norm = lower_bound + (upper_bound - lower_bound) * (D - D_min) / (D_max - D_min)
    return D_norm


def standardize_mrc_data(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray, device: str = "cpu"
) -> np.ndarray:
    """Standardizes the data using scikit-learn's StandardScaler.

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray or cupy.ndarray
        MRC file or MRC data (numpy array or cupy array).
    device : str, optional
        Device to use for computation. Default is cpu.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Array of standardized MRC data.

    Notes
    -----
    Standardization: z = (x - u) / s
    u: mean of the data
    s: standard deviation of the data

    """
    D = _ensure_array(mrc_data, device)

    if device == "gpu" and IS_GPU_AVAILABLE:
        scaler = StandardScaler_GPU()
    else:
        scaler = StandardScaler()
    D_standardized = scaler.fit_transform(D.reshape(-1, 1)).reshape(D.shape)

    return D_standardized


def threshold_mrc_data(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray,
    threshold: float,
    device: str = "cpu",
) -> np.ndarray:
    """
    Keeps only the MRC data above the threshold value, setting the other data to the min value.

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray or cupy.ndarray
        MRC file or MRC data (numpy array or cupy array).
    threshold : float
        Threshold value.
    device : str, optional
        Device to use for computation. Default is cpu.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Thresholded MRC data.

    Notes
    ------
    This function is agnostic towards data normalization.
    """
    D = _ensure_array(mrc_data, device)
    thresholded_D = np.where(D > threshold, D, np.min(D))
    return thresholded_D


def identify_threshold_ratio(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray, ratio: float, device: str = "cpu"
) -> np.ndarray:
    """
    Given the MRC file, finds the density value that corresponds to the requested threshold ratio from 0 to 1.

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray or cupy.ndarray
        MRC file or MRC data (numpy array or cupy array).
    ratio : float
        Threshold ratio from 0 to 1. 1 = max density (remove all points), 0 = min density (keep all points).
    device : str, optional
        Device to use for computation. Default is cpu.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Density value corresponding to the requested threshold ratio.
    """
    assert 0 <= ratio <= 1
    D = _ensure_array(mrc_data, device)
    sorted_densities = np.sort(D.flatten())
    threshold_index = int(round(ratio * len(sorted_densities)))
    if threshold_index == len(sorted_densities):
        raise ValueError(
            f"Threshold index is at the end of the sorted densities. Please choose a lower ratio."
        )
    return sorted_densities[threshold_index]


def generate_point_cloud_from_mrc_data(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray,
    threshold: float,
    device: str = "cpu",
) -> np.ndarray:
    """
    Generates an array of X, Y, and Z coordinates of the data that is above a threshold.

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray or cupy.ndarray
        MRC file or MRC data (numpy array or cupy array).
    threshold : float
        Threshold value.
    device : str, optional
        Device to use for computation. Default is cpu.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Array of X, Y, and Z coordinates of the data that is above the threshold.

    Notes
    -----
    This function is agnostic towards data normalization.
    """
    D = _ensure_array(mrc_data, device)
    x, y, z = np.where(D > threshold)
    xyz_data = np.transpose(np.stack([x, y, z]))
    return xyz_data


def augment_mrc_data(
    mrc_data: mrcfile.mrcfile.MrcFile | np.ndarray,
    offset_percent: float = 10,
    n: int = 1,
) -> list[np.ndarray] | np.ndarray:
    """
    Generates N augmented MRC datasets that are injected with random noise from a uniform distribution bounded by the offset percentage.

    Parameters
    ----------
    mrc_data : mrcfile.mrcfile.MrcFile or numpy.ndarray
        MRC file or MRC data (numpy array).
    offset_percent : float, optional
        Offset percentage for random noise. Default is 10.
    n : int, optional
        Number of augmented datasets to generate. Default is 1.

    Returns
    -------
    list[numpy.ndarray] or numpy.ndarray
        List of augmented MRC datasets or a single augmented dataset.
    """
    D = _ensure_array(mrc_data)
    assert 0.0 <= offset_percent <= 100
    assert n >= 1

    random_uniform = np.random.uniform

    if n == 1:
        random_offsets = (
            random_uniform(low=-offset_percent, high=offset_percent, size=D.shape) / 100
        )
        adjusted_data = D * (1 + random_offsets)
        return adjusted_data
    else:
        adjusted_data_list = []
        for _ in range(n):
            random_offsets = (
                random_uniform(low=-offset_percent, high=offset_percent, size=D.shape)
                / 100
            )
            adjusted_data = D * (1 + random_offsets)
            adjusted_data_list.append(adjusted_data)
        return adjusted_data_list


def calculate_point_cloud_mass(
    point_cloud: np.ndarray, voxel_dim: float, density: float
) -> float:
    """
    Calculates the mass of an array of voxel xyz coordinates.

    Parameters
    ----------
    point_cloud : numpy.ndarray or cupy.ndarray
        Array of voxel xyz coordinates.
    voxel_dim : float
        Voxel dimension.
    density : float
        Density value.

    Returns
    -------
    float
        Estimated mass of the point cloud.

    Notes
    -----
    The estimated mass = n_voxels * voxel_dim * density.
    For example, mass in kDa = 10 voxels * 1 nm3/ voxel * 0.82 kDa / nm3.
    The units for `voxel_dim` and `density` must be consistent.
    Use 0.82 kDa / nm3 for large globular proteins.
    """
    n_voxels = point_cloud.shape[0]
    return n_voxels * voxel_dim * density


def cluster_data(xyz_data, DBSCAN_epsilon, DBSCAN_min_samples):
    """
    Clusters data using DBSCAN from sklearn or cuml.

    Parameters
    ----------
    xyz_data : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    DBSCAN_epsilon : float
        The maximum distance between two samples for one to be considered as in the neighborhood of the other.
    DBSCAN_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.

    Returns
    -------
    sklearn.cluster.DBSCAN or cuml.cluster.DBSCAN: clustering results stored in an object

    See Also
    -------
    sklearn.cluster.DBSCAN
    cuml.cluster.DBSCAN
    cluster_data_hdbscan
    cluster_data_optics
    """
    if IS_GPU_AVAILABLE and isinstance(xyz_data, cp.ndarray):
        clusterer = DBSCAN_GPU
    else:
        clusterer = DBSCAN
    try:
        model = clusterer(
            eps=DBSCAN_epsilon, min_samples=DBSCAN_min_samples
        )  # apply coarse-graining (DBSCAN)
        model.fit_predict(xyz_data)
        return model
    except:
        print(
            'Perhaps, "Fatal Python error: _PyErr_NormalizeException: Cannot recover from MemoryErrors while normalizing exceptions. Python runtime state: initialized" at tahoma for empiar inverted mrc'
        )
        print(f"Failed xyz_data.shape:{xyz_data.shape}")
        print(f"DBSCAN_epsilon:{DBSCAN_epsilon}")
        print(f"DBSCAN_min_samples:{DBSCAN_min_samples}")
        return False


def cluster_data_hdbscan(xyz_data, HDBSCAN_min_cluster_size, HDBSCAN_min_samples):
    """
    Clusters data using HDBSCAN from hdbscan or cuml.

    Parameters
    ----------
    xyz_data : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    HDBSCAN_min_cluster_size : int
        The minimum number of points in a cluster.
    HDBSCAN_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.

    Returns
    -------
    hdbscan.HDBSCAN or cuml.cluster.HDBSCAN: clustering results stored in an object

    See Also
    -------
    hdbscan.HDBSCAN
    cuml.cluster.HDBSCAN
    cluster_data
    cluster_data_optics
    """
    if IS_GPU_AVAILABLE and isinstance(xyz_data, cp.ndarray):
        clusterer = HDBSCAN_GPU
    else:
        clusterer = HDBSCAN
    try:
        model = clusterer(
            min_cluster_size=HDBSCAN_min_cluster_size,
            min_samples=HDBSCAN_min_samples,
            allow_single_cluster=False,
            core_dist_n_jobs=-1,
            approx_min_span_tree=False,
            cluster_selection_method="eom",
        )  # apply clustering (HDBSCAN)
        model.fit(xyz_data)
        return model
    except Exception as e:
        print(f"Error while performing HDBSCAN: {e}")
        print(f"Failed xyz_data.shape:{xyz_data.shape}")
        print(f"HDBSCAN_min_cluster_size:{HDBSCAN_min_cluster_size}")
        print(f"HDBSCAN_min_samples:{HDBSCAN_min_samples}")
        return False


def cluster_data_optics(
    xyz_data,
    OPTICS_min_samples,
    OPTICS_min_cluster_size=None,
    OPTICS_max_epsilon=np.inf,
):
    """
    Clusters data using OPTICS from sklearn.

    Parameters
    ----------
    xyz_data : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    OPTICS_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.
    OPTICS_min_cluster_size : int, optional
        The minimum number of points in a cluster. Default is None.
    OPTICS_max_epsilon : float, optional
        The maximum distance between points considered by OPTICS. Default is numpy.inf.

    Returns
    -------
    sklearn.cluster.OPTICS: clustering results stored in an object

    See Also
    -------
    sklearn.cluster.OPTICS
    cluster_data_hdbscan
    cluster_data
    """
    try:
        model = OPTICS(
            min_cluster_size=OPTICS_min_cluster_size,
            min_samples=OPTICS_min_samples,
            max_eps=OPTICS_max_epsilon,
            n_jobs=-1,
        )  # apply clustering (OPTICS)
        model.fit(xyz_data)
        return model
    except Exception as e:
        print(f"Error while performing OPTICS: {e}")
        print(f"Failed xyz_data.shape:{xyz_data.shape}")
        print(f"OPTICS_min_cluster_size:{OPTICS_min_cluster_size}")
        print(f"OPTICS_min_samples:{OPTICS_min_samples}")
        return False


def voxel_coarsening(voxel_dim, point_cloud, labels, density=None, averaged=False):
    """
    Coarsens largest cluster into voxel_dim^3 chunks.

    Parameters
    ----------
    voxel_dim : int
        Chunk dimension for coarsening.
    point_cloud : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    labels : numpy.ndarray or cupy.ndarray
        Array of cluster labels.
    density : numpy.ndarray or cupy.ndarray, optional
        Array of density values. Default is None.
    averaged : bool, optional
        If True, the coarsened point cloud will be the average of the points in each chunk, using the density as weights, if provided. Default is False.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Coordinates of the centroid of each chunk that contain more than half of the non-empty voxels.
    """
    max_x_dim, max_y_dim, max_z_dim = np.max(point_cloud, axis=0) + 1
    coarsened_dim_x = int(np.ceil(max_x_dim / voxel_dim))
    coarsened_dim_y = int(np.ceil(max_y_dim / voxel_dim))
    coarsened_dim_z = int(np.ceil(max_z_dim / voxel_dim))
    largest_cluster = np.argmax(np.bincount(labels[labels != -1]))
    dense_bins_count = np.zeros((coarsened_dim_x, coarsened_dim_y, coarsened_dim_z))
    label_mask = labels == largest_cluster
    isolated_points = point_cloud[label_mask]
    coarsened_points = isolated_points // voxel_dim
    coarsened_points = coarsened_points.astype(int)
    if density is not None:
        isolated_density = density[label_mask]
        coarse_dense_count = np.zeros_like(dense_bins_count)
    if not averaged:
        for point in coarsened_points:
            dense_bins_count[point[0], point[1], point[2]] += 1
        x, y, z = np.where(dense_bins_count >= (voxel_dim**3) / 2)
        coarsened_point_cloud = np.transpose(np.stack([x, y, z])) * 3 + 1
    else:
        dense_bins = np.zeros((coarsened_dim_x, coarsened_dim_y, coarsened_dim_z, 3))
        if density is not None:
            np.add.at(
                dense_bins,
                (
                    coarsened_points[:, 0],
                    coarsened_points[:, 1],
                    coarsened_points[:, 2],
                ),
                isolated_points * isolated_density[:, np.newaxis],
            )
            np.add.at(
                dense_bins_count,
                (
                    coarsened_points[:, 0],
                    coarsened_points[:, 1],
                    coarsened_points[:, 2],
                ),
                isolated_density,
            )
            np.add.at(
                coarse_dense_count,
                (
                    coarsened_points[:, 0],
                    coarsened_points[:, 1],
                    coarsened_points[:, 2],
                ),
                1,
            )
        else:
            for coarse_point, raw_point in zip(coarsened_points, isolated_points):
                dense_bins[coarse_point[0], coarse_point[1], coarse_point[2]] += (
                    raw_point
                )
                dense_bins_count[coarse_point[0], coarse_point[1], coarse_point[2]] += 1
        dense_bin_points = np.where(dense_bins_count > 0)
        coarsened_point_cloud = (
            dense_bins[dense_bin_points]
            / dense_bins_count[dense_bin_points][:, np.newaxis]
        )
        if density is not None:
            coarsened_density = (
                dense_bins_count[dense_bin_points]
                / coarse_dense_count[dense_bin_points]
            )
            return coarsened_point_cloud, coarsened_density
    return coarsened_point_cloud


def get_cluster_centroids(xyz_data, model):
    """
    Coarse grain density model using cluster centroids.

    Parameters
    ----------
    xyz_data : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    model : sklearn.cluster.DBSCAN or cuml.cluster.DBSCAN
        DBSCAN clustering results stored in an object.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Array of cluster centroids. A[0] = [centroid_x0, centroid_y0, centroid_z0].
    """
    samples_w_lbls = np.concatenate((xyz_data, model.labels_[:, np.newaxis]), axis=1)
    if -1 in set(model.labels_):  # if noise detected
        coarse_model = np.zeros(
            (len(set(model.labels_)) - 1, 3)
        )  # remove last label which is noise
        for i in range(len(set(model.labels_)) - 1):
            # https://stackoverflow.com/questions/55604239/find-which-points-belong-to-a-cluster-in-dbscan-in-python
            tmp_T = np.transpose(
                samples_w_lbls[np.in1d(samples_w_lbls[:, -1], np.asarray([i]))]
            )
            x_mean = np.mean(tmp_T[0])
            y_mean = np.mean(tmp_T[1])
            z_mean = np.mean(tmp_T[2])
            coarse_model[i] = [x_mean, y_mean, z_mean]
    else:
        coarse_model = np.zeros((len(set(model.labels_)), 3))
        for i in range(len(set(model.labels_))):
            # https://stackoverflow.com/questions/55604239/find-which-points-belong-to-a-cluster-in-dbscan-in-python
            tmp_T = np.transpose(
                samples_w_lbls[np.in1d(samples_w_lbls[:, -1], np.asarray([i]))]
            )
            x_mean = np.mean(tmp_T[0])
            y_mean = np.mean(tmp_T[1])
            z_mean = np.mean(tmp_T[2])
            coarse_model[i] = [x_mean, y_mean, z_mean]
    return coarse_model


def plot_clustering_results(xyz_data, coarse_model, figsize=3):
    """
    Creates a 3D scatter plot containing both the xyz data and the cluster centroids.

    Note: should rotate afterwards for better visualization.

    Parameters
    ----------
    xyz_data : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    coarse_model : numpy.ndarray or cupy.ndarray
        Array of cluster centroids. A[0] = [centroid_x0, centroid_y0, centroid_z0].
    figsize : int, optional
        Size of the figure. Default is 3.

    Returns
    -------
    matplotlib.figure.Figure
        3D scatter plot figure.
    """
    fig = plt.figure(figsize=(figsize, figsize))
    plt.title("clustering results")
    ax = fig.add_subplot(projection="3d")
    ax.scatter(
        xyz_data[:, 0], xyz_data[:, 1], xyz_data[:, 2], c="purple", s=2, alpha=0.3
    )
    ax.scatter(
        coarse_model[:, 0],
        coarse_model[:, 1],
        coarse_model[:, 2],
        c="k",
        s=10,
        alpha=0.9,
    )
    return fig


def create_cugraph_from_point_cloud(point_cloud, proximity_px):
    """
    Creates a cuGraph object from a point cloud with edges assigned based on proximity.

    Parameters
    ----------
    point_cloud : numpy.ndarray or cupy.ndarray
        Array of x, y, z coordinates.
    proximity_px : float
        Pairwise cutoff distance for assigning edges to nodes, in pixels.
    savepath : str, optional
        Path to save the cuGraph object as a GML file. If None, don't save.
        Default is None.

    Returns
    -------
    cugraph.Graph
        Unweighted, undirected cuGraph graph object.

    Notes
    -----
    The cuGraph graph is unweighted and undirected.
    Uses the KDTree algorithm to find the neighbors under the cutoff distance.
    """
    if IS_GPU_AVAILABLE:
        if isinstance(point_cloud, cp.ndarray):
            point_cloud = cp.asnumpy(point_cloud)
        tree = KDTree(point_cloud)
        edgelist = cp.array(tree.query_pairs(proximity_px, output_type="ndarray"))
        edgelist_cudf = cudf.DataFrame(edgelist, columns=["source", "destination"])
        edgelist_cudf = cg.symmetrize_df(edgelist_cudf, "source", "destination")
        G = cg.Graph()
        G.from_cudf_edgelist(edgelist_cudf, source="source", destination="destination")
        return G
    else:
        raise ValueError("GPU support is not available.")


def create_igraph_from_point_cloud(point_cloud, proximity_px, savepath=None):
    """
    Creates an igraph object from a point cloud with edges assigned based on proximity.

    Parameters
    ----------
    point_cloud : numpy.ndarray
        Array of x, y, z coordinates.
    proximity_px : float
        Pairwise cutoff distance for assigning edges to nodes, in pixels.
    savepath : str, optional
        Path to save the igraph object as a GML file. If None, don't save.
        Default is None.

    Returns
    -------
    igraph.Graph
        Unweighted, undirected igraph graph object.

    Notes
    -----
    The igraph graph is unweighted and undirected.
    Uses the KDTree algorithm to find the neighbors under the cutoff distance.
    """
    tree = KDTree(point_cloud)
    edgelist = tree.query_pairs(proximity_px, output_type="ndarray")
    G = ig.Graph(edgelist)
    if savepath:
        G.write(savepath, format="gml")
    return G


def create_and_save_graph(coarse_model, proximity_px, out_fname, save=True):
    """
    Creates a Networkx graph from the coarse grained model (cluster centroids) and saves it as a graph XML file (.gexf).

    Parameters
    ----------
    coarse_model : numpy.ndarray or cupy.ndarray
        Array of cluster centroids. A[0] = [centroid_x0, centroid_y0, centroid_z0].
    proximity_px : float
        Pairwise cutoff distance for assigning edges to nodes, in pixels.
    out_fname : str
        Filename for output.
    save : bool, optional
        Flag to save file (True) or not (False). Default is True.

    Returns
    -------
    networkx.Graph
        Networkx graph representation of coarse model (cluster centroids).
    """
    tree = KDTree(coarse_model)
    edgelist = tree.query_pairs(proximity_px, output_type="ndarray")

    # Arsam Firoozfar "It appears that the current code using nx is slightly faster than using igraph to directly create graph from matrix." 3/31/2023

    G = nx.from_edgelist(edgelist)

    if save:
        nx.write_gexf(G, f"{out_fname}.gexf")
    return G


def add_Gaussian_noise(mrc, loc=0.0, scale=1.0):
    """
    Adds Gaussian white noise to the data in an mrc file.

    Parameters
    ----------
    mrc : mrcfile.mrcfile.MrcFile
        MRC object to add noise to.
    loc : float, optional
        Mean of Gaussian distribution. Default is 0.
    scale : float, optional
        Standard deviation of Gaussian distribution. Default is 1.

    Returns
    -------
    numpy.ndarray or cupy.ndarray
        Data with noise added.
    """
    D = mrc.data
    noise = np.random.normal(loc=loc, scale=scale, size=D.shape)
    D_w_noise = D + noise
    return D_w_noise


def leave_2D_density(mrc_data, middle_z_index, xyz_data):
    """
    Extract a 2D slice from 3D density data at a specified z-index.

    This function creates a binary mask of the original 3D density array,
    preserving voxels at the specified z-index plane and zeroing all others.
    Used for 2D projection analysis or slice-based feature extraction.

    Parameters
    ----------
    mrc_data : numpy.ndarray
        3D array of density values from MRC file.
    middle_z_index : int
        Z-axis index of the slice to preserve (0-based indexing).
    xyz_data : numpy.ndarray
        Nx3 array of [x, y, z] coordinates indicating voxels above threshold.

    Returns
    -------
    numpy.ndarray
        3D binary array with same shape as mrc_data, where voxels at
        middle_z_index matching xyz_data coordinates are set to 1, all others 0.

    Notes
    -----
    Primarily used for debugging or visualization of 2D slices from 3D volumes.
    The function name may be misleading - it creates a binary mask rather than
    extracting density values.

    TODO: Consider renaming to `create_2D_slice_mask` for clarity.
    """
    three_D_ndarray_left = np.zeros(
        (mrc_data.shape[0], mrc_data.shape[1], mrc_data.shape[2])
    )
    for k in range(xyz_data.shape[0]):
        # print (xyz_data[k]) # 0 0 64
        # print (type(xyz_data[k])) #<class 'numpy.ndarray'>
        x_to_keep = xyz_data[k][0]  # 0
        y_to_keep = xyz_data[k][1]  # 0
        z_to_keep = xyz_data[k][2]  # 64
        three_D_ndarray_left[x_to_keep, y_to_keep, z_to_keep] = mrc_data[
            x_to_keep, y_to_keep, z_to_keep
        ]
    return three_D_ndarray_left


def main(args):
    """
    Takes a 3D density (.mrc), applies threshold, coarse-grains data, and converts it into a graph network.
    Outputs a .png file of the coarse grained model, and a .gexf graph xml file.

    Parameters
    ----------
    args : argparse.Namespace
        Argument parser object containing the following attributes:
        - fname : str
            .mrc filename (white density w/ black background).
        - mrc_t : float
            Unnormalized pixel intensity threshold level.
        - eps : float
            DBSCAN epsilon (inter cluster distance).
        - ms : int
            DBSCAN min samples (minimum number of samples in cluster).
        - d_cut : float
            Pairwise distance cutoff for assigning edges to graph, in pixels.
    """
    fname = Path(args.fname)
    mrc_t = args.mrc_t
    DBSCAN_epsilon = args.eps  # try 1
    DBSCAN_min_samples = args.ms  # try 4
    d_cut = args.d_cut  # try 8
    out_fname = fname.with_suffix("")

    mrc = load_density_file(fname)
    xyz_data = normalize_and_threshold_data(mrc, mrc_t)

    # print (f"xyz_data:{xyz_data}")

    model = cluster_data(xyz_data, DBSCAN_epsilon, DBSCAN_min_samples)
    coarse_model = get_cluster_centroids(xyz_data, model)
    G = create_and_save_graph(coarse_model, d_cut, out_fname)
    fig = plot_clustering_results(xyz_data, coarse_model)


if __name__ == "__main__":
    # example: >> python density2graph.py fname.mrc 0.5 1 4 8
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fname",
        help="tomogram .mrc filename (white density w/ black background)",
        type=str,
    )
    parser.add_argument(
        "t", help="pixel intensity threshold cutoff (unormalized)", type=float
    )
    parser.add_argument(
        "eps", help="DBSCAN epsilon (inter cluster distance) in pixels", type=float
    )
    parser.add_argument(
        "ms", help="DBSCAN min samples (minimum number of samples in cluster)", type=int
    )
    parser.add_argument(
        "d_cut",
        help="pairwise distance cutoff for assigning edges in pixels",
        type=float,
    )
    args = parser.parse_args()
    main(args)
