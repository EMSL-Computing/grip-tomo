# August George, 2022, PNNL

import argparse
import mrcfile
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import networkx as nx
from scipy.spatial.distance import pdist,squareform
from pathlib import Path


def load_density_file(fname):
    """load a .mrc file using the mrcfile package

    Args:
        fname ([str]): filename / filepath

    Returns:
        [mrcfile object]: MRC data 
    """
    # load .mrc tomogram file as a MRC object which has header and data properties. 
    # see: https://mrcfile.readthedocs.io/en/latest/usage_guide.html 
    mrc = mrcfile.mmap(fname, mode=u'r')  # memory mapped mode for large files
    return mrc


def normalize_and_threshold_data(mrc, t, noise_stdev=0.0, norm_T=False):
    """normalizes threshold value and densities then applies a cutoff threshold

    Args:
        mrc ([mrcfile object]): mrc data
        t ([float]): raw (unormalized) pixel intensity cutoff threshold
        noise_stdev ([float]): Standard deviation of Gaussian noise (mean=0) to add. Default is 0 (no noise added)
        norm_T ([bool]): threshold value is normalized (True) or not normalized (False). Default is False.

    Returns:
        [numpy array]: array of x,y,z coordinates which are above the cutoff threshold. A[0] = [x0,y0,z0]
    """
    # load and normalize data, normalize threshold value
    if noise_stdev == 0:
        D = mrc.data
    else:
        assert(noise_stdev>=0)
        D = add_Gaussian_noise(mrc,scale=noise_stdev)
           
    D_min = np.min(D)
    D_max = np.max(D)
    D_norm = (D - D_min)/(D_max-D_min)  # normalize to 0,1: (x_i-x_min) / (x_max - x_min)
    if norm_T ==True:
        t_norm = t 
    else:
        t_norm = (t-D_min) / (D_max-D_min)

    # get x,y,z coordinates above threshold
    x,y,z = np.where(D_norm > t_norm)
    xyz_data = np.transpose(np.stack([x,y,z]))
    return xyz_data


def cluster_data(xyz_data, DBSCAN_epsilon, DBSCAN_min_samples):
    """Clusters data using DBSCAN from sklearn

    Args:
        xyz_data ([numpy array]):  A[0] = [x0,y0,z0]
        DBSCAN_epsilon ([float]): DBSCAN epsilon value (in pixels)
        DBSCAN_min_samples ([int]): DBSCAN min_samples

    Returns:
        [sklearn DBSCAN cluster object]: clustering results stored in an object
    """
    model = DBSCAN(eps=DBSCAN_epsilon, min_samples=DBSCAN_min_samples)  # apply coarse-graining (DBSCAN)
    model.fit_predict(xyz_data)
    return model


def get_cluster_centroids(xyz_data, model):
    """Coarse grain density model using cluster centroids

    Args:
        xyz_data ([numpy array]): A[0] = [x0,y0,z0]
        model ([sklearn DBSCAN cluster object]): clustering results stored in an object

    Returns:
        [numpy array]: array of cluster centroids, A[0] = [centroid_x0, centroid_y0, centroid_z0]
    """
    samples_w_lbls = np.concatenate((xyz_data,model.labels_[:,np.newaxis]),axis=1)
    if -1 in set(model.labels_):  # if noise detected
        coarse_model = np.zeros((len(set(model.labels_))-1,3))  # remove last label which is noise
        for i in range(len(set(model.labels_))-1):
        # https://stackoverflow.com/questions/55604239/find-which-points-belong-to-a-cluster-in-dbscan-in-python
            tmp_T = np.transpose(samples_w_lbls[np.in1d(samples_w_lbls[:,-1], np.asarray([i]))])
            x_mean = np.mean(tmp_T[0])
            y_mean = np.mean(tmp_T[1])
            z_mean = np.mean(tmp_T[2])
            coarse_model[i] = [x_mean,y_mean,z_mean]
    else:
        coarse_model = np.zeros((len(set(model.labels_)),3))
        for i in range(len(set(model.labels_))):
            # https://stackoverflow.com/questions/55604239/find-which-points-belong-to-a-cluster-in-dbscan-in-python
            tmp_T = np.transpose(samples_w_lbls[np.in1d(samples_w_lbls[:,-1], np.asarray([i]))])
            x_mean = np.mean(tmp_T[0])
            y_mean = np.mean(tmp_T[1])
            z_mean = np.mean(tmp_T[2])
            coarse_model[i] = [x_mean,y_mean,z_mean]
    return coarse_model


def plot_clustering_results(xyz_data, coarse_model, figsize=3):
    """creates a 3D scatter plot containing both the xyz data and the cluster centroids.

    Note: should rotate afterwards for better visualization.

    Args:
        xyz_data ([numpy array]):  A[0] = [x0,y0,z0]
        coarse_model ([numpy array]): array of cluster centroids, A[0] = [centroid_x0, centroid_y0, centroid_z0]

    Returns:
        [matplotlib figure object]: 3d scatter plot figure
    """
    fig = plt.figure(figsize=(figsize, figsize))
    plt.title('clustering results')
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xyz_data[:,0], xyz_data[:,1], xyz_data[:,2], c='purple', s=2, alpha=0.3)
    ax.scatter(coarse_model[:,0], coarse_model[:,1], coarse_model[:,2], c='k', s=10, alpha=0.9)
    return fig
    

def create_and_save_graph(coarse_model, proximity_px, out_fname, save=True):
    """creates a Networkx graph from the coarse grained model (cluster centroids) and saves it as a graph XML file (.gexf)

    Args:
        coarse_model ([numpy array]): array of cluster centroids, A[0] = [centroid_x0, centroid_y0, centroid_z0]
        proximity_px ([float]): pairwise cutoff distance for assigning edges to nodes, in pixels.
        out_fname ([string]): filename for output
        save ([boolean]): flag to save file (True) or not (False)

    Returns:
        [networkx graph object]: graph representation of coarse model (cluster centroids)
    """
    d_matrix = squareform(pdist(coarse_model, 'euclid'))  
    d_matrix_thresh = np.where(d_matrix>proximity_px, 0, d_matrix)  # if the distance is > t, replace it with 0 (i.e. remove edge)
    G = nx.convert_matrix.from_numpy_matrix(d_matrix_thresh)   # convert distance matrix to networkx graph object
    if save:
        nx.write_gexf(G,f'{out_fname}.gexf')
    return G


def add_Gaussian_noise(mrc, loc=0.0, scale=1.0):
    """Adds Gaussian white noise to the data in an mrc file

    ref: https://numpy.org/doc/stable/reference/random/generated/numpy.random.normal.html 

    Args:
        mrc (mrcfile object): mrc object to add noise to
        loc (float, optional): mean of Gaussian distribution. Defaults to 0.
        scale (float, optional): standard deviation of Gaussian distribution. Defaults to 1.

    Returns:
        [numpy array]: data with noise added 
    """
    D = mrc.data
    noise = np.random.normal(loc=loc, scale=scale, size=D.shape)
    D_w_noise = D + noise 
    return D_w_noise


def main(args):
    """Takes a 3D density (.mrc), applies threshold, coarse-grains data, and converts it into a graph network. 
    Outputs a .png file of the coarse grained model, and a .gexf graph xml file. 

    Args:
        args ([argument parser object]):

            - args.fname: .mrc filename (white density w/ black background)

            - args.t: unormalized pixel intensity threshold level

            - args.eps: DBSCAN epsilon (inter cluster distance)

            - args.ms: DBSCAN min samples (minimum number of samples in cluster)

            - args.d_cut: pairwise distance cutoff for assigning edges to graph, in pixels
    """
    fname = Path(args.fname)
    t = args.t
    DBSCAN_epsilon = args.eps  # try 1
    DBSCAN_min_samples = args.ms  # try 4
    d_cut = args.d_cut  # try 8
    out_fname = fname.with_suffix('')

    mrc = load_density_file(fname)
    xyz_data = normalize_and_threshold_data(mrc,t)
    model = cluster_data(xyz_data,DBSCAN_epsilon,DBSCAN_min_samples)
    coarse_model = get_cluster_centroids(xyz_data,model)
    G = create_and_save_graph(coarse_model,d_cut,out_fname)
    fig = plot_clustering_results(xyz_data,coarse_model)


if __name__ == '__main__':

    # example: >> python density2graph.py fname.mrc 0.5 1 4 8
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", help="tomogram .mrc filename (white density w/ black background)", type=str)
    parser.add_argument("t",  help="pixel intensity threshold cutoff (unormalized)", type=float)
    parser.add_argument("eps", help="DBSCAN epsilon (inter cluster distance) in pixels", type=float)
    parser.add_argument("ms", help="DBSCAN min samples (minimum number of samples in cluster)", type=float)
    parser.add_argument("d_cut",  help="pairwise distance cutoff for assigning edges in pixels", type=float)
    args = parser.parse_args()
    main(args)
