import time
import networkx as nx
import numpy as np
import pandas as pd
from pathlib import Path
import multiprocessing
import igraph as ig


import warnings
import time

warnings.filterwarnings("ignore")


def calc_bc(G, return_dict):
    """
    Parallel subprocess function to calculate the betweenness centrality.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    return_dict : dict
        Dictionary to store the result.

    Returns
    -------
    dict
        Betweenness centrality dictionary from multiple processes.
    """
    return_dict[1] = np.max(list(nx.betweenness_centrality(G).values()))


def calc_shortest_pthlen(G, return_dict):
    """
    Parallel subprocess function to calculate the average shortest path length.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    return_dict : dict
        Dictionary to store the result.

    Returns
    -------
    dict
        Average shortest path length dictionary from multiple processes.
    """
    return_dict[2] = nx.average_shortest_path_length(G)


def benchmark_igraph_graph_features(G):
    """
    Benchmarks several graph network features using igraph. If not connected, largest subgraph is used.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.

    Returns
    -------
    dict
        Dictionary containing the benchmarked features and their computation times.
    """

    # G = ig.Graph.Load(G, format='gexf')
    if type(G) != nx.classes.graph.Graph:
        G = nx.read_gexf(G)
    try:
        assert nx.is_connected(G)
    except:
        # graph isn't connect --> use largest connected component
        # see ref: https://stackoverflow.com/questions/26105764/how-do-i-get-the-giant-component-of-a-networkx-graph
        G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
        G = G.subgraph(G_cc[0])
    """ 
    Maybe rewrite the gexf files to a file format just in this case 
    so igraph can read it? To remove the step of having to load from networkx?
    """
    G_time = {}
    start_time = time.process_time()
    G = ig.Graph.from_networkx(G)
    G_time["conversion"] = time.process_time() - start_time
    G_time["backend"] = "igraph"
    G_feat = {}
    start_time = time.process_time()

    G_feat["n nodes"] = float(ig.Graph.vcount(G))  # number of nodes
    G_time["n nodes"] = time.process_time() - start_time
    G_time["n nodes backend"] = "igraph"

    start_time = time.process_time()
    G_feat["n edges"] = float(ig.Graph.ecount(G))  # number of edges
    G_time["n edges"] = time.process_time() - start_time
    G_time["n edges backend"] = "igraph"

    start_time = time.process_time()
    G_feat["density"] = ig.Graph.density(
        G
    )  # how close the network is to a 'complete graph' where each node is connected
    G_time["density"] = time.process_time() - start_time
    G_time["density backend"] = "igraph"

    start_time = time.process_time()
    G_feat["diameter"] = ig.Graph.diameter(
        G
    )  # the farthest distance (e.g. number of edges) between two nodes in the graph
    G_time["diameter"] = time.process_time() - start_time
    G_time["diameter backend"] = "igraph"

    start_time = time.process_time()
    clustering = ig.Graph.transitivity_local_undirected(G, mode="zero")
    if sum(clustering) == 0 or len(clustering) == 0:
        G_feat["avg clustering"] = sum(clustering)
    else:
        G_feat["avg clustering"] = sum(clustering) / len(clustering)
    G_time["avg clustering"] = time.process_time() - start_time
    G_time["avg clustering backend"] = "igraph"

    start_time = time.process_time()
    G_feat["max closeness centrality"] = np.max(
        list(ig.Graph.closeness(G))
    )  # max closeness centraility. high closeness --> short distance to all other nodes
    G_time["max closeness centrality"] = time.process_time() - start_time
    G_time["max closeness centrality backend"] = "igraph"

    start_time = time.process_time()
    G_feat["max eigenvector centrality"] = np.max(
        ig.Graph.eigenvector_centrality(G, scale=False)
    )  # eigenvector centraility. how 'important'/'influential' a node is
    G_time["max eigenvector centrality"] = time.process_time() - start_time
    G_time["max eigenvector centrality backend"] = "igraph"

    start_time = time.process_time()
    G_feat["degree assortativity"] = (
        ig.Graph.assortativity_degree(G) + 1
    ) / 2  # tendency of a node to be connected with other nodes of the same degree, normalized to 0 and 1.
    G_time["degree assortativity"] = time.process_time() - start_time
    G_time["degree assortativity backend"] = "igraph"

    start_time = time.process_time()
    G_feat["max clique number"] = len(
        max(ig.Graph.largest_cliques(G), key=len)
    )  # largest clique (i.e. an induced subgraph that is complete) size
    G_time["max clique number"] = time.process_time() - start_time
    G_time["max clique number backend"] = "igraph"

    start_time = time.process_time()
    G_feat["n communities"] = ig.Graph.community_fastgreedy(
        G
    ).optimal_count  # number of communities
    G_time["n communities"] = time.process_time() - start_time
    G_time["n communities backend"] = "igraph"

    start_time = time.process_time()
    G_feat["avg path length"] = float(ig.Graph.average_path_length(G))
    G_time["avg path length"] = time.process_time() - start_time
    G_time["avg path length backend"] = "igraph"

    start_time = time.process_time()
    G_feat["max betweenness centrality"] = round(
        np.max(ig.Graph.betweenness(G))
        / (((G_feat["n nodes"] - 1) * (G_feat["n nodes"] - 2)) / 2),
        15,
    )
    G_time["max betweenness centrality"] = time.process_time() - start_time
    G_time["max betweenness centrality backend"] = "igraph"

    return G_time


def deduce_backend(backends_list, backend=None):
    """
    Deduce the backend used for the graph features calculation.

    Parameters
    ----------
    backends_list : list
        List of available backends.
    backend : str, optional
        Preferred backend to use. Default is None.

    Returns
    -------
    str
        Backend to use.
    """
    if backend in backends_list:
        return backend
    else:
        return "serial"


# import nx_cugraph as nxcg
# import nx_parallel as nxp


def benchmark_nx_graph_features(G, backend="serial"):
    """
    Benchmarks several graph network features using networkx. If not connected, largest subgraph is used. Uses multiprocessing for parallelism.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    backend : str, optional
        Backend to use for computation. Default is 'serial'.

    Returns
    -------
    dict
        Dictionary containing the benchmarked features and their computation times.
    """

    try:
        assert nx.is_connected(G)
    except:
        # graph isn't connect --> use largest connected component
        # see ref: https://stackoverflow.com/questions/26105764/how-do-i-get-the-giant-component-of-a-networkx-graph
        G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
        G = G.subgraph(G_cc[0])

    G_time = {}  # graph features time dictionary

    if backend == "cugraph":
        start_time = time.process_time()
        G_conv = nxcg.from_networkx(G)
        G_time["conversion"] = time.process_time() - start_time
    elif backend == "parallel":
        start_time = time.process_time()
        G_conv = nxp.ParallelGraph(G)
        G_time["conversion"] = time.process_time() - start_time
    elif backend == "serial":
        G_conv = G
        G_time["conversion"] = 0
    else:
        raise ValueError("backend must be 'cugraph', 'parallel', or 'serial'")

    G_time["backend"] = backend

    G_feat = {}  # graph features dictionary
    start_time = time.process_time()

    G_feat["n nodes"] = float(G.number_of_nodes())  # number of nodes
    G_time["n nodes"] = time.process_time() - start_time
    print(
        f"Execution time for {backend} : n nodes {G_time['n nodes']} : {G_feat['n nodes']}"
    )
    G_time["n nodes backend"] = "serial"

    start_time = time.process_time()
    G_feat["n edges"] = float(G.number_of_edges())  # number of edges
    G_time["n edges"] = time.process_time() - start_time
    print(
        f"Execution time for {backend} : n edges {G_time['n edges']} : {G_feat['n edges']}"
    )
    G_time["n edges backend"] = "serial"

    start_time = time.process_time()
    G_feat["density"] = nx.density(
        G
    )  # how close the network is to a 'complete graph' where each node is connected
    G_time["density"] = time.process_time() - start_time
    print(
        f"Execution time for {backend} : density {G_time['density']} : {G_feat['density']}"
    )
    G_time["density backend"] = "serial"

    start_time = time.process_time()
    if deduce_backend(nx.diameter.backends, backend) == backend:
        G_feat["diameter"] = nx.diameter(
            G_conv
        )  # the farthest distance (e.g. number of edges) between two nodes in the graph
        G_time["diameter"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : diameter {G_time['diameter']} : {G_feat['diameter']}"
        )
        G_time["diameter backend"] = backend
    elif backend == "cugraph":
        # Calculate shortest path length, find max
        shortest_paths = nxcg.shortest_path(G_conv)
        G_feat["diameter"] = np.max(shortest_paths)
        G_time["diameter"] = time.process_time() - start_time
    else:
        G_feat["diameter"] = 0
        # G_feat['diameter'] = nx.diameter(G)  # the farthest distance (e.g. number of edges) between two nodes in the graph
        G_time["diameter"] = 0
        G_time["diameter backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.average_clustering.backends, backend) == backend:
        G_feat["avg clustering"] = nx.average_clustering(
            G_conv
        )  # the (averaged) fraction of possible triangles through a node.
        G_time["avg clustering"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : avg clustering {G_time['avg clustering']} : {G_feat['avg clustering']}"
        )
        G_time["avg clustering backend"] = backend
    else:
        G_feat["avg clustering"] = 0
        # G_feat['avg clustering'] = nx.average_clustering(G)  # the (averaged) fraction of possible triangles through a node.
        G_time["avg clustering"] = 0
        G_time["avg clustering backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.closeness_centrality.backends, backend) == backend:
        G_feat["max closeness centrality"] = np.max(
            list(nx.closeness_centrality(G_conv).values())
        )  # max closeness centraility. high closeness --> short distance to all other nodes
        G_time["max closeness centrality"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : max closeness centrality {G_time['max closeness centrality']} : {G_feat['max closeness centrality']}"
        )
        G_time["max closeness centrality backend"] = backend
    else:
        G_feat["max closeness centrality"] = 0
        # G_feat['max closeness centrality'] = np.max(list(nx.closeness_centrality(G).values()))  # max closeness centraility. high closeness --> short distance to all other nodes
        G_time["max closeness centrality"] = 0
        G_time["max closeness centrality backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.eigenvector_centrality.backends, backend) == backend:
        G_feat["max eigenvector centrality"] = np.max(
            list(nx.eigenvector_centrality(G_conv, max_iter=10000).values())
        )  # eigenvector centraility. how 'important'/'influential' a node is
        G_time["max eigenvector centrality"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : max eigenvector centrality {G_time['max eigenvector centrality']} : {G_feat['max eigenvector centrality']}"
        )
        G_time["max eigenvector centrality backend"] = backend
    else:
        G_feat["max eigenvector centrality"] = 0
        # G_feat['max eigenvector centrality'] = np.max(list(nx.eigenvector_centrality(G, max_iter=10000).values()))  # eigenvector centraility. how 'important'/'influential' a node is
        G_time["max eigenvector centrality"] = 0
        G_time["max eigenvector centrality backend"] = "N/A"

    start_time = time.process_time()
    if (
        deduce_backend(nx.degree_pearson_correlation_coefficient.backends, backend)
        == backend
    ):
        G_feat["degree assortativity"] = (
            nx.degree_pearson_correlation_coefficient(G_conv) + 1
        ) / 2  # tendency of a node to be connected with other nodes of the same degree, normalized to 0 and 1.
        G_time["degree assortativity"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : degree assortativity {G_time['degree assortativity']} : {G_feat['degree assortativity']}"
        )
        G_time["degree assortativity backend"] = backend
    else:
        G_feat["degree assortativity"] = 0
        # G_feat['degree assortativity'] = (nx.degree_pearson_correlation_coefficient(G) + 1)/2  # tendency of a node to be connected with other nodes of the same degree, normalized to 0 and 1.
        G_time["degree assortativity"] = 0
        G_time["degree assortativity backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.find_cliques.backends, backend) == backend:
        G_feat["max clique number"] = len(
            max(nx.find_cliques(G_conv), key=len)
        )  # largest clique (i.e. an induced subgraph that is complete) size
        G_time["max clique number"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : max clique number {G_time['max clique number']} : {G_feat['max clique number']}"
        )
        G_time["max clique number backend"] = backend
    else:
        G_feat["max clique number"] = 0
        # G_feat['max clique number'] = len(max(nx.find_cliques(G),key=len))  # largest clique (i.e. an induced subgraph that is complete) size
        G_time["max clique number"] = 0
        G_time["max clique number backend"] = "N/A"

    start_time = time.process_time()
    if (
        deduce_backend(
            nx.algorithms.community.modularity_max.greedy_modularity_communities.backends,
            backend,
        )
        == backend
    ):
        G_feat["n communities"] = len(
            nx.algorithms.community.modularity_max.greedy_modularity_communities(G_conv)
        )  # number of communities
        G_time["n communities"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : n communities {G_time['n communities']} : {G_feat['n communities']}"
        )
        G_time["n communities backend"] = backend
    else:
        G_feat["n communities"] = 0
        # G_feat['n communities'] = len(nx.algorithms.community.modularity_max.greedy_modularity_communities(G))  # number of communities
        G_time["n communities"] = 0
        G_time["n communities backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.average_shortest_path_length.backends, backend) == backend:
        G_feat["avg path length"] = nx.average_shortest_path_length(G_conv)
        G_time["avg path length"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : avg path length {G_time['avg path length']} : {G_feat['avg path length']}"
        )
        G_time["avg path length backend"] = backend
    else:
        G_feat["avg path length"] = 0
        # G_feat['avg path length'] = nx.average_shortest_path_length(G)
        G_time["avg path length"] = 0
        G_time["avg path length backend"] = "N/A"

    start_time = time.process_time()
    if deduce_backend(nx.betweenness_centrality.backends, backend) == backend:
        G_feat["max betweenness centrality"] = np.max(
            list(nx.betweenness_centrality(G_conv).values())
        )
        G_time["max betweenness centrality"] = time.process_time() - start_time
        print(
            f"Execution time for {backend} : max betweenness centrality {G_time['max betweenness centrality']} : {G_feat['max betweenness centrality']}"
        )
        G_time["max betweenness centrality backend"] = backend
    else:
        G_feat["max betweenness centrality"] = 0
        # G_feat['max betweenness centrality'] = (np.max(list(nx.betweenness_centrality(G).values())))
        G_time["max betweenness centrality"] = 0
        G_time["max betweenness centrality backend"] = "N/A"
    return G_time


def calc_graph_features(G, backend=None):
    """
    Calculates several graph network features using networkx. If not connected, largest subgraph is used. Uses multiprocessing for parallelism.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    backend : str, optional
        Backend to use for computation. Default is None.

    Returns
    -------
    dict
        Dictionary containing the calculated features.
    """
    try:
        assert nx.is_connected(G)
    except:
        # graph isn't connect --> use largest connected component
        # see ref: https://stackoverflow.com/questions/26105764/how-do-i-get-the-giant-component-of-a-networkx-graph
        G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
        G = G.subgraph(G_cc[0])

    G_feat = {}  # graph features dictionary
    start_time = time.process_time()

    G_feat["n nodes"] = float(G.number_of_nodes())  # number of nodes
    print("Execution time for 'n nodes':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["n edges"] = float(G.number_of_edges())  # number of edges
    print("Execution time for 'n edges':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["density"] = nx.density(
        G
    )  # how close the network is to a 'complete graph' where each node is connected
    print("Execution time for 'density':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["diameter"] = nx.diameter(
        G
    )  # the farthest distance (e.g. number of edges) between two nodes in the graph
    print("Execution time for 'diameter':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["avg clustering"] = nx.average_clustering(
        G
    )  # the (averaged) fraction of possible triangles through a node.
    print("Execution time for 'avg clustering':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["max closeness centrality"] = np.max(
        list(nx.closeness_centrality(G).values())
    )  # max closeness centraility. high closeness --> short distance to all other nodes
    print(
        "Execution time for 'max closeness centrality':",
        time.process_time() - start_time,
    )

    start_time = time.process_time()
    G_feat["max eigenvector centrality"] = np.max(
        list(nx.eigenvector_centrality(G, max_iter=10000).values())
    )  # eigenvector centraility. how 'important'/'influential' a node is
    print(
        "Execution time for 'max eigenvector centrality':",
        time.process_time() - start_time,
    )

    start_time = time.process_time()
    G_feat["degree assortativity"] = (
        nx.degree_pearson_correlation_coefficient(G) + 1
    ) / 2  # tendency of a node to be connected with other nodes of the same degree, normalized to 0 and 1.
    print(
        "Execution time for 'degree assortativity':", time.process_time() - start_time
    )

    start_time = time.process_time()
    G_feat["max clique number"] = len(
        max(nx.find_cliques(G), key=len)
    )  # largest clique (i.e. an induced subgraph that is complete) size
    print("Execution time for 'max clique number':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["n communities"] = len(
        nx.algorithms.community.modularity_max.greedy_modularity_communities(G)
    )  # number of communities
    print("Execution time for 'n communities':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["avg path length"] = nx.average_shortest_path_length(G)
    print("Execution time for 'avg path length':", time.process_time() - start_time)

    start_time = time.process_time()
    G_feat["max betweenness centrality"] = np.max(
        list(nx.betweenness_centrality(G).values())
    )
    print(
        "Execution time for 'max betweenness centrality':",
        time.process_time() - start_time,
    )
    return G_feat


def igraph_calc_graph_features(G):
    """
    Calculates several graph network features using igraph. If not connected, largest subgraph is used.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.

    Returns
    -------
    dict
        Dictionary containing the calculated features.
    """
    # G = ig.Graph.Load(G, format='gexf')
    if type(G) != nx.classes.graph.Graph:
        G = nx.read_gexf(G)
    try:
        assert nx.is_connected(G)
    except:
        # graph isn't connect --> use largest connected component
        # see ref: https://stackoverflow.com/questions/26105764/how-do-i-get-the-giant-component-of-a-networkx-graph
        G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
        G = G.subgraph(G_cc[0])
    """ 
    Maybe rewrite the gexf files to a file format just in this case 
    so igraph can read it? To remove the step of having to load from networkx?
    """
    G = ig.Graph.from_networkx(G)
    G_feat = {}
    G_feat["n nodes"] = float(ig.Graph.vcount(G))
    G_feat["n edges"] = float(ig.Graph.ecount(G))
    G_feat["density"] = ig.Graph.density(G)
    G_feat["diameter"] = ig.Graph.diameter(G)
    clustering = ig.Graph.transitivity_local_undirected(G, mode="zero")
    if sum(clustering) == 0 or len(clustering) == 0:
        G_feat["avg clustering"] = sum(clustering)
    else:
        G_feat["avg clustering"] = sum(clustering) / len(clustering)
    G_feat["max closeness centrality"] = np.max(list(ig.Graph.closeness(G)))
    G_feat["max eigenvector centrality"] = np.max(
        ig.Graph.eigenvector_centrality(G, scale=False)
    )
    G_feat["degree assortativity"] = (
        ig.Graph.assortativity_degree(G) + 1
    ) / 2  # normalize to 0-1
    G_feat["max clique number"] = len(max(ig.Graph.largest_cliques(G), key=len))
    G_feat["n communities"] = ig.Graph.community_fastgreedy(G).optimal_count
    G_feat["avg path length"] = float(ig.Graph.average_path_length(G))
    G_feat["max betweenness centrality"] = round(
        np.max(ig.Graph.betweenness(G))
        / (((G_feat["n nodes"] - 1) * (G_feat["n nodes"] - 2)) / 2),
        15,
    )

    return G_feat


def similarity_measure(x1, x2):
    """
    Calculates the similarity between two feature values.
    Similarity is defined as 1 minus the relative distance between features (x1 and x2).

    Parameters
    ----------
    x1 : float
        Feature value from graph 1 (must range between 0 and 1).
    x2 : float
        Feature value from graph 2 (must range between 0 and 1).

    Returns
    -------
    float
        Relative similarity between the two features.
    """

    return 1 - (np.abs(x1 - x2) / max(x1, x2))


def calc_similarity_score(G1_dict, G2_dict, feature_list):
    """
    Calculates the similarity score of two graphs.

    Parameters
    ----------
    G1_dict : dict or pandas.DataFrame
        Graph 1 features dictionary or dataframe. Must be able to use a key to access values.
    G2_dict : dict or pandas.DataFrame
        Graph 2 features dictionary or dataframe. Must be able to use a key to access values.
    feature_list : list
        List of graph features to compare. Must be keys in the graph features dictionary.

    Returns
    -------
    float
        Similarity score (0 to 1) where 1 indicates identical graphs.
    """
    s_list = []
    for feat in feature_list:
        f1 = G1_dict[feat]
        f2 = G2_dict[feat]
        s_tmp = similarity_measure(f1, f2)
        s_list.append(s_tmp)
    s = np.sum(s_list) / len(s_list)
    return s


def process_graphs(graph_fnames, package=2):
    """
    Takes a list of graph files, calculates their features, and returns them as a dataframe.

    Parameters
    ----------
    graph_fnames : list
        List of graph filenames to process.
    package : int, optional
        Integer value to choose between networkx (1) and igraph (2). Default is 2 (igraph).

    Returns
    -------
    pandas.DataFrame
        Dataframe containing graph features for each graph in the filename list.
    """
    pool = multiprocessing.Pool()
    tmp_feat_list = []
    tmp_df = ""
    if package == 1:  # networkx
        for fname in graph_fnames:
            G_tmp = nx.read_gexf(fname)
            tmp_feat_list.append(calc_graph_features(G_tmp))
        tmp_df = pd.DataFrame(tmp_feat_list)
    elif package == 2:  # igraph
        result = pool.map_async(igraph_calc_graph_features, graph_fnames)
        pool.close()
        pool.join()
        # for fname in graph_fnames:
        # tmp_feat_list.append(igraph_calc_graph_features(fname))
        # print (f"type(result.get()):{type(result.get())}")
        # <class 'list'>
        tmp_df = pd.DataFrame(result.get())
    tmp_df["name"] = [Path(i).stem for i in graph_fnames]
    assert not tmp_df.empty
    return tmp_df


def classify_graphs(class_file_list, sample_file_list, feature_list, package=2):
    """
    Classifies a similarity score from a list of Class and Sample graphs.

    Parameters
    ----------
    class_file_list : list
        List of control/reference graph files (classes).
    sample_file_list : list
        List of non-control/non-reference graph files (samples).
    feature_list : list
        List of features to use for similarity score. Must be valid keys in the graph features dictionary.
    package : int, optional
        Integer value to choose between networkx (1) and igraph (2). Default is 2 (igraph).

    Returns
    -------
    pandas.DataFrame
        Dataframe where each column is a class and each row is the similarity score of the sampled graph.
    """
    if package == 1:
        class_feat_df = process_graphs(class_file_list, 1)
        sample_feat_df = process_graphs(sample_file_list, 1)
    elif package == 2:
        class_feat_df = process_graphs(class_file_list, 2)
        sample_feat_df = process_graphs(sample_file_list, 2)
    class_list = []
    for index, row in sample_feat_df.iterrows():  # for each sample graph
        tmp_graph_feat = row
        class_dict = {}
        class_dict["name"] = tmp_graph_feat["name"]
        for index2, row2 in class_feat_df.iterrows():  # for each class graph
            tmp_class_feat = row2
            s_tmp = calc_similarity_score(tmp_graph_feat, tmp_class_feat, feature_list)
            class_dict[tmp_class_feat["name"]] = s_tmp
        class_list.append(class_dict)
    class_similarity_df = pd.DataFrame(class_list)
    return class_similarity_df


def process_similarity_df(class_similarity_df):
    """
    Generates y_true and y_pred based on the similarity score dataframe.

    y_true is a list where each index is a class and each value is the class value.
    y_pred is a list where each index is a sample and each value is the maximum similarity score for that sample.

    Note
    ----
    This assumes the correct classification is along the diagonal of the similarity matrix/dataframe.

    Parameters
    ----------
    class_similarity_df : pandas.DataFrame
        Dataframe where each column is a class graph and each row is a sample graph. A_ij is the similarity score between graphs i and j. The exception is one column 'name' which contains the names of the sampled graphs for each row.

    Returns
    -------
    tuple of lists
        y_true, y_pred
    """
    num_df = class_similarity_df.drop(columns="name")
    y_true = [i for i in range(len(class_similarity_df.columns) - 1)]
    y_pred = list(num_df.idxmax(axis=0))
    return y_true, y_pred
