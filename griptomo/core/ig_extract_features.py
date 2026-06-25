import igraph as ig
import numpy as np
import time


def igraph_calc_graph_features(G):
    """
    Calculate various graph features for an igraph.Graph object.

    Parameters
    ----------
    G : igraph.Graph
        The input graph.

    Returns
    -------
    dict
        A dictionary containing the calculated graph features:
        - "total n nodes": Total number of nodes in the graph.
        - "total n edges": Total number of edges in the graph.
        - "n nodes": Number of nodes in the largest connected component.
        - "n edges": Number of edges in the largest connected component.
        - "density": Density of the graph.
        - "diameter": Diameter of the graph.
        - "avg clustering": Average clustering coefficient of the graph.
        - "max closeness centrality": Maximum closeness centrality in the graph.
        - "max eigenvector centrality": Maximum eigenvector centrality in the graph.
        - "degree assortativity": Degree assortativity of the graph.
        - "max clique number": Size of the largest clique in the graph.
        - "n communities": Number of communities detected using the fast greedy algorithm.
        - "avg path length": Average path length in the graph.
        - "max betweenness centrality": Maximum betweenness centrality in the graph, normalized by the number of possible node pairs.

    Raises
    ------
    ValueError
        If the input is not an igraph.Graph object.
    """
    if not isinstance(G, ig.Graph):
        raise ValueError("Input must be an igraph.Graph object")

    G_feat = {}
    G_feat["total n nodes"] = G.vcount()
    G_feat["total n edges"] = G.ecount()

    if not G.is_connected():
        # graph isn't connected --> use largest connected component
        # see ref: https://igraph.org/python/doc/igraph.GraphBase-class.html#components
        G = G.clusters().giant()

    G_feat["n nodes"] = float(G.vcount())
    G_feat["n edges"] = float(G.ecount())
    G_feat["density"] = G.density()
    G_feat["diameter"] = G.diameter()
    clustering = G.transitivity_local_undirected(mode="zero")
    if sum(clustering) == 0 or len(clustering) == 0:
        G_feat["avg clustering"] = sum(clustering)
    else:
        G_feat["avg clustering"] = sum(clustering) / len(clustering)
    G_feat["max closeness centrality"] = np.max(list(G.closeness()))
    G_feat["max eigenvector centrality"] = np.max(G.eigenvector_centrality(scale=False))
    G_feat["degree assortativity"] = G.assortativity_degree()
    G_feat["max clique number"] = len(max(G.largest_cliques(), key=len))
    G_feat["n communities"] = G.community_fastgreedy().optimal_count
    G_feat["avg path length"] = float(G.average_path_length())
    G_feat["max betweenness centrality"] = round(
        (
            np.max(G.betweenness())
            / (((G_feat["n nodes"] - 1) * (G_feat["n nodes"] - 2)) / 2)
        ),
        15,
    )

    return G_feat


def igraph_calc_graph_features_timed(G, skip_clique_num=False):
    """
    Calculate various graph features and the time taken to compute each feature for a given igraph.Graph object.

    Parameters
    ----------
    G : igraph.Graph
        The input graph for which features are to be calculated.
    skip_clique_num : bool, optional
        If True, the calculation of the largest clique number will be skipped. Default is False.

    Returns
    -------
    tuple
        A tuple containing two dictionaries:
        - G_feat (dict): A dictionary with keys as feature names and values as the computed feature values.
            Features include:
            - "total n nodes": Total number of nodes in the graph.
            - "total n edges": Total number of edges in the graph.
            - "n nodes": Number of nodes in the largest connected component.
            - "n edges": Number of edges in the largest connected component.
            - "density": Density of the graph.
            - "diameter": Diameter of the graph.
            - "avg clustering": Average clustering coefficient of the graph.
            - "max closeness centrality": Maximum closeness centrality in the graph.
            - "max eigenvector centrality": Maximum eigenvector centrality in the graph.
            - "degree assortativity": Degree assortativity of the graph, normalized between 0 and 1.
            - "max clique number": Size of the largest clique in the graph.
            - "n communities": Number of communities detected using the fast greedy algorithm.
            - "avg path length": Average path length in the graph.
            - "max betweenness centrality": Maximum betweenness centrality in the graph, normalized.

        - G_time (dict): A dictionary with keys as feature names and values as the time, in seconds, taken to compute each feature.
    """
    if not isinstance(G, ig.Graph):
        raise ValueError("Input must be an igraph.Graph object")

    G_feat = {}
    G_time = {}
    start_time = time.perf_counter()
    G_feat["total n nodes"] = G.vcount()
    G_time["total n nodes"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["total n edges"] = G.ecount()
    G_time["total n edges"] = time.perf_counter() - start_time

    if not G.is_connected():
        # graph isn't connected --> use largest connected component
        # see ref: https://igraph.org/python/doc/igraph.GraphBase-class.html#components
        start_time = time.perf_counter()
        G = G.clusters().giant()
        G_time["connected component"] = time.perf_counter() - start_time
    else:
        G_time["connected component"] = 0

    start_time = time.perf_counter()
    G_feat["n nodes"] = float(G.vcount())
    G_time["n nodes"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["n edges"] = float(G.ecount())
    G_time["n edges"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["density"] = G.density()
    G_time["density"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["diameter"] = G.diameter()
    G_time["diameter"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    clustering = G.transitivity_local_undirected(mode="zero")
    G_time["clustering"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    if sum(clustering) == 0 or len(clustering) == 0:
        G_feat["avg clustering"] = sum(clustering)
    else:
        G_feat["avg clustering"] = sum(clustering) / len(clustering)
    G_time["avg clustering"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["max closeness centrality"] = np.max(list(G.closeness()))
    G_time["max closeness centrality"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["max eigenvector centrality"] = np.max(G.eigenvector_centrality(scale=False))
    G_time["max eigenvector centrality"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["degree assortativity"] = (
        G.assortativity_degree() + 1
    ) / 2  # normalized between 0 and 1.
    G_time["degree assortativity"] = time.perf_counter() - start_time

    if not skip_clique_num:
        start_time = time.perf_counter()
        G_feat["max clique number"] = len(max(G.largest_cliques(), key=len))
        G_time["max clique number"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["n communities"] = G.community_fastgreedy().optimal_count
    G_time["n communities"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["avg path length"] = float(G.average_path_length())
    G_time["avg path length"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["max betweenness centrality"] = round(
        (
            np.max(G.betweenness())
            / (((G_feat["n nodes"] - 1) * (G_feat["n nodes"] - 2)) / 2)
        ),
        15,
    )
    G_time["max betweenness centrality"] = time.perf_counter() - start_time

    return G_feat, G_time
