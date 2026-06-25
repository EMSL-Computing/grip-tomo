from __future__ import annotations

import time
from typing import TYPE_CHECKING

import igraph as ig

if TYPE_CHECKING:  # pragma: no cover - hints only evaluated by type-checkers.
    import cudf  # type: ignore
    import cugraph  # type: ignore

try:  # pragma: no cover - GPU-only dependencies may be absent in CPU environments.
    import cugraph as cg  # type: ignore
    import cudf  # type: ignore
    import cupy as cp  # type: ignore
except ImportError as exc:  # pragma: no cover - exercised when RAPIDS is unavailable.
    cg = None  # type: ignore
    cudf = None  # type: ignore
    cp = None  # type: ignore
    _GPU_IMPORT_ERROR: Exception | None = exc
else:
    _GPU_IMPORT_ERROR = None


def cugraph_calc_graph_features(
    G_ig: ig.Graph, G_cg: "cugraph.Graph", skip_clique_num: bool = False
):
    """
    Calculate various graph features for a cugraph.Graph object.

    Parameters
    ----------
    G_ig : igraph.Graph
        The input graph for which features are to be calculated.
    G_cg : cugraph.Graph
        The input graph for which features are to be calculated.
    skip_clique_num : bool, optional
        If True, the calculation of the largest clique number will be skipped. Default is False.

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
        - "degree assortativity": Degree assortativity of the graph (normalized to [0, 1]).
        - "max clique number": Size of the largest clique in the graph.
        - "n communities": Number of communities detected using the louvian algorithm.
        - "avg path length": Average path length in the graph.
        - "max betweenness centrality": Maximum betweenness centrality in the graph, normalized by the number of possible node pairs.

    Raises
    ------
    ValueError
        If the input is not a cugraph.Graph object.

    Notes
    -----
    - The graph is assumed to be undirected and unweighted.
    - If the graph is not connected, the largest connected component is used for feature calculation.
    - The diameter is calculated using a 2-step BFS.
    - GPU implementations are used for all calculations except for clique_num, which is not supported on GPU.
    - n_communities is calculated using the Louvain algorithm instead of the fast greedy algorithm in igraph.
    """
    if cg is None or cudf is None or cp is None:
        raise RuntimeError(
            "cugraph, cudf, and cupy must be installed to use GPU feature extractors"
        ) from _GPU_IMPORT_ERROR

    if not isinstance(G_cg, cg.Graph):
        raise ValueError("Input must be a cugraph.Graph object")

    if not isinstance(G_ig, ig.Graph):
        raise ValueError("Input must be an igraph.Graph object")

    G_feat = {}
    G_feat["total n nodes"] = G_cg.number_of_nodes()
    print("total n nodes", G_feat["total n nodes"])
    G_feat["total n edges"] = G_cg.number_of_edges()
    print("total n edges", G_feat["total n edges"])

    connected_components_df = cg.weakly_connected_components(G_cg)
    if connected_components_df["labels"].nunique() > 1:
        largest_cc = connected_components_df["labels"].value_counts().idxmax()
        G_cg = cg.subgraph(G_cg, connected_components_df["labels"] == largest_cc)
        G_ig = G_ig.clusters().giant()

    G_feat["n nodes"] = float(G_cg.number_of_nodes())
    print("n nodes", G_feat["n nodes"])
    G_feat["n edges"] = float(G_cg.number_of_edges())
    print("n edges", G_feat["n edges"])
    G_feat["density"] = (
        2
        * G_cg.number_of_edges()
        / (G_cg.number_of_nodes() * (G_cg.number_of_nodes() - 1))
    )
    print("density", G_feat["density"])

    # Perform 2-step BFS to calculate diameter
    # Assumed graph is unweighted and undirected
    first_vertex_id = G_cg.nodes().iloc[0]
    print("first_vertex_id", first_vertex_id)
    bfs_result = cg.bfs(G_cg, first_vertex_id, return_predecessors=False)
    furthest_node = bfs_result["vertex"].iloc[
        cp.argmax(cp.array(bfs_result["distance"]))
    ]
    bfs_result_from_furthest = cg.bfs(G_cg, furthest_node, return_predecessors=False)
    G_feat["diameter"] = float(bfs_result_from_furthest["distance"].max())

    # Avg clustering coefficient
    edges_cudf = G_cg.edges()
    src_indices = cp.array(edges_cudf["source"])
    dst_indices = cp.array(edges_cudf["destination"])
    # Don't count selfloops
    not_selfloops = src_indices != dst_indices
    src_indices = src_indices[not_selfloops]
    dst_indices = dst_indices[not_selfloops]
    if src_indices.size == 0:
        G_feat["avg clustering"] = 0
    else:
        degrees = cp.bincount(src_indices, minlength=G_cg.number_of_nodes())
        degrees += cp.bincount(dst_indices, minlength=G_cg.number_of_nodes())
        triangles = cg.triangle_count(G_cg)
        degrees = degrees[triangles["vertex"]]
        counts = cp.array(triangles["counts"])
        denominators = degrees * (degrees - 1)
        clustering = 2 * counts / denominators
        clustering = cp.where(denominators, clustering, 0)
        assert clustering.ndim == 1
        G_feat["avg clustering"] = float(clustering.mean())

    # Max closeness centrality
    max_centrality = 0
    avg_shortest_path = 0
    for vertex_id in G_cg.nodes().values_host:
        shortest_paths = cg.bfs(G_cg, vertex_id)["distance"]
        shortest_path_sum = cp.sum(shortest_paths)
        avg_shortest_path += shortest_path_sum / (
            G_cg.number_of_nodes() * (G_cg.number_of_nodes() - 1)
        )
        closeness = (G_cg.number_of_nodes() - 1) / shortest_path_sum
        if closeness > max_centrality:
            max_centrality = closeness
    G_feat["max closeness centrality"] = max_centrality
    G_feat["max eigenvector centrality"] = cg.eigenvector_centrality(G_cg)[
        "eigenvector_centrality"
    ].max()

    # Degree assortativity
    edges_cudf = G_cg.edges()
    # print(edges_cudf)
    degrees = G_cg.degree()
    degrees_cudf = cudf.merge(
        edges_cudf, degrees, left_on="source", right_on="vertex", how="left"
    )
    degrees_cudf.rename(columns={"degree": "degree_source"}, inplace=True)
    degrees_cudf = cudf.merge(
        degrees_cudf, degrees, left_on="destination", right_on="vertex", how="left"
    )
    degrees_cudf.rename(columns={"degree": "degree_destination"}, inplace=True)
    print(edges_cudf)
    print(degrees_cudf)
    degree_assortativity = cp.corrcoef(
        degrees_cudf["degree_source"], degrees_cudf["degree_destination"]
    )[0, 1]
    G_feat["degree assortativity"] = float((degree_assortativity + 1) / 2)

    if not skip_clique_num:
        raise NotImplementedError(
            "Calculating the largest clique number is not supported on GPU"
        )
        # G_feat["max clique number"] = G_ig.clique_number()

    G_feat["n communities"] = cg.louvain(G_cg)[0]["partition"].nunique()

    G_feat["avg path length"] = avg_shortest_path
    G_feat["max betweenness centrality"] = round(
        cg.betweenness_centrality(G_cg)["betweenness_centrality"].max(), 15
    )

    return G_feat


def cugraph_calc_graph_features_timed(
    G_ig: ig.Graph, G_cg: "cugraph.Graph", skip_clique_num: bool = False
):
    """
    Calculate various graph features and the time taken to compute each feature for a given cugraph.Graph object.

    Parameters
    ----------
    G_ig : igraph.Graph
        The input graph for which features are to be calculated.
    G_cg : cugraph.Graph
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
            - "degree assortativity": Degree assortativity of the graph (normalized to [0, 1]).
            - "max clique number": Size of the largest clique in the graph.
            - "n communities": Number of communities detected using the louvian algorithm.
            - "avg path length": Average path length in the graph.
            - "max betweenness centrality": Maximum betweenness centrality in the graph, normalized by the number of possible node pairs.

        - G_time (dict): A dictionary with keys as feature names and values as the time, in seconds, taken to compute each feature.

    Raises
    ------
    ValueError
        If the input is not a cugraph.Graph object.

    Notes
    -----
    - The graph is assumed to be undirected and unweighted.
    - If the graph is not connected, the largest connected component is used for feature calculation.
    - The diameter is calculated using a 2-step BFS.
    - GPU implementations are used for all calculations except for clique_num, which is not supported on GPU.
    - n_communities is calculated using the Louvain algorithm instead of the fast greedy algorithm in igraph.
    """
    if cg is None or cudf is None or cp is None:
        raise RuntimeError(
            "cugraph, cudf, and cupy must be installed to use GPU feature extractors"
        ) from _GPU_IMPORT_ERROR

    if not isinstance(G_cg, cg.Graph):
        raise ValueError("Input must be a cugraph.Graph object")

    if not isinstance(G_ig, ig.Graph):
        raise ValueError("Input must be an igraph.Graph object")

    G_feat = {}
    G_time = {}
    start_time = time.perf_counter()
    G_feat["total n nodes"] = G_cg.number_of_nodes()
    G_time["total n nodes"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["total n edges"] = G_cg.number_of_edges()
    G_time["total n edges"] = time.perf_counter() - start_time

    connected_components_df = cg.weakly_connected_components(G_cg)
    if connected_components_df["labels"].nunique() > 1:
        start_time = time.perf_counter()
        largest_cc = connected_components_df["labels"].value_counts().idxmax()
        G_cg = cg.subgraph(G_cg, connected_components_df["labels"] == largest_cc)
        G_time["connected component"] = time.perf_counter() - start_time
        G_ig = G_ig.clusters().giant()
    else:
        G_time["connected component"] = 0

    start_time = time.perf_counter()
    G_feat["n nodes"] = float(G_cg.number_of_nodes())
    G_time["n nodes"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["n edges"] = float(G_cg.number_of_edges())
    G_time["n edges"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["density"] = (
        2
        * G_cg.number_of_edges()
        / (G_cg.number_of_nodes() * (G_cg.number_of_nodes() - 1))
    )
    G_time["density"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    first_vertex_id = G_cg.nodes().iloc[0]
    bfs_result = cg.bfs(G_cg, first_vertex_id, return_predecessors=False)
    furthest_node = bfs_result["vertex"].iloc[
        cp.argmax(cp.array(bfs_result["distance"]))
    ]
    bfs_result_from_furthest = cg.bfs(G_cg, furthest_node, return_predecessors=False)
    G_feat["diameter"] = float(bfs_result_from_furthest["distance"].max())
    G_time["diameter"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    edges_cudf = G_cg.edges()
    src_indices = cp.array(edges_cudf["source"])
    dst_indices = cp.array(edges_cudf["destination"])
    # Don't count selfloops
    not_selfloops = src_indices != dst_indices
    src_indices = src_indices[not_selfloops]
    dst_indices = dst_indices[not_selfloops]
    if src_indices.size == 0:
        G_feat["avg clustering"] = 0
    else:
        degrees = cp.bincount(src_indices, minlength=G_cg.number_of_nodes())
        degrees += cp.bincount(dst_indices, minlength=G_cg.number_of_nodes())
        triangles = cg.triangle_count(G_cg)
        degrees = degrees[triangles["vertex"]]
        counts = cp.array(triangles["counts"])
        denominators = degrees * (degrees - 1)
        clustering = 2 * counts / denominators
        clustering = cp.where(denominators, clustering, 0)
        assert clustering.ndim == 1
        G_feat["avg clustering"] = float(clustering.mean())
    G_time["avg clustering"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    max_centrality = 0
    avg_shortest_path = 0
    for vertex_id in G_cg.nodes().values_host:
        shortest_paths = cg.bfs(G_cg, vertex_id)["distance"]
        shortest_path_sum = cp.sum(shortest_paths)
        avg_shortest_path += shortest_path_sum / (
            G_cg.number_of_nodes() * (G_cg.number_of_nodes() - 1)
        )
        closeness = (G_cg.number_of_nodes() - 1) / shortest_path_sum
        if closeness > max_centrality:
            max_centrality = closeness
    G_feat["max closeness centrality"] = max_centrality
    G_time["max closeness centrality"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["max eigenvector centrality"] = cg.eigenvector_centrality(G_cg)[
        "eigenvector_centrality"
    ].max()
    G_time["max eigenvector centrality"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    edges_cudf = G_cg.edges()
    degrees = G_cg.degree()
    degrees_cudf = cudf.merge(
        edges_cudf, degrees, left_on="source", right_on="vertex", how="left"
    )
    degrees_cudf.rename(columns={"degree": "degree_source"}, inplace=True)
    degrees_cudf = cudf.merge(
        degrees_cudf, degrees, left_on="destination", right_on="vertex", how="left"
    )
    degrees_cudf.rename(columns={"degree": "degree_destination"}, inplace=True)
    degree_assortativity = cp.corrcoef(
        degrees_cudf["degree_source"], degrees_cudf["degree_destination"]
    )[0, 1]
    G_feat["degree assortativity"] = float((degree_assortativity + 1) / 2)
    G_time["degree assortativity"] = time.perf_counter() - start_time

    if not skip_clique_num:
        raise NotImplementedError(
            "Calculating the largest clique number is not supported on GPU"
        )
        # G_feat["max clique number"] = G_ig.clique_number()
        # G_time["max clique number"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["n communities"] = cg.louvain(G_cg)[0]["partition"].nunique()
    G_time["n communities"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["avg path length"] = avg_shortest_path
    G_time["avg path length"] = time.perf_counter() - start_time

    start_time = time.perf_counter()
    G_feat["max betweenness centrality"] = round(
        cg.betweenness_centrality(G_cg)["betweenness_centrality"].max(), 15
    )
    G_time["max betweenness centrality"] = time.perf_counter() - start_time

    return G_feat, G_time
