import networkx as nx
import numpy as np
import pandas as pd
from pathlib import Path
import multiprocessing


def calc_bc(G, return_dict):
    """Parallel subprocess function to calculate the betweenness centrality.

    Args:
        G ([networkx graph]): graph

    Returns:
        [dictionary]: betweeness centrality dictionary from multiple processes
    """
    return_dict[1] = (np.max(list(nx.betweenness_centrality(G).values())))


def calc_shortest_pthlen(G, return_dict):
    """Parallel subprocess function to calculate the average shortest path length.

    Args:
        G ([networkx graph]): graph

    Returns:
        [dictionary]: average shortest path length dictionary from multiple processes
    """
    return_dict[2] = nx.average_shortest_path_length(G)


def calc_graph_features(G):
    """Calculates several graph network features. If not connected, largest subgraph is used. Uses multiprocessing for parallelsim.

    Args:
        G ([networkx graph]): graph

    Returns:
        [dictionary]: features dictionary
    """
    try: 
        assert(nx.is_connected(G))  
    except:
        # graph isn't connect --> use largest connected component
        # see ref: https://stackoverflow.com/questions/26105764/how-do-i-get-the-giant-component-of-a-networkx-graph
        G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
        G = G.subgraph(G_cc[0])

    # manager for sharing the dictionary that will store the return values
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    # call and pass the values to functions that are in different processes
    process1 = multiprocessing.Process(target = calc_bc, args=(G, return_dict))
    process2 = multiprocessing.Process(target = calc_shortest_pthlen, args=(G, return_dict))
    process1.start()
    process2.start()

    G_feat = {}  # graph features dictionary
    G_feat['n nodes'] = float(G.number_of_nodes())  # number of nodes
    G_feat['n edges'] = float(G.number_of_edges())  # number of edges
    G_feat['density'] = nx.density(G)  # how close the network is to a 'complete graph' where each node is connected
    G_feat['diameter'] = nx.diameter(G)  # the farthest distance (e.g. number of edges) between two nodes in the graph
    G_feat['avg clustering'] = nx.average_clustering(G)  # the (averaged) fraction of possible triangles through a node.
    G_feat['max closeness centrality'] = np.max(list(nx.closeness_centrality(G).values()))  # max closeness centraility. high closeness --> short distance to all other nodes
    G_feat['max eigenvector centrality'] = np.max(list(nx.eigenvector_centrality(G, max_iter=10000).values()))  # eigenvector centraility. how 'important'/'influential' a node is
    G_feat['degree assortativity'] = nx.degree_pearson_correlation_coefficient(G)  # tendency of a node to be connected with other nodes of the same degree
    G_feat['max clique number'] = nx.graph_clique_number(G)  # largest clique (i.e. an induced subgraph that is complete) size 
    G_feat['n communities'] = len(nx.algorithms.community.modularity_max.greedy_modularity_communities(G))  # number of communities

    # finish multiple processes
    process2.join()
    process1.join()
    G_feat['avg path length'] = return_dict[2]
    G_feat['max betweenness centrality'] = return_dict[1]
    return G_feat


def similarity_measure(x1, x2):
    """ calculates the similarity between two feature values.
    similarity = 1 - the relative distance between features (x1 and x2) 

    Args:
        x1 ([float]): feature from graph 1 (must range between 0,1)
        x2 ([float]): feature from graph 2 (must range between 0,1)

    Returns:
        [float]: returns the relative similarity between 2 features
    """
  
    return 1-(np.abs(x1 - x2)/max(x1,x2))


def calc_similarity_score(G1_dict, G2_dict, feature_list):
    """calculates the similarity score of two graphs

    Args:
        G1_dict ([dict] or [Pandas datafrane]): graph 1 features dictionary or dataframe. must be able to use a key to access values
        G2_dict ([dict] or [Pandas datafrane]): graph 2 features dictionary or dataframe. must be able to use a key to access values
        features_list ([list]): list of graph features to compare. must be keys in graph features dictionary (above)

    Returns:
        [float]: similarity score (0,1) where 1 is an identical graph. 
    """
    s_list = []
    for feat in feature_list:
        f1 = G1_dict[feat]
        f2 = G2_dict[feat]
        s_tmp = similarity_measure(f1,f2)
        s_list.append(s_tmp)
    s = np.sum(s_list)/len(s_list)
    return s  


def process_graphs(graph_fnames):
    """ take a list of graph files, calculate their features, and return as a dataframe

    Args:
        graph_fnames ([list]): list of graph filenames to process

    Returns:
        [pandas dataframe]: dataframe containing graph features for each graph in filename list
    """
    tmp_feat_list = []
    for fname in graph_fnames:
        G_tmp = nx.read_gexf(fname)
        tmp_feat_list.append(calc_graph_features(G_tmp))
    tmp_df = pd.DataFrame(tmp_feat_list)
    tmp_df['name'] = [Path(i).stem for i in graph_fnames]
    assert(not tmp_df.empty)
    return tmp_df


def classify_graphs(class_file_list, sample_file_list, feature_list):
    """Classifies a similarity score from a list of Class and Sample graphs

    Args:
        class_file_list ([list]):  list of control/reference graph files (classes)
        sample_file_list ([list]): list of non-control/non-reference graph files (samples)
        feature_list ([list]): list of which features to use for similarity score. must be a valid key to the graph features dictionary/dataframe (above)

    Returns:
        [pandas dataframe]: each colummn is a class and each row is the similarity score of the sampled graph
    """
    class_feat_df = process_graphs(class_file_list)
    sample_feat_df = process_graphs(sample_file_list)

    class_list = []
    for index, row in sample_feat_df.iterrows():  # for each sample graph
        tmp_graph_feat = row
        class_dict = {}
        class_dict['name'] = tmp_graph_feat['name']
        for index2, row2 in class_feat_df.iterrows():  # for each class graph
            tmp_class_feat = row2
            s_tmp = calc_similarity_score(tmp_graph_feat, tmp_class_feat,feature_list)
            class_dict[tmp_class_feat['name']] = s_tmp
        class_list.append(class_dict)
    class_similarity_df = pd.DataFrame(class_list)
    return class_similarity_df


def process_similarity_df(class_similarity_df):
    """Generates y_true and y_pred based on the similarity score dataframe. 

    y_true is a list where each index is a class and each value is the class value. E.g., class 1 is y_true[1] = 1, class 2 is y_true[2]=2, etc.

    y_pred is a list where each index is a sample and each value is the maximum similarity score for that sample.

    Note: This assumes the correct classification is along the diagonal of the similarity matrix/dataframe. 

    ref: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.classification_report.html#sklearn.metrics.classification_report

    Args:
        class_similarity_df ([pandas dataframe]): each column is a class graph and each row is a sample graph. A_ij is the similarity score between graphs i and j. The exception is one column 'name' which contains the names of the sampled graphs for each row. 

    Returns:
        ([tuple of lists]): y_true, y_pred
    """
    num_df = class_similarity_df.drop(columns='name')
    y_true = [i for i in range(len(class_similarity_df.columns)-1)]
    y_pred = list(num_df.idxmax(axis=0))
    return y_true, y_pred
