from parsl import python_app


@python_app
def convert_density_file_to_centroids(
    protein_name,
    mrc_file_path,
    DBSCAN_epsilon,
    DBSCAN_min_samples,
    density_threshold,
    output_dir,
):
    """
    Convert a density file to centroids using DBSCAN clustering.

    Parameters
    ----------
    protein_name : str
        Name of the protein.
    mrc_file_path : str
        Path to the MRC density file.
    DBSCAN_epsilon : float
        The maximum distance between two samples for one to be considered as in the neighborhood of the other.
    DBSCAN_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.
    density_threshold : float
        The threshold value for density.
    output_dir : str
        Directory where the output will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - centroids (numpy.ndarray): The centroids of the clusters.
        - centroid_dir (str): The directory where the centroids are saved.
    """
    import os
    import numpy as np
    import csv

    from griptomo.core import density2graph as d2g

    fname = os.path.basename(mrc_file_path)
    density_name = fname[:-4]

    centroid_dir = os.path.join(
        output_dir,
        protein_name,
        density_name,
        f"density_{density_threshold}",
        f"min_samples_{DBSCAN_min_samples}",
        f"epsilon_{DBSCAN_epsilon}",
    )
    if os.path.exists(centroid_dir):
        print(f"{centroid_dir} already exists.")
    else:
        os.makedirs(centroid_dir)
        print(f"{centroid_dir} does not exist. Creating it.")

    centroid_filename = os.path.join(
        centroid_dir,
        "centroids.npy",
    )

    if os.path.exists(centroid_filename):
        print(f"{centroid_filename} already exists.")
        centroids = np.load(centroid_filename)
    else:
        print(f"{centroid_filename} does not exist. Generating coarse model.")
        ### load density file, normalize the data and threshold, and extract the x,y,z coordinates of the remaining pixels
        mrc = d2g.load_density_file(mrc_file_path)  # load density file
        D_norm = d2g.normalize_mrc_data(mrc.data)
        total_voxel_count = mrc.data.size
        D_thresholded = d2g.threshold_mrc_data(D_norm, density_threshold)
        xyz_data = d2g.generate_point_cloud_from_mrc_data(
            D_thresholded, density_threshold
        )
        point_cloud_node_count = len(xyz_data)

        ### perform clustering and get cluster centers
        model = d2g.cluster_data(
            xyz_data, DBSCAN_epsilon, DBSCAN_min_samples
        )  # cluster thresholded data using DBSCAN
        noisy_node_count = np.sum(model.labels_ == -1)
        non_noise_cluster_count = len(np.unique(model.labels_))
        if noisy_node_count > 0:
            non_noise_cluster_count -= 1
        centroids = d2g.get_cluster_centroids(
            xyz_data, model
        )  # coarse grain model by getting cluster centroids
        # Save centroids in output folder.
        np.save(centroid_filename, centroids)
        # Save centroids_run_info in output folder.
        run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
        with open(run_info_filename, "w", newline="") as csvfile:
            fieldnames = [
                "Total Voxel Count",
                "Point Cloud Node Count",
                "Non Noise Cluster Count",
                "Noisy Node Count",
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow(
                {
                    "Total Voxel Count": total_voxel_count,
                    "Point Cloud Node Count": point_cloud_node_count,
                    "Non Noise Cluster Count": non_noise_cluster_count,
                    "Noisy Node Count": noisy_node_count,
                }
            )
    return (centroids, centroid_dir)


@python_app
def extract_graph_features_from_centroids(
    centroids,
    cutoff,
    output_dir,
    force_overwrite=False,
):
    """
    Extracts graph features from centroid coordinates and saves them to a CSV file.

    Parameters
    ----------
    centroids : str or numpy.ndarray
        Path to a file containing centroid coordinates or a numpy array of centroid coordinates.
    cutoff : float
        Distance threshold for creating edges in the graph.
    output_dir : str
        Directory where the output CSV file will be saved.
    force_overwrite : bool, optional
        If True, overwrite the existing output file. Defaults to False.

    Returns
    -------
    str
        Path to the output CSV file containing the graph features.
    """
    import os

    feature_output_dir = os.path.join(output_dir, f"cutoff_{cutoff}")

    if os.path.exists(feature_output_dir):
        print(f"{feature_output_dir} already exists.")
    else:
        os.makedirs(feature_output_dir)
        print(f"{feature_output_dir} does not exist. Creating it.")

    feature_output_filename = os.path.join(feature_output_dir, "graph_features.csv")

    if os.path.exists(feature_output_filename) and not force_overwrite:
        print(f"{feature_output_filename} already exists and force_overwrite is False.")
    else:
        import numpy as np
        import csv

        from griptomo.core import ig_extract_features as g2c
        from griptomo.core.density2graph import create_igraph_from_point_cloud

        # If the centroids are not already loaded, load them.
        if isinstance(centroids, str):
            centroids = np.load(centroids)

        graph = create_igraph_from_point_cloud(centroids, cutoff)
        graph_features = g2c.igraph_calc_graph_features(graph)

        if os.path.exists(feature_output_filename):
            print(f"{feature_output_filename} already exists but will be overwritten.")
        else:
            print(f"{feature_output_filename} does not exist. Writing graph features.")
        with open(feature_output_filename, "w", newline="") as csvfile:
            fieldnames = list(graph_features.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(graph_features)

    return feature_output_filename


@python_app
def convert_density_file_to_centroids_timed(
    protein_name,
    mrc_file_path,
    DBSCAN_epsilon,
    DBSCAN_min_samples,
    density_threshold,
    output_dir,
):
    """
    Convert a density file to centroids using DBSCAN clustering.
    Additionally, save the timings of each step to a CSV file.

    Parameters
    ----------
    protein_name : str
        Name of the protein.
    mrc_file_path : str
        Path to the MRC density file.
    DBSCAN_epsilon : float
        The maximum distance between two samples for one to be considered as in the neighborhood of the other.
    DBSCAN_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.
    density_threshold : float
        The threshold value for density.
    output_dir : str
        Directory where the output will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - centroids (numpy.ndarray): The centroids of the clusters.
        - centroid_dir (str): The directory where the centroids are saved.
    """
    import time

    times = {}

    start_time = time.perf_counter()
    import os
    import numpy as np
    import csv

    from griptomo.core import density2graph as d2g

    times["module import"] = time.perf_counter() - start_time

    fname = os.path.basename(mrc_file_path)
    density_name = fname[:-4]

    centroid_dir = os.path.join(
        output_dir,
        protein_name,
        density_name,
        f"density_{density_threshold}",
        f"min_samples_{DBSCAN_min_samples}",
        f"epsilon_{DBSCAN_epsilon}",
    )
    if os.path.exists(centroid_dir):
        print(f"{centroid_dir} already exists.")
    else:
        os.makedirs(centroid_dir)
        print(f"{centroid_dir} does not exist. Creating it.")

    centroid_filename = os.path.join(
        centroid_dir,
        "centroids.npy",
    )

    if os.path.exists(centroid_filename):
        print(f"{centroid_filename} already exists.")
        start_time = time.perf_counter()
        centroids = np.load(centroid_filename)
        times["load centroids"] = time.perf_counter() - start_time
    else:
        print(f"{centroid_filename} does not exist. Generating coarse model.")
        ### load density file, normalize the data and threshold, and extract the x,y,z coordinates of the remaining pixels
        start_time = time.perf_counter()
        mrc = d2g.load_density_file(mrc_file_path)  # load density file
        times["load density file"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        D_norm = d2g.normalize_mrc_data(mrc.data)
        times["normalize data"] = time.perf_counter() - start_time
        total_voxel_count = mrc.data.size

        start_time = time.perf_counter()
        D_thresholded = d2g.threshold_mrc_data(D_norm, density_threshold)
        times["threshold data"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        xyz_data = d2g.generate_point_cloud_from_mrc_data(
            D_thresholded, density_threshold
        )
        times["generate point cloud"] = time.perf_counter() - start_time

        point_cloud_node_count = len(xyz_data)

        start_time = time.perf_counter()
        ### perform clustering and get cluster centers
        model = d2g.cluster_data(
            xyz_data, DBSCAN_epsilon, DBSCAN_min_samples
        )  # cluster thresholded data using DBSCAN
        times["cluster data"] = time.perf_counter() - start_time

        noisy_node_count = np.sum(model.labels_ == -1)
        non_noise_cluster_count = len(np.unique(model.labels_))
        if noisy_node_count > 0:
            non_noise_cluster_count -= 1

        start_time = time.perf_counter()
        centroids = d2g.get_cluster_centroids(
            xyz_data, model
        )  # coarse grain model by getting cluster centroids
        times["get cluster centroids"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        np.save(centroid_filename, centroids)
        times["save centroids"] = time.perf_counter() - start_time

        run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
        with open(run_info_filename, "w", newline="") as csvfile:
            fieldnames = [
                "Total Voxel Count",
                "Point Cloud Node Count",
                "Non Noise Cluster Count",
                "Noisy Node Count",
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow(
                {
                    "Total Voxel Count": total_voxel_count,
                    "Point Cloud Node Count": point_cloud_node_count,
                    "Non Noise Cluster Count": non_noise_cluster_count,
                    "Noisy Node Count": noisy_node_count,
                }
            )
            # Print timing information to a file
            times_filename = centroid_filename.replace("centroids.npy", "times.csv")
            with open(times_filename, "w", newline="") as csvfile:
                fieldnames = list(times.keys())
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(times)
    return (centroids, centroid_dir)


@python_app
def extract_graph_features_from_centroids_timed(
    centroids,
    cutoff,
    output_dir,
    force_overwrite=False,
    skip_clique_num=False,
):
    """
    Extracts graph features from centroid coordinates and saves them to a CSV file.
    Additionally, save the times of each step to a separate CSV file.

    Parameters
    ----------
    centroids : str or numpy.ndarray
        Path to a file containing centroid coordinates or a numpy array of centroid coordinates.
    cutoff : float
        Distance threshold for creating edges in the graph.
    output_dir : str
        Directory where the output CSV file will be saved.
    force_overwrite : bool, optional
        If True, overwrite the existing output file. Defaults to False.
    skip_clique_num : bool, optional
        If True, skip the calculation of the number of maximal cliques. Defaults to False.

    Returns
    -------
    str
        Path to the output CSV file containing the graph features.
    """
    import os

    feature_output_dir = os.path.join(output_dir, f"cutoff_{cutoff}")

    if os.path.exists(feature_output_dir):
        print(f"{feature_output_dir} already exists.")
    else:
        os.makedirs(feature_output_dir)
        print(f"{feature_output_dir} does not exist. Creating it.")

    feature_output_filename = os.path.join(feature_output_dir, "graph_features.csv")

    if os.path.exists(feature_output_filename) and not force_overwrite:
        print(f"{feature_output_filename} already exists and force_overwrite is False.")
    else:
        import time

        times = {}

        start_time = time.perf_counter()
        import numpy as np
        import csv

        from griptomo.core import ig_extract_features as g2c
        from griptomo.core.density2graph import create_igraph_from_point_cloud

        times["module import"] = time.perf_counter() - start_time

        # If the centroids are not already loaded, load them.
        if isinstance(centroids, str):
            centroids = np.load(centroids)

        graph = create_igraph_from_point_cloud(centroids, cutoff)
        times["convert to igraph"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        graph_features, graph_times = g2c.igraph_calc_graph_features_timed(
            graph, skip_clique_num=skip_clique_num
        )
        times["extract features"] = time.perf_counter() - start_time

        if os.path.exists(feature_output_filename):
            print(f"{feature_output_filename} already exists but will be overwritten.")
        else:
            print(f"{feature_output_filename} does not exist. Writing graph features.")
        with open(feature_output_filename, "w", newline="") as csvfile:
            fieldnames = list(graph_features.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(graph_features)
        graph_times_filename = feature_output_filename.replace(
            "graph_features.csv", "graph_times.csv"
        )
        with open(graph_times_filename, "w", newline="") as csvfile:
            fieldnames = list(graph_times.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(graph_times)
        times_filename = feature_output_filename.replace(
            "graph_features.csv", "times.csv"
        )
        with open(times_filename, "w", newline="") as csvfile:
            fieldnames = list(times.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(times)

    return feature_output_filename


@python_app
def convert_density_file_to_point_cloud_hdbscan_timed(
    protein_name,
    mrc_file_path,
    HDBSCAN_min_cluster_size,
    HDBSCAN_min_samples,
    coarsening_dim,
    density_threshold,
    output_dir,
    averaged=False,
    density_coarsening=False,
    precoarsened=False,
    standardize=False,
    threshold_ratio=False,
):
    """
    Convert a density file to centroids using HDBSCAN clustering.
    Additionally, save the timings of each step to a CSV file.

    Parameters
    ----------
    protein_name : str
        Name of the protein.
    mrc_file_path : str
        Path to the MRC density file.
    HDBSCAN_min_cluster_size : float
        The minimum number of points in a cluster.
    HDBSCAN_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.
    coarsening_dim : int
    The dimensionality passed in to griptomo.core.density2graph.voxel_coarsening's voxel_dim argument to be used for performing coarsening on the thresholded point cloud.
    density_threshold : float
        The threshold value for density.
    output_dir : str
        Directory where the output will be saved.
    averaged : bool, optional
        If True, average the points in the coarsening box instead of returning the center of the box. Defaults to False.
    density_coarsening : bool, optional
        If True, use the density values of each voxel in the box to provide weight when averaging during coarsening. Defaults to False.
    precoarsened : bool, optional
        If True, pre-coarsen the data before extracting the cluster centroids. Defaults to False.
    standardize : bool, optional
        If True, standardize the data before thresholding. Defaults to False.
    threshold_ratio : bool, optional
        If True, use the threshold ratio to identify the threshold. Defaults to False.

    Returns
    -------
    tuple
        A tuple containing:
        - centroids (numpy.ndarray): The centroids of the clusters.
        - centroid_dir (str): The directory where the centroids are saved.
    """
    import time

    times = {}

    start_time = time.perf_counter()
    import os
    import numpy as np
    import csv

    from griptomo.core import density2graph as d2g

    times["module import"] = time.perf_counter() - start_time

    fname = os.path.basename(mrc_file_path)
    density_name = fname[:-4]

    centroid_dir = os.path.join(
        output_dir,
        protein_name,
        density_name,
        f"density_{density_threshold}",
        f"min_samples_{HDBSCAN_min_samples}",
        f"min_cluster_size_{HDBSCAN_min_cluster_size}",
        f"coarsening_dim_{coarsening_dim}",
    )
    if os.path.exists(centroid_dir):
        print(f"{centroid_dir} already exists.")
    else:
        print(f"{centroid_dir} does not exist. Creating it.")
        os.makedirs(centroid_dir)

    centroid_filename = os.path.join(
        centroid_dir,
        "centroids.npy",
    )

    if os.path.exists(centroid_filename):
        print(f"{centroid_filename} already exists.")
        start_time = time.perf_counter()
        centroids = np.load(centroid_filename)
        times["load centroids"] = time.perf_counter() - start_time
    else:
        print(f"{centroid_filename} does not exist. Generating coarse model.")
        ### load density file, normalize the data and threshold, and extract the x,y,z coordinates of the remaining pixels
        start_time = time.perf_counter()
        mrc = d2g.load_density_file(mrc_file_path)  # load density file
        times["load density file"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        if standardize:
            D_standardized = d2g.standardize_mrc_data(mrc.data)
            times["standardize data"] = time.perf_counter() - start_time
            start_time = time.perf_counter()
            D_norm = d2g.normalize_mrc_data(D_standardized)
            times["normalize data"] = time.perf_counter() - start_time
        else:
            D_norm = d2g.normalize_mrc_data(mrc.data)
            times["normalize data"] = time.perf_counter() - start_time
        total_voxel_count = mrc.data.size

        start_time = time.perf_counter()
        if threshold_ratio:
            threshold = d2g.identify_threshold_ratio(D_norm, density_threshold)
        else:
            threshold = density_threshold
        start_time = time.perf_counter()
        xyz_data = d2g.generate_point_cloud_from_mrc_data(D_norm, threshold)
        times["generate point cloud"] = time.perf_counter() - start_time

        if density_coarsening:
            xyz_densities = D_norm[xyz_data[:, 0], xyz_data[:, 1], xyz_data[:, 2]] + 1
            xyz_densities = xyz_densities / 2
            fine_mean_density = np.mean(xyz_densities)
            fine_median_density = np.median(xyz_densities)
            fine_max_density = np.max(xyz_densities)
            fine_min_density = np.min(xyz_densities)
            fine_std_density = np.std(xyz_densities)
            fine_sum_density = np.sum(xyz_densities)
            xyz_densities_filename = os.path.join(
                centroid_dir,
                "xyz_densities.npy",
            )
            np.save(xyz_densities_filename, xyz_densities)
        else:
            xyz_densities = None

        point_cloud_node_count = len(xyz_data)

        if precoarsened:
            start_time = time.perf_counter()
            coarsened = d2g.voxel_coarsening(
                coarsening_dim,
                xyz_data,
                np.zeros(xyz_data.shape[0], dtype=np.int64),
                density=xyz_densities,
                averaged=averaged,
            )
            if xyz_densities is not None:
                xyz_densities = coarsened[1]
                xyz_data = coarsened[0]
            else:
                xyz_data = coarsened
            times["precoarsening"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        ### perform clustering and get cluster centers
        model = d2g.cluster_data_hdbscan(
            xyz_data, HDBSCAN_min_cluster_size, HDBSCAN_min_samples
        )  # cluster thresholded data using DBSCAN
        times["cluster data"] = time.perf_counter() - start_time

        noisy_node_count = np.sum(model.labels_ == -1)
        non_noise_cluster_count = len(np.unique(model.labels_))
        if noisy_node_count > 0:
            non_noise_cluster_count -= 1
        if non_noise_cluster_count == 0:
            print("No clusters found for {centroid_dir}")
            run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
            with open(run_info_filename, "w", newline="") as csvfile:
                fieldnames = ["Total Voxel Count", "Point Cloud Node Count"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                    }
                )
            # Print timing information to a file
            times_filename = centroid_filename.replace("centroids.npy", "times.csv")
            with open(times_filename, "w", newline="") as csvfile:
                fieldnames = list(times.keys())
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(times)
            raise ValueError("No clusters found.")
        largest_cluster = np.argmax(np.bincount(model.labels_[model.labels_ != -1]))
        largest_cluster_size = np.sum(model.labels_ == largest_cluster)

        xyz_filename = os.path.join(
            centroid_dir,
            "xyz_data.npy",
        )
        labels_filename = os.path.join(
            centroid_dir,
            "labels.npy",
        )
        np.save(xyz_filename, xyz_data)
        np.save(labels_filename, model.labels_)

        if not precoarsened:
            start_time = time.perf_counter()
            coarsened = d2g.voxel_coarsening(
                coarsening_dim,
                xyz_data,
                model.labels_,
                density=xyz_densities,
                averaged=averaged,
            )  # coarse grain model by getting cluster centroids
            if density_coarsening:
                centroids, centroid_densities = coarsened
            else:
                centroids = coarsened
            coarsened_size = centroids.shape[0]
            times["voxel coarsening"] = time.perf_counter() - start_time
        else:
            centroids = xyz_data[model.labels_ == largest_cluster]
            centroid_densities = (
                xyz_densities[model.labels_ == largest_cluster]
                if xyz_densities is not None
                else None
            )
            coarsened_size = centroids.shape[0]
        if density_coarsening:
            coarse_mean_density = np.mean(centroid_densities)
            coarse_median_density = np.median(centroid_densities)
            coarse_max_density = np.max(centroid_densities)
            coarse_min_density = np.min(centroid_densities)
            coarse_std_density = np.std(centroid_densities)
            coarse_sum_density = np.sum(centroid_densities)
        else:
            centroids = coarsened

        coarsened_size = centroids.shape[0]
        times["voxel coarsening"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        if density_coarsening:
            centroid_densities_filename = os.path.join(
                centroid_dir,
                "centroid_densities.npy",
            )
            np.save(centroid_densities_filename, centroid_densities)
        np.save(centroid_filename, centroids)
        times["save centroids"] = time.perf_counter() - start_time

        run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
        with open(run_info_filename, "w", newline="") as csvfile:
            fieldnames = [
                "Total Voxel Count",
                "Point Cloud Node Count",
                "Non Noise Cluster Count",
                "Noisy Node Count",
                "Largest Cluster Size",
                "Coarsened Size",
            ]
            if density_coarsening:
                fieldnames += [
                    "Fine Mean Density",
                    "Fine Median Density",
                    "Fine Max Density",
                    "Fine Min Density",
                    "Fine Std Density",
                    "Fine Sum Density",
                    "Coarse Mean Density",
                    "Coarse Median Density",
                    "Coarse Max Density",
                    "Coarse Min Density",
                    "Coarse Std Density",
                    "Coarse Sum Density",
                ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            if density_coarsening:
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                        "Non Noise Cluster Count": non_noise_cluster_count,
                        "Noisy Node Count": noisy_node_count,
                        "Largest Cluster Size": largest_cluster_size,
                        "Coarsened Size": coarsened_size,
                        "Fine Mean Density": fine_mean_density,
                        "Fine Median Density": fine_median_density,
                        "Fine Max Density": fine_max_density,
                        "Fine Min Density": fine_min_density,
                        "Fine Std Density": fine_std_density,
                        "Fine Sum Density": fine_sum_density,
                        "Coarse Mean Density": coarse_mean_density,
                        "Coarse Median Density": coarse_median_density,
                        "Coarse Max Density": coarse_max_density,
                        "Coarse Min Density": coarse_min_density,
                        "Coarse Std Density": coarse_std_density,
                        "Coarse Sum Density": coarse_sum_density,
                    }
                )
            else:
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                        "Non Noise Cluster Count": non_noise_cluster_count,
                        "Noisy Node Count": noisy_node_count,
                        "Largest Cluster Size": largest_cluster_size,
                        "Coarsened Size": coarsened_size,
                    }
                )
        # Print timing information to a file
        times_filename = centroid_filename.replace("centroids.npy", "times.csv")
        with open(times_filename, "w", newline="") as csvfile:
            fieldnames = list(times.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(times)
    return (centroids, centroid_dir)


@python_app
def convert_density_file_to_point_cloud_optics_timed(
    protein_name,
    mrc_file_path,
    OPTICS_min_samples,
    coarsening_dim,
    density_threshold,
    output_dir,
    OPTICS_min_cluster_size=None,
    OPTICS_max_epsilon=None,
    averaged=False,
    density_coarsening=False,
):
    """
    Convert a density file to centroids using OPTICS clustering.
    Additionally, save the timings of each step to a CSV file.

    Parameters
    ----------
    protein_name : str
        Name of the protein.
    mrc_file_path : str
        Path to the MRC density file.
    OPTICS_min_samples : int
        The number of samples in a neighborhood for a point to be considered as a core point.
    coarsening_dim : int
    The dimensionality passed in to griptomo.core.density2graph.voxel_coarsening's voxel_dim argument to be used for performing coarsening on the thresholded point cloud.
    density_threshold : float
        The threshold value for density.
    output_dir : str
        Directory where the output will be saved.
    OPTICS_min_cluster_size : float, optional
        The minimum number of points in a cluster. Defaults to None.
    OPTICS_max_epsilon : float, optional
        The maximum distance between points considered by OPTICS. Defaults to None, which uses np.inf. A lower value may reduce computation time.
    averaged : bool, optional
        If True, average the points in the coarsening box instead of returning the center of the box. Defaults to False.
    density_coarsening : bool, optional
        If True, use the density values of each voxel in the box to provide weight when averaging during coarsening. Defaults to False.

    Returns
    -------
    tuple
        A tuple containing:
        - centroids (numpy.ndarray): The centroids of the clusters.
        - centroid_dir (str): The directory where the centroids are saved.
    """
    import time

    times = {}

    start_time = time.perf_counter()
    import os
    import numpy as np
    import csv

    from griptomo.core import density2graph as d2g

    times["module import"] = time.perf_counter() - start_time

    OPTICS_max_epsilon = (
        OPTICS_max_epsilon if OPTICS_max_epsilon is not None else np.inf
    )

    fname = os.path.basename(mrc_file_path)
    density_name = fname[:-4]

    centroid_dir = os.path.join(
        output_dir,
        protein_name,
        density_name,
        f"density_{density_threshold}",
        f"min_samples_{OPTICS_min_samples}",
        f"max_eps_{OPTICS_max_epsilon}",
        f"coarsening_dim_{coarsening_dim}",
    )
    if os.path.exists(centroid_dir):
        print(f"{centroid_dir} already exists.")
    else:
        print(f"{centroid_dir} does not exist. Creating it.")
        os.makedirs(centroid_dir)

    centroid_filename = os.path.join(
        centroid_dir,
        "centroids.npy",
    )

    if os.path.exists(centroid_filename):
        print(f"{centroid_filename} already exists.")
        start_time = time.perf_counter()
        centroids = np.load(centroid_filename)
        times["load centroids"] = time.perf_counter() - start_time
    else:
        print(f"{centroid_filename} does not exist. Generating coarse model.")
        ### load density file, normalize the data and threshold, and extract the x,y,z coordinates of the remaining pixels
        start_time = time.perf_counter()
        mrc = d2g.load_density_file(mrc_file_path)  # load density file
        times["load density file"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        D_norm = d2g.normalize_mrc_data(mrc.data)
        times["normalize data"] = time.perf_counter() - start_time
        total_voxel_count = mrc.data.size

        start_time = time.perf_counter()
        D_thresholded = d2g.threshold_mrc_data(D_norm, density_threshold)
        times["threshold data"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        xyz_data = d2g.generate_point_cloud_from_mrc_data(
            D_thresholded, density_threshold
        )
        times["generate point cloud"] = time.perf_counter() - start_time

        point_cloud_node_count = len(xyz_data)

        start_time = time.perf_counter()
        ### perform clustering and get cluster centers
        model = d2g.cluster_data_optics(
            xyz_data, OPTICS_min_samples, OPTICS_min_cluster_size, OPTICS_max_epsilon
        )  # cluster thresholded data using DBSCAN
        times["cluster data"] = time.perf_counter() - start_time

        noisy_node_count = np.sum(model.labels_ == -1)
        non_noise_cluster_count = len(np.unique(model.labels_))
        if noisy_node_count > 0:
            non_noise_cluster_count -= 1
        if non_noise_cluster_count == 0:
            print("No clusters found for {centroid_dir}")
            run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
            with open(run_info_filename, "w", newline="") as csvfile:
                fieldnames = ["Total Voxel Count", "Point Cloud Node Count"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                    }
                )
            # Print timing information to a file
            times_filename = centroid_filename.replace("centroids.npy", "times.csv")
            with open(times_filename, "w", newline="") as csvfile:
                fieldnames = list(times.keys())
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(times)
            raise ValueError("No clusters found.")
        largest_cluster = np.argmax(np.bincount(model.labels_[model.labels_ != -1]))
        largest_cluster_size = np.sum(model.labels_ == largest_cluster)

        if density_coarsening:
            xyz_densities = D_norm[xyz_data[:, 0], xyz_data[:, 1], xyz_data[:, 2]] + 1
            xyz_densities = xyz_densities / 2
            fine_mean_density = np.mean(xyz_densities)
            fine_median_density = np.median(xyz_densities)
            fine_max_density = np.max(xyz_densities)
            fine_min_density = np.min(xyz_densities)
            fine_std_density = np.std(xyz_densities)
            fine_sum_density = np.sum(xyz_densities)
            xyz_densities_filename = os.path.join(
                centroid_dir,
                "xyz_densities.npy",
            )
            np.save(xyz_densities_filename, xyz_densities)
        else:
            xyz_densities = None

        xyz_filename = os.path.join(
            centroid_dir,
            "xyz_data.npy",
        )
        labels_filename = os.path.join(
            centroid_dir,
            "labels.npy",
        )
        np.save(xyz_filename, xyz_data)
        np.save(labels_filename, model.labels_)

        start_time = time.perf_counter()
        coarsened = d2g.voxel_coarsening(
            coarsening_dim,
            xyz_data,
            model.labels_,
            density=xyz_densities,
            averaged=averaged,
        )  # coarse grain model by getting cluster centroids
        if density_coarsening:
            centroids, centroid_densities = coarsened
            coarse_mean_density = np.mean(centroid_densities)
            coarse_median_density = np.median(centroid_densities)
            coarse_max_density = np.max(centroid_densities)
            coarse_min_density = np.min(centroid_densities)
            coarse_std_density = np.std(centroid_densities)
            coarse_sum_density = np.sum(centroid_densities)
        else:
            centroids = coarsened

        coarsened_size = centroids.shape[0]
        times["voxel coarsening"] = time.perf_counter() - start_time

        start_time = time.perf_counter()
        if density_coarsening:
            centroid_densities_filename = os.path.join(
                centroid_dir,
                "centroid_densities.npy",
            )
            np.save(centroid_densities_filename, centroid_densities)
        np.save(centroid_filename, centroids)
        times["save centroids"] = time.perf_counter() - start_time

        run_info_filename = centroid_filename.replace(".npy", "_run_info.csv")
        with open(run_info_filename, "w", newline="") as csvfile:
            fieldnames = [
                "Total Voxel Count",
                "Point Cloud Node Count",
                "Non Noise Cluster Count",
                "Noisy Node Count",
                "Largest Cluster Size",
                "Coarsened Size",
            ]
            if density_coarsening:
                fieldnames += [
                    "Fine Mean Density",
                    "Fine Median Density",
                    "Fine Max Density",
                    "Fine Min Density",
                    "Fine Std Density",
                    "Fine Sum Density",
                    "Coarse Mean Density",
                    "Coarse Median Density",
                    "Coarse Max Density",
                    "Coarse Min Density",
                    "Coarse Std Density",
                    "Coarse Sum Density",
                ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            if density_coarsening:
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                        "Non Noise Cluster Count": non_noise_cluster_count,
                        "Noisy Node Count": noisy_node_count,
                        "Largest Cluster Size": largest_cluster_size,
                        "Coarsened Size": coarsened_size,
                        "Fine Mean Density": fine_mean_density,
                        "Fine Median Density": fine_median_density,
                        "Fine Max Density": fine_max_density,
                        "Fine Min Density": fine_min_density,
                        "Fine Std Density": fine_std_density,
                        "Fine Sum Density": fine_sum_density,
                        "Coarse Mean Density": coarse_mean_density,
                        "Coarse Median Density": coarse_median_density,
                        "Coarse Max Density": coarse_max_density,
                        "Coarse Min Density": coarse_min_density,
                        "Coarse Std Density": coarse_std_density,
                        "Coarse Sum Density": coarse_sum_density,
                    }
                )
            else:
                writer.writerow(
                    {
                        "Total Voxel Count": total_voxel_count,
                        "Point Cloud Node Count": point_cloud_node_count,
                        "Non Noise Cluster Count": non_noise_cluster_count,
                        "Noisy Node Count": noisy_node_count,
                        "Largest Cluster Size": largest_cluster_size,
                        "Coarsened Size": coarsened_size,
                    }
                )
        # Print timing information to a file
        times_filename = centroid_filename.replace("centroids.npy", "times.csv")
        with open(times_filename, "w", newline="") as csvfile:
            fieldnames = list(times.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(times)
    return (centroids, centroid_dir)
