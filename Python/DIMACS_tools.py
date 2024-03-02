import json
import time
from sklearn.cluster import AgglomerativeClustering
from scipy.sparse import coo_matrix, csgraph
from scipy.sparse.linalg import eigsh
from scipy.sparse import coo_matrix
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import os
import networkx as nx
import numpy as np

def pareto_width(high_cost, low_cost):
    return (high_cost - low_cost) / low_cost

def sparse_adjacency_matrix(edges_list):
    rows = edges_list[:, 0].astype(int) - 1
    columns = edges_list[:, 1].astype(int) - 1
    return coo_matrix((edges_list[:, 2].astype(int), (rows, columns))).tocsr()

def adjacency_matrix(vertices_count, edges_list):
    adj_matrix = np.zeros((vertices_count, vertices_count))
    for i in range(edges_list.shape[0]):
        start_node = int(edges_list[i, 0])
        end_node = int(edges_list[i, 1])
        cost = int(edges_list[i, 2])
        adj_matrix[start_node, end_node] = cost

    return adj_matrix

def generate_multiple_correlated_graph(distance_filename, time_filename, coords_filename,
                                       clusters_count, new_distance_gr_filename, new_time_gr_filename, new_coords_filename,
                                       clusters_metafile, min_lon=None, max_lon=None, min_lat=None, max_lat=None):
    # Reading the graphs
    c1_graph, vertices_count = load_graph(distance_filename)
    c2_graph, _ = load_graph(time_filename)
    edges_count = c1_graph.shape[0]
    coordinates = read_coords_file(coords_filename)

    # Computing the original correlation between the two cost functions
    original_corr = pearson_correlation(c1_graph[:, 2], c2_graph[:, 2])

    # Filtering only the desired area of interest
    if min_lon is not None:
        inside_nodes = np.array(np.where(
            (coordinates[:, 0] >= min_lon) & (coordinates[:, 0] <= max_lon) &
            (coordinates[:, 1] >= min_lat) & (coordinates[:, 1] <= max_lat)))

        inside_nodes = inside_nodes.reshape(inside_nodes.shape[1])

        coordinates = coordinates[inside_nodes, :]

        inside_nodes += 1

        valid_start_nodes = np.where(np.in1d(c1_graph[:, 0], inside_nodes))[0]
        valid_end_nodes = np.where(np.in1d(c1_graph[:, 1], inside_nodes))[0]
        intersection = np.intersect1d(valid_start_nodes, valid_end_nodes)
        c1_graph = c1_graph[intersection, :]
        c2_graph = c2_graph[intersection, :]

        vertices_count = coordinates.shape[0]
        edges_count = c1_graph.shape[0]
        vertices_set = set(coordinates[:, 2])

    # Computing the geographic east and north extent
    east_width = max(coordinates[:, 0]) - min(coordinates[:, 0])
    north_width = max(coordinates[:, 1]) - min(coordinates[:, 1])
    sw_corner_east = min(coordinates[:, 0])
    sw_corner_north = min(coordinates[:, 1])

    clusters_center_x = sw_corner_east + np.random.uniform(0, 1, clusters_count) * east_width
    clusters_center_y = sw_corner_north + np.random.uniform(0, 1, clusters_count) * north_width

    clusters_radius = np.random.uniform(east_width / 30, east_width / 20, clusters_count)

    updated_nodes = np.array([]).reshape(-1, 2)
    clusters_correlation = 0
    cluster_id = 0
    for cluster_id in range(clusters_count):
        print(f'Cluster {cluster_id + 1} / {clusters_count}...')
        # Retrieving the list of nodes indices that reside inside the bounding circle
        nodes = points_inside_circle(coordinates[:, 0:2], [clusters_center_x[cluster_id], clusters_center_y[cluster_id]],
                                     clusters_radius[cluster_id])

        # Keeping track of nodes to be updated
        nodes_with_corr = np.hstack((nodes.reshape(-1, 1), np.ones((len(nodes), 1)) * clusters_correlation))
        updated_nodes = np.concatenate((updated_nodes, nodes_with_corr))

        # Updating the weight of all the edges that leave the desired nodes
        edges_to_be_updated = np.array([])
        for node in nodes:
            edges_to_be_updated = np.concatenate((edges_to_be_updated, np.argwhere(c1_graph[:, 0] == node).flatten()))

        if len(edges_to_be_updated) > 0:
            distance_cost, time_cost = generate_correlated_vectors(len(edges_to_be_updated), target_correlation=clusters_correlation)
            c1_graph[edges_to_be_updated.astype(int), 2] = distance_cost
            c2_graph[edges_to_be_updated.astype(int), 2] = time_cost

        vertices_set = vertices_set - set(coordinates[nodes, 2])

    # Exporting new files
    export_gr_file(c1_graph, vertices_count, edges_count, new_distance_gr_filename,
                   'Distance graph for multiple correlated clusters')
    export_gr_file(c2_graph, vertices_count, edges_count, new_time_gr_filename,
                   'Time graph for multiple correlated clusters')
    export_coords_file(coordinates, vertices_count, new_coords_filename,
                       'Filtered NY coords file')

    clusters_mapping = np.zeros((len(vertices_set), 2))
    clusters_mapping[:, 0] = cluster_id
    clusters_mapping[:, 1] = np.array(list(vertices_set)).astype(int)

    export_clusters_file(clusters_mapping, clusters_metafile)

def load_graph(filename):
    edge_ind = 0
    with open(filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
                edges_count = int(edges_count_str)
                graph = np.zeros((edges_count, 3))
                continue

            # Inputting each edge data
            if line[0] == 'a':
                _, start_node, end_node, cost = line.rstrip('\n').split(' ')
                graph[edge_ind, 0] = int(start_node)
                graph[edge_ind, 1] = int(end_node)
                graph[edge_ind, 2] = cost
                edge_ind += 1

    return graph, vertices_count

def points_inside_circle(points, center, radius):
    """
    Identify points inside a circle.

    Parameters:
    - points: Nx2 array of (x, y) coordinates of points
    - center: (x, y) coordinates of the circle's center
    - radius: Radius of the circle

    Returns:
    - List of indices of points inside the circle
    """
    distances = np.linalg.norm(points - center, axis=1)
    indices_inside = np.where(distances <= radius)[0]
    return indices_inside

def create_sparse_adj_matrix(edge_list, num_nodes):
    # Create a COO sparse matrix from edge_list
    row_indices = [edge[0] for edge in edge_list]
    col_indices = [edge[1] for edge in edge_list]
    weights = [edge[2] for edge in edge_list]

    sparse_adj_matrix = coo_matrix((weights, (row_indices, col_indices)), shape=(num_nodes, num_nodes))
    return sparse_adj_matrix

def spectral_clustering(edge_list, num_clusters):
    # Extract nodes and create a mapping from node indices to consecutive integers
    nodes = list(set(node for edge in edge_list for node in edge))
    node_indices = {node: i for i, node in enumerate(nodes)}

    # Create a sparse adjacency matrix
    sparse_adj_matrix = create_sparse_adj_matrix(edge_list, len(nodes))

    # Compute the Laplacian matrix
    laplacian_matrix = csgraph.laplacian(sparse_adj_matrix, normed=False)

    # Perform a partial eigendecomposition using sparse eigensolver
    num_eigenvectors = min(num_clusters + 1, len(nodes) - 1)  # Choose the number of eigenvectors to compute
    eigenvalues, eigenvectors = eigsh(laplacian_matrix, k=num_eigenvectors, which='SM')  # 'SM' for smallest magnitude

    print('After eigsh, invoking SpectralClustering...')
    # Use the eigenvectors for spectral clustering
    spectral_model = SpectralClustering(n_clusters=num_clusters, affinity='nearest_neighbors', eigen_solver='arpack')
    cluster_labels = spectral_model.fit_predict(eigenvectors[:, 1:num_clusters+1])  # Exclude the first eigenvector (constant)

    return nodes, cluster_labels

def graph_agglomerative_clustering(distance_filename, time_filename, eps):
    c1_graph, vertices_count = load_graph(distance_filename)
    c2_graph, _ = load_graph(time_filename)

    # c1_adj_mat = sparse_adjacency_matrix(c1_graph)
    # c2_adj_mat = sparse_adjacency_matrix(c2_graph)

    # visited_vertices = np.zeros(vertices_count)
    # vertices_to_test = np.arange(vertices_count)

    # current_cluster_id = 0

    # while len(vertices_to_test) > 0:
    #     vertex = vertices_to_test[0]
    #     vertices_to_test = np.delete(vertices_to_test, 0)

    print('Building graph for networkx...')
    c1_graph = c1_graph[0:10000, :]
    G = nx.DiGraph()
    for i in range(c1_graph.shape[0]):
        start_node = int(c1_graph[i, 0])
        end_node = int(c1_graph[i, 1])
        cost = int(c1_graph[i, 2])
        G.add_edge(start_node, end_node, weight=cost)

    print('Running Floyd-Warshall...')
    t1 = time.time()
    result_matrix = nx.floyd_warshall(G)
    t2 = time.time()
    print(f'Floyd-Warshall took {t2-t1} [sec]')


def read_coords_file(coords_filename):
    # Reading the coords file
    with open(coords_filename, "r") as f:
        f.readline()
        row_id = 0
        for line in f:
            if line[0] == 'p':
                _, _, _, _, vertices_count = line.rstrip('\n').split(' ')
                coordinates = np.zeros((int(vertices_count), 3))
            if line[0] == 'v':
                _, node, x, y = line.rstrip('\n').split(' ')
                coordinates[row_id, :] = [int(x), int(y), int(node)]
                row_id += 1

    return coordinates

def plot_graph(adjacency_filename, coords_filename, clusters_metafile=None):
    # Reading the adjacency graph structure
    edge_ind = 0
    with open(adjacency_filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
                edges_count = int(edges_count_str)
                graph = np.zeros((edges_count, 2))
                continue

            # Inputting each edge data
            if line[0] == 'a':
                _, graph[edge_ind, 0], graph[edge_ind, 1], \
                    _ = line.rstrip('\n').split(' ')
                edge_ind += 1

    # Reading the coords file
    coordinates = read_coords_file(coords_filename)

    # Reading the clusters metafile
    correlated_nodes = np.array([])
    if clusters_metafile is not None:
        with open(clusters_metafile, "r") as f:
            f.readline()
            for line in f:
                cluster_id, node = line.rstrip('\n').split(' ')
                data = np.array([int(cluster_id), int(node)])
                if correlated_nodes.size == 0:
                    correlated_nodes = data
                else:
                    correlated_nodes = np.vstack((correlated_nodes, data))

        # Retrieving the coordinates row corresponding to the clusters' nodes ids
        ind = np.where(np.isin(coordinates[:, 2], correlated_nodes[:, 1]))
        ind = np.array(ind[0])

    marker_size = 1
    plt.scatter(coordinates[:, 0], coordinates[:, 1], color='red', marker='.', s=marker_size, label='0 Correlation')
    if clusters_metafile is not None:
        plt.scatter(coordinates[ind, 0], coordinates[ind, 1],
                    color='blue', marker='.', s=marker_size, label=f'0.96 Correlation')

    # Plotting the boundary nodes
    boundary_nodes_files = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\boundary_nodes.txt'
    bnd_nodes_id = np.array([])
    with open(boundary_nodes_files, "r") as f:
        f.readline()
        for line in f:
            node = line.rstrip('\n').split(' ')
            if bnd_nodes_id.size == 0:
                bnd_nodes_id = np.array(int(node[0]))
            else:
                bnd_nodes_id = np.vstack((bnd_nodes_id, int(node[0])))

    # Retrieving the coordinates row corresponding to the clusters' nodes ids
    bnd_ind = np.where(np.isin(coordinates[:, 2], bnd_nodes_id))
    bnd_ind = np.array(bnd_ind[0])
    plt.scatter(coordinates[bnd_ind, 0], coordinates[bnd_ind, 1],
                color='green', marker='.', s=50, label='Boundary Nodes')

    plt.grid()
    plt.axis('equal')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.legend()
    plt.title('NY DIMACS graph with multiple correlations regions', fontsize=14, fontweight='bold')
    plt.show()

def export_gr_file(graph, vertices_count, edges_count, filename, header_string):
    with open(filename, 'w') as file:
        # Write header
        file.write('c ' + header_string + '\n')

        # Write problem line (optional)
        file.write("p sp {} {}\n".format(vertices_count, edges_count))

        # Write edge lines
        for i in range(edges_count):
            file.write("a {} {} {}\n".format(graph[i, 0].astype(int), graph[i, 1].astype(int), graph[i, 2].astype(int)))

def export_coords_file(coordinates, vertices_count, new_coords_filename, header_string):
    with open(new_coords_filename, 'w') as file:
        # Write header
        file.write('c ' + header_string + '\n')

        # Write problem line (optional)
        file.write("p aux sp co {}\n".format(vertices_count))

        # Write nodes' coordinates
        for i in range(vertices_count):
            file.write("v {} {} {}\n".format(coordinates[i, 2].astype(int), coordinates[i, 0].astype(int), coordinates[i, 1].astype(int)))

def export_clusters_file(clusters, filename):
    with open(filename, 'w') as file:
        for i in range(clusters.shape[0]):
            file.write("{} {}\n".format(clusters[i, 0].astype(int), clusters[i, 1].astype(int)))

def testbench_A(input_filename, output_dir, corr_vec, samples_per_corr, path2ApexExe):
    # Reading the input graph structure
    edge_ind = 0
    with open(input_filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
                edges_count = int(edges_count_str)
                graph = np.zeros((edges_count, 3))
                continue

            # Inputting each edge data
            if line[0] == 'a':
                _, graph[edge_ind, 0], graph[edge_ind, 1], \
                    _ = line.rstrip('\n').split(' ')
                edge_ind += 1

        # Scanning through requested target correlations
        for target_corr in corr_vec:
            # Generating two costs vector with desired linear correlation
            distance_cost, time_cost = generate_correlated_vectors(edges_count, target_corr)

            # Shifting and scaling the data so all value will be non-negative
            min_value = 0
            max_value = 10e3
            distance_cost = np.clip(((abs(min(distance_cost)) + distance_cost) * 1000).astype(int), min_value, max_value)
            time_cost = np.clip(((abs(min(time_cost)) + time_cost) * 1000).astype(int), min_value, max_value)

            # debug printing
            new_correlation = pearson_correlation(distance_cost, time_cost)
            print(f'Current correlation is {new_correlation}')

            # Generating the new distance graph file
            new_distance_filename = output_dir + f'\Distance_corr_{int(round(target_corr*100))}.gr'
            with open(new_distance_filename, 'w') as file:
                # Write header
                file.write(f'c Distance graph adjusted to correlation = {target_corr}\n')

                # Write problem line (optional)
                file.write("p sp {} {}\n".format(vertices_count, edges_count))

                # Write edge lines
                for i in range(edges_count):
                    file.write("a {} {} {}\n".format(graph[i, 0].astype(int), graph[i, 1].astype(int), distance_cost[i]))

            # Generating the new time graph file
            new_time_filename = output_dir + f'\Time_corr_{int(round(target_corr*100))}.gr'
            with open(new_time_filename, 'w') as file:
                # Write header
                file.write(f'c Time graph adjusted to correlation = {target_corr}\n')

                # Write problem line (optional)
                file.write("p sp {} {}\n".format(vertices_count, edges_count))

                # Write edge lines
                for i in range(edges_count):
                    file.write("a {} {} {}\n".format(graph[i, 0].astype(int), graph[i, 1].astype(int), time_cost[i]))

            # Performing multiple A*pex invocations for the same correlation
            for i in range(samples_per_corr):
                print(f'Iteration {i}:')
                ready = False
                while not ready:
                    startNode = np.random.randint(min(graph[:, 0]), max(graph[:, 0]))
                    goalNode = np.random.randint(min(graph[:, 0]), max(graph[:, 0]))
                    if abs(startNode - goalNode) > vertices_count * 0.5:
                        ready = True
                log_file = f'{output_dir}\log_corr_{int(round(target_corr*100))}_iter_{i+1}.json'
                command_line = f'{path2ApexExe}  -m {new_distance_filename} {new_time_filename} \
                               -e 0 -s {startNode} -g {goalNode} -a Apex -o output.txt -l {log_file}'
                os.system(command_line)

def fix_json_file(filename):
    tmp_filename = filename + '_'
    blocks_number = 0
    current_block = 1
    with open(filename, "r") as f_in:
        with open(tmp_filename, "w") as f_out:
            f_in.readline()
            for line in f_in:
                if line.find('start_time') != -1:
                    blocks_number += 1

    if blocks_number == 1:
        return

    with open(filename, "r") as f_in:
        with open(tmp_filename, "w") as f_out:
            f_in.readline()
            for line in f_in:
                if current_block == blocks_number:
                    f_out.write(line)
                elif line[0] == ']':
                    current_block += 1

    os.rename(filename, filename + '_old')
    os.rename(tmp_filename, filename)

def analyze_corr_vs_optimality_ratio(output_dir, correlations, samples):
    results = {}
    for corr in correlations:
        results[corr] = np.array([])
        for sample in range(1, samples + 1):
            log_file = f'{output_dir}\log_corr_{int(round(corr * 100))}_iter_{sample}.json'
            with open(log_file, 'r') as file:
                log = json.load(file)
            n_solutions = log[0]['finish_info']['amount_of_solutions']
            solutions = np.zeros((n_solutions, 2))
            for i in range(n_solutions):
                solutions[i, :] = log[0]['finish_info']['solutions'][i]['full_cost']

            # Computing the width in respect to the 1st cost function
            best_cost2_sol_ind = np.argmin(solutions[:, 1])
            cost1_at_best_cost2 = solutions[best_cost2_sol_ind, 0]
            best_cost1 = min(solutions[:, 0])
            needed_epsilon_cost1 = (cost1_at_best_cost2 - best_cost1) / best_cost1

            # Computing the width in respect to the 2nd cost function
            best_cost1_sol_ind = np.argmin(solutions[:, 0])
            cost2_at_best_cost1 = solutions[best_cost1_sol_ind, 1]
            best_cost2 = min(solutions[:, 1])
            needed_epsilon_cost2 = (cost2_at_best_cost1 - best_cost2) / best_cost2

            needed_epsilon = max([needed_epsilon_cost1, needed_epsilon_cost2])

            if needed_epsilon < 5:
                results[corr] = np.append(results[corr], needed_epsilon)

    for corr in results.keys():
        plt.scatter(np.ones(results[corr].shape) * corr, results[corr], color='grey', marker='o', s=5)

    plt.plot(correlations, [max(results[corr]) for corr in results.keys()], color='red', linewidth=3)
    # plt.plot(correlations, [np.percentile(results[corr], 99) for corr in results.keys()], color='red', linewidth=3)

    # Adding labels and title
    plt.xlabel('Pearson Correlation Coefficient')
    plt.ylabel(r'Needed Approximation Factor ($\varepsilon$)')
    plt.title('Pareto-Optimal front width vs. Correlation')
    plt.grid()
    plt.show()

def pearson_correlation(vector1, vector2):
    # Standardize the vectors
    std_vector1 = (vector1 - np.mean(vector1)) / np.std(vector1)
    std_vector2 = (vector2 - np.mean(vector2)) / np.std(vector2)

    # Calculate the current correlation
    return np.corrcoef(std_vector1, std_vector2)[0, 1]

def generate_correlated_vectors(size, target_correlation):
    # Generate uncorrelated random variables
    random_variables = np.random.randn(2, size)

    # Create the covariance matrix from the target correlation
    covariance_matrix = np.array([[1, target_correlation], [target_correlation, 1]])

    # Perform Cholesky decomposition on the covariance matrix
    cholesky_matrix = np.linalg.cholesky(covariance_matrix)

    # Transform the uncorrelated variables to be correlated
    correlated_variables = np.dot(cholesky_matrix, random_variables)

    # Extract the correlated vectors, adding bias to avoid negative values
    vector1 = correlated_variables[0, :] - min(correlated_variables[0, :])
    vector2 = correlated_variables[1, :] - min(correlated_variables[1, :])

    return vector1, vector2

def generate_correlated_graph(gr_filename, out_distance_filename, out_time_filename, target_corr):
    # Reading the input graph
    edge_ind = 0
    with open(gr_filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
                edges_count = int(edges_count_str)
                graph = np.zeros((edges_count, 3))
                continue

            # Inputting each edge data
            if line[0] == 'a':
                _, graph[edge_ind, 0], graph[edge_ind, 1], \
                    _ = line.rstrip('\n').split(' ')
                edge_ind += 1

        # Generating two costs vector with desired linear correlation
        distance_cost, time_cost = generate_correlated_vectors(edges_count, target_corr)

        # Shifting and scaling the data so all value will be non-negative
        distance_cost = ((abs(min(distance_cost)) + distance_cost) * 1000).astype(int)
        time_cost = ((abs(min(time_cost)) + time_cost) * 1000).astype(int)

        # debug check new correlation
        new_correlation = pearson_correlation(distance_cost, time_cost)
        print(f'Correlation check is {new_correlation}')

        # Generating the new distance graph file
        with open(out_distance_filename, 'w') as file:
            # Write header
            file.write(f'c Distance graph adjusted to correlation = {target_corr}\n')

            # Write problem line (optional)
            file.write("p sp {} {}\n".format(vertices_count, edges_count))

            # Write edge lines
            for i in range(edges_count):
                file.write("a {} {} {}\n".format(graph[i, 0].astype(int), graph[i, 1].astype(int), distance_cost[i]))

        # Generating the new time graph file
        with open(out_time_filename, 'w') as file:
            # Write header
            file.write(f'c Time graph adjusted to correlation = {target_corr}\n')

            # Write problem line (optional)
            file.write("p sp {} {}\n".format(vertices_count, edges_count))

            # Write edge lines
            for i in range(edges_count):
                file.write("a {} {} {}\n".format(graph[i, 0].astype(int), graph[i, 1].astype(int), time_cost[i]))
'''
def plot_coords_graph(filename):
    with open(filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
'''

def testbench_B(distance_filename, time_filename, output_dir, samples, path2ApexExe):
    # Reading the input graph structure
    edge_ind = 0
    with open(distance_filename, "r") as f:
        f.readline()
        for line in f:
            # Inputting the graph dimensions and allocating data structures
            if line[0] == 'p':
                _, _, vertices_count_str, edges_count_str = line.rstrip('\n').split(' ')
                vertices_count = int(vertices_count_str)
                edges_count = int(edges_count_str)
                graph = np.zeros((edges_count, 3))
                continue

            # Inputting each edge data
            if line[0] == 'a':
                _, graph[edge_ind, 0], graph[edge_ind, 1], \
                    _ = line.rstrip('\n').split(' ')
                edge_ind += 1

        # Performing multiple A*pex invocations for the same correlation
        for i in range(samples):
            print(f'Iteration {i}:')
            ready = False
            while not ready:
                startNode = np.random.randint(min(graph[:, 0]), max(graph[:, 0]))
                goalNode = np.random.randint(min(graph[:, 0]), max(graph[:, 0]))
                ready = True
                # if abs(startNode - goalNode) > vertices_count * 0.5:
                #     ready = True
            log_file = f'{output_dir}\log_iter_{i+1}.json'
            command_line = f'{path2ApexExe}  -m {distance_filename} {time_filename} \
                           -e 0 -s {startNode} -g {goalNode} -a Apex -o output.txt -l {log_file}'
            os.system(command_line)

            # Analyzing the Pareto-set width
            with open(log_file, 'r') as file:
                log = json.load(file)
            n_solutions = log[0]['finish_info']['amount_of_solutions']
            solutions = np.zeros((n_solutions, 2))
            for i in range(n_solutions):
                solutions[i, :] = log[0]['finish_info']['solutions'][i]['full_cost']

            # Computing the width in respect to the 1st cost function
            best_cost2_sol_ind = np.argmin(solutions[:, 1])
            cost1_at_best_cost2 = solutions[best_cost2_sol_ind, 0]
            best_cost1 = min(solutions[:, 0])
            needed_epsilon_cost1 = (cost1_at_best_cost2 - best_cost1) / best_cost1

            # Computing the width in respect to the 2nd cost function
            best_cost1_sol_ind = np.argmin(solutions[:, 0])
            cost2_at_best_cost1 = solutions[best_cost1_sol_ind, 1]
            best_cost2 = min(solutions[:, 1])
            needed_epsilon_cost2 = (cost2_at_best_cost1 - best_cost2) / best_cost2

            needed_epsilon = max([needed_epsilon_cost1, needed_epsilon_cost2])
            print('=========================================================================================')
            print(f'-------------> Epsilons = ({round(needed_epsilon_cost1, 2)},{round(needed_epsilon_cost2, 2)}), Needed Epsilon is {round(needed_epsilon, 2)}')
            print('=========================================================================================')

class LogAnalysis(object):
    def __init__(self, filename):
        with open(filename, 'r') as file:
            log = json.load(file)
        n_solutions = log[0]['finish_info']['amount_of_solutions']
        self.solutions = np.zeros((n_solutions, 2))
        for i in range(n_solutions):
            self.solutions[i, :] = log[0]['finish_info']['solutions'][i]['full_cost']

    def plot_pareto_set(self):
        plt.scatter(self.solutions[:, 0] / max(self.solutions[:, 0]), self.solutions[:, 1] / max(self.solutions[:, 1]), color='blue', marker='o')

        # Adding labels and title
        plt.xlabel('Distance cost')
        plt.ylabel('Time cost')
        plt.title('Set of Pareto-Optimal solutions')

        plt.grid()

        plt.show()



