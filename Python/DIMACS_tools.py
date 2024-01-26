import json
import time
from scipy.sparse import coo_matrix, csgraph
from scipy.sparse.linalg import eigsh
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import os
import networkx as nx
import numpy as np

def pareto_width(high_cost, low_cost):
    return (high_cost - low_cost) / low_cost

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

def graph_compress(distance_filename, time_filename, eps):
    c1_graph, vertices_count = load_graph(distance_filename)
    c2_graph, _ = load_graph(time_filename)

    # Truncating the graph
    n_edges = 1000
    c1_graph = c1_graph[0:n_edges, :]
    c2_graph = c2_graph[0:n_edges, :]

    # Building adjacency matrix
    N = c1_graph.shape[0]
    adj_matrix = np.zeros((N, N))
    for i in range(N):
        start_node = int(c1_graph[i, 0])
        end_node = int(c1_graph[i, 1])
        cost = int(c1_graph[i, 2])
        adj_matrix[start_node, end_node] = cost

    t1 = time.time()
    U, S, Vh = np.linalg.svd(adj_matrix, full_matrices=False)
    t2 = time.time()
    delta = t2-t1
    print(delta)
    '''
    # Creating the ratio graph
    mean_c1 = np.mean(c1_graph[:, 2])
    mean_c2 = np.mean(c2_graph[:, 2])
    ratio_graph = np.zeros(c1_graph.shape)
    ratio_graph[:, 0:2] = c1_graph[:, 0:2]
    ratio_graph[:, 2] = (c1_graph[:, 2] - mean_c1) / (c2_graph[:, 2] - mean_c2)

    # Plotting
    plt.hist(ratio_graph[:, 2], bins=1000, color='skyblue', edgecolor='black')

    # Adding labels and title
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title('Histogram Example')

    # Display the plot
    plt.show()

    # Call Spectral Clustering
    num_clusters = 2

    # Perform spectral clustering
    nodes, cluster_labels = spectral_clustering(ratio_graph, num_clusters)
    '''
    print('Building graph for networkx...')
    G = nx.DiGraph()
    for i in range(c1_graph.shape[0]):
        start_node = int(c1_graph[i, 0])
        end_node = int(c1_graph[i, 1])
        cost = int(c1_graph[i, 2])
        G.add_edge(start_node, end_node, weight=cost)

    print('Running Floyd-Warshall...')
    result_matrix = nx.floyd_warshall(G)

    print('done')

def plot_graph(adjacency_filename, coords_filename):
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
    coordinates = {}
    with open(coords_filename, "r") as f:
        f.readline()
        for line in f:
            if line[0] == 'v':
                _, node, x, y = line.rstrip('\n').split(' ')
                coordinates[int(node)] = np.array([x,y]).astype(int)

    coordinates_matrix = np.array(list(coordinates.values()))
    plt.scatter(coordinates_matrix[:, 0], coordinates_matrix[:, 1], color='blue', marker='.', s=5)
    plt.grid()
    plt.axis('equal')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()

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

            # Finding the best cost1 solution
            best_cost1_sol_ind = np.argmin(solutions[:, 0])
            cost2_at_best_cost1 = solutions[best_cost1_sol_ind, 1]
            best_cost2 = min(solutions[:, 1])
            needed_epsilon = (cost2_at_best_cost1 - best_cost2) / best_cost2
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

    # Extract the correlated vectors
    vector1 = correlated_variables[0, :]
    vector2 = correlated_variables[1, :]

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



