#include "Utils.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <limits>
#include <chrono>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#define MAX_PATH_LENGTH 1e6

int get_edge_id(int source, int target, AdjacencyMatrix graph)
{
    std::vector<Edge> Edges = graph[source];
    for(int edge_id = 0; edge_id < Edges.size(); edge_id++)
    {
        if (Edges[edge_id].target == target)
        {
            return edge_id;
        }
    }
    
    return -1;
}


AdjacencyMatrix::AdjacencyMatrix(size_t graph_size, std::vector<Edge>& edges, bool inverse)
    : matrix((graph_size + 1), std::vector<Edge>()), graph_size(graph_size) {

    num_of_objectives = edges[0].cost.size();

    for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
        if (inverse) {
            this->add(iter->inverse());
        }
        else {
            this->add((*iter));
        }
    }
}


size_t AdjacencyMatrix::get_num_of_objectives() const {
    return num_of_objectives;
}

void AdjacencyMatrix::add(Edge edge) {
    (this->matrix[edge.source]).push_back(edge);
}


size_t AdjacencyMatrix::size() const { return this->graph_size; }


const std::vector<Edge>& AdjacencyMatrix::operator[](size_t vertex_id) const {
    return this->matrix.at(vertex_id);
}



void split_string(std::string string, std::string delimiter, std::vector<std::string>& results)
{
    size_t first_delimiter;

    while ((first_delimiter = string.find_first_of(delimiter)) != string.npos) {
        if (first_delimiter > 0) {
            results.push_back(string.substr(0, first_delimiter));
        }
        string = string.substr(first_delimiter + 1);
    }

    if (string.length() > 0) {
        results.push_back(string);
    }
}

bool load_gr_files(std::vector<std::string> gr_files, std::vector<Edge>& edges_out, size_t& graph_size) {
    size_t          max_node_num = 0;
    for (auto gr_file : gr_files) {
        std::ifstream file(gr_file.c_str());

        if (file.is_open() == false) {
            std::cerr << "cannot open the gr file " << gr_file << std::endl;
            return false;
        }

        std::string line;
        int idx_edge = 0;
        while (file.eof() == false) {
            std::getline(file, line);

            if (line == "") {
                break;
            }

            std::vector<std::string> decomposed_line;
            split_string(line, " ", decomposed_line);

            std::string type = decomposed_line[0];
            if ((std::strcmp(type.c_str(), "c") == 0) || (std::strcmp(type.c_str(), "p") == 0)) {
                continue; //comment or problem lines, not part of the graph
            }

            if (std::strcmp(type.c_str(), "a") == 0) { //arc
                if (idx_edge <= (int)edges_out.size() - 1) {
                // original line was : if (idx_edge < (int)edges_out.size() - 1) {
                    if (
                        (stoul(decomposed_line[1]) != edges_out[idx_edge].source) ||
                        (stoul(decomposed_line[2]) != edges_out[idx_edge].target)) {
                        // arc_sign src dest should be same in both files
                        std::cerr << "file inconsistency" << std::endl;
                        return false;
                    }
                    edges_out[idx_edge].cost.push_back(std::stoul(decomposed_line[3]));
                }
                else {
                    Edge e(std::stod(decomposed_line[1]),
                        std::stod(decomposed_line[2]),
                        { std::stod(decomposed_line[3]) });
                    /*
                    Edge e(std::stoul(decomposed_line[1]),
                        std::stoul(decomposed_line[2]),
                        { std::stoul(decomposed_line[3]) });
                    */
                    edges_out.push_back(e);
                    max_node_num = std::max({ max_node_num, e.source, e.target });
                }
            }
            idx_edge++;
        }
        file.close();
    }
    graph_size = max_node_num;
    return true;
}

bool load_clusters_mapping(std::string filename, ClusterMapping& clusters_map)
{
    std::ifstream file(filename.c_str());

    if (file.is_open() == false) 
    {
        std::cerr << "cannot open the gr file " << filename << std::endl;
        return false;
    }

    std::string line;
    int current_cluster_id = -1;
    std::vector<int> current_cluster_nodes;
    bool new_cluster_data = false;

    // debug
    int max_nodes = 20;

    while (file.eof() == false)
    {
        std::getline(file, line);

        std::vector<std::string> decomposed_line;
        split_string(line, " ", decomposed_line);

        // Asserting that the string line holds data
        if (decomposed_line.size() == 0)
        {
            continue;
        }
        
        if (stoul(decomposed_line[0]) != current_cluster_id)
        {
            current_cluster_id = stoul(decomposed_line[0]);
            if (current_cluster_nodes.size() > 0)
            {
                // Start of the new cluster: adding the previous cluster data to the general vector and progressing to next cluster
                clusters_map.push_back(current_cluster_nodes);
            }
            current_cluster_nodes.clear();
            current_cluster_nodes.push_back(stoul(decomposed_line[1]));
        }
        else
        {
            if(0 * current_cluster_nodes.size() < max_nodes)
            { 
                // Appending additional data to the current cluster
                current_cluster_nodes.push_back(stoul(decomposed_line[1]));
            }
            // Appending additional data to the current cluster
            //current_cluster_nodes.push_back(stoul(decomposed_line[1]));
        }
    }
    // Adding the last data row
    clusters_map.push_back(current_cluster_nodes);
    
    return true;
}

bool build_lookup_table(std::unordered_set<int>& lut, std::vector<int> cluster_nodes)
{
    for (int node_id : cluster_nodes)
    {
        lut.insert(node_id);
    }
    return true;
}

bool get_boundary_nodes(std::vector<int>& boundary_nodes, AdjacencyMatrix graph, const std::unordered_set<int>lookupTable)
{
    for (int node_id : lookupTable)
    {
        // Scanning over all outgoing edges
        for (Edge outGoingEdge : graph[node_id])
        {
            // Retrieving all the target nodes that are connected to the current node
            int targetNode = outGoingEdge.target;
            if (lookupTable.find(targetNode) == lookupTable.end()) 
            {
                // The target node is not part of the cluster, thus the current node is a boundary node
                boundary_nodes.push_back(node_id);
                break;
            }
        }
    }

    return true;
}

std::vector<int> reconstructPath_floydWarshall(int start, int end, std::vector<std::vector<int>>& next)
{
    std::vector<int> path;
    path.push_back(start);

    while (start != end) {
        start = next[start][end];
        path.push_back(start);
    }

    return path;
}

void floydWarshall(std::vector<std::unordered_map<int, int>>& graph, std::vector<std::vector<int>>& next) {
    int V = graph.size();

    // Initialize next array for path reconstruction
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            next[i][j] = j;
        }
    }

    // Apply Floyd-Warshall algorithm
    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) 
            {
                if (graph[i].find(k) != graph[i].end() && graph[k].find(j) != graph[k].end() &&
                    graph[i][k] + graph[k][j] < graph[i][j]) {
                    graph[i][j] = graph[i][k] + graph[k][j];
                    next[i][j] = next[i][k]; // Update next node on the shortest path
                }
            }
        }
    }
}
void reconstructPath_Dijkstra(int source, int target, std::vector<int> predecessors, 
    std::vector<int>& path, int* path_length)
{
    *path_length = 0;
    for (int v = target; v != source; v = predecessors[v]) 
    {
        path[*path_length] = v;
        *path_length = *path_length + 1;
        if (*path_length >= MAX_PATH_LENGTH)
        {
            std::cout << "Warning! reconstructPath_Dijkstra is not working properly!" << std::endl;
        }
    }
    path[*path_length] = source;
    *path_length = *path_length + 1;
}

/*
void reconstructPath_Dijkstra(int source, int target, std::vector<int> predecessors,
    std::vector<int>& path)
{
    for (int v = target; v != source; v = predecessors[v])
    {
        auto it = path.begin();
        path.insert(it, v);
        if (path.size() > 1000000)
        {
            std::cout << "Warning! reconstructPath_Dijkstra is not working properly!" << std::endl;
        }
    }
    auto it = path.begin();
    path.insert(it, source);
}
*/

bool is_path_contractable(CrossObjectiveCost crossCostsMat, std::vector<double> approx_factor)
{
    double width;
    int n_objectives = approx_factor.size();
    for (int objective_id = 0; objective_id < n_objectives; objective_id++)
    {
        // Same row share the same objective function under minimization
        double min_cost = crossCostsMat[objective_id][objective_id];
        for (int cross_obj_id = 0; cross_obj_id < n_objectives; cross_obj_id++)
        {
            if (objective_id != cross_obj_id)
            {
                width = ((double)crossCostsMat[cross_obj_id][objective_id] - min_cost) / min_cost;

                if (width < approx_factor[objective_id])
                {
                    // This condition needs to hold only once for the entire path to be contractable
                    return true;
                }
            }
        }
    }
    return false;
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All-Pairs Shortest-Paths solver based on Dijkstra algorithm O(VlogV + E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

bool all_pairs_shortest_paths(std::vector<int> cluster_nodes, std::vector<int> boundary_nodes,
    AdjacencyMatrix entire_graph, AllPairsCostsTensor& costs, std::vector<double> approx_factor)
{
    int cluster_nodes_count = cluster_nodes.size();
    int boundary_nodes_count = boundary_nodes.size();
    int objectives_count = entire_graph.get_num_of_objectives();

    // Asserting that number of boundary nodes is not greater than the total number of cluster's nodes
    if (boundary_nodes_count > cluster_nodes_count)
    {
        std::cout << "Warning! Boundary nodes count is bigger than cluster nodes count" << std::endl;
    }

    // Hash table between node id to adjacency matrix id
    std::unordered_map<int, int> node_id_to_adj_mat_id;
    for (int k = 0; k < cluster_nodes_count; k++)
    {
        node_id_to_adj_mat_id[cluster_nodes[k]] = k;
    }

    // Hash table between node id to cross-costs tensor matrix id
    std::unordered_map<int, int> node_id_to_tensor_mat_id;
    for (int k = 0; k < boundary_nodes_count; k++)
    {
        node_id_to_tensor_mat_id[boundary_nodes[k]] = k;
    }

    // Hash table between node id and its node successors
    std::unordered_map<int, std::unordered_map<int, int>> succesor_node_to_edge_id;
        
    // Allocating memory for the cross-costs tensor
    std::cout << "Costs output structure memory allocation ... ";
    costs.resize(boundary_nodes_count); // Resize the matrix to have 2 rows
    for (int i = 0; i < boundary_nodes_count; i++)
    {
        costs[i].resize(boundary_nodes_count); // Resize each row to have 2 elements
    }

    for (int i = 0; i < boundary_nodes_count; i++)
    {
        for (int j = 0; j < boundary_nodes_count; j++)
        {
            costs[i][j].resize(objectives_count);
            for (int k = 0; k < objectives_count; k++)
            {
                costs[i][j][k].resize(objectives_count, -1);
            }
        }
    }
    std::vector<std::vector<int>> pathLengths(boundary_nodes_count, std::vector<int>(boundary_nodes_count));
    std::cout << "Done" << std::endl;

    // Defining the graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;

    std::vector<int> path(MAX_PATH_LENGTH);
    int path_length;

    std::cout << "Allocating memory for distances and predecessors structures ... ";
    // Define the vector to store distances
    std::vector<std::vector<double>> distances(cluster_nodes_count, std::vector<double>(cluster_nodes_count));

    // Define the vector to store predecessors
    std::vector<std::vector<int>> predecessors(cluster_nodes_count, std::vector<int>(cluster_nodes_count));

    std::cout << "Done" << std::endl;

    // Computing all-pairs shortest path for every objective function
    for (int objective_id = 0; objective_id < objectives_count; objective_id++)
    {
        std::cout << "Constructing graph for objective " << (objective_id + 1) << " ... ";
        // Define the graph object
        Graph g(cluster_nodes_count);

        // Building the adjacency matrix for the specific objective function
        for (int i = 0; i < cluster_nodes_count; i++)
        {
            std::vector<Edge> edges = entire_graph[cluster_nodes[i]];
            std::unordered_map<int, int> edges_hash;
            for(int k = 0; k < edges.size() ; k++)
            {
                // Check if both target node is part of the cluster
                if (node_id_to_adj_mat_id.find(edges[k].target) != node_id_to_adj_mat_id.end())
                {
                    // Key found
                    int adj_mat_source_id = i;
                    int adj_mat_target_id = node_id_to_adj_mat_id[edges[k].target];
                    double cost = edges[k].cost[objective_id];
                    boost::add_edge(adj_mat_source_id, adj_mat_target_id, cost, g);
                    if (cost < 0)
                    {
                        std::cout << "Found negative edge weight of cost= " << cost << " in start vertex = " << i << std::endl;
                    }
                    if (objective_id == 0)
                    {
                        // Hash should be populated only once
                        edges_hash[adj_mat_target_id] = k;
                    }
                }
            }
            if (objective_id == 0)
            {
                // Hash should be populated only once
                succesor_node_to_edge_id[i] = edges_hash;
            }
            
        }
        std::cout << "Done" << std::endl;

        // Define a property map for storing edge weights
        typedef boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMap;
        EdgeWeightMap weightMap = boost::get(boost::edge_weight, g);

        // Call the All-Pairs Dijkstra algorithm
        std::cout << "All-Pairs Dijkstra invoked for objective number " << (objective_id + 1) << " ... ";
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<double> dist(boost::num_vertices(g));
        std::vector<int> pred(boost::num_vertices(g));

        for (int source = 0; source < boundary_nodes_count; ++source) 
        {
            int source_id = node_id_to_adj_mat_id[boundary_nodes[source]];
            
            boost::dijkstra_shortest_paths(g, source_id,
                boost::distance_map(boost::make_iterator_property_map(dist.begin(), boost::get(boost::vertex_index, g)))
                .predecessor_map(boost::make_iterator_property_map(pred.begin(), boost::get(boost::vertex_index, g))));

            for (int target = 0; target < cluster_nodes_count; ++target)
            {
                int target_id = node_id_to_adj_mat_id[cluster_nodes[target]];

                // Checking that target is reachable from source by asserting that
                // the target's predecessor is NOT itself (marks unreachable node)
                if (pred[target_id] != target_id)
                {
                    distances[source_id][target_id] = dist[target_id];
                    predecessors[source_id][target_id] = pred[target_id];
                }
                else
                {
                    // target is not reachable, marking by a negative number
                    distances[source_id][target_id] = -1;
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Done" << std::endl;
        std::cout << "Finished in " << duration.count() << " [sec]" << std::endl;

        // Updating the cross-costs tensor data structure
        std::cout << "Updating cross-costs tensor matrix ... ";
        start = std::chrono::high_resolution_clock::now();
        for (int source = 0; source < boundary_nodes_count; source++)
        {
            int adj_mat_source_id = node_id_to_adj_mat_id[boundary_nodes[source]];

            for (int target = 0; target < boundary_nodes_count; target++)
            {
                if (source == target)
                {
                    // Avoid equal start and target nodes
                    continue;
                }
                        
                int adj_mat_target_id = node_id_to_adj_mat_id[boundary_nodes[target]];

                // Checking that source and target are reachable
                bool target_reachable = distances[adj_mat_source_id][adj_mat_target_id] >= 0;

                // Updating the optimal cost considering the current objective id
                costs[source][target][objective_id][objective_id] =
                        distances[adj_mat_source_id][adj_mat_target_id];

                // Updating the rest of the costs (not under the optimization)
                for (int cross_obj_id = 0; cross_obj_id < objectives_count; cross_obj_id++)
                {
                    if (objective_id != cross_obj_id)
                    {
                        if (target_reachable)
                        {
                            costs[source][target][objective_id][cross_obj_id] = 0;
                            reconstructPath_Dijkstra(adj_mat_source_id, adj_mat_target_id, 
                                predecessors[adj_mat_source_id], path, &path_length);

                            pathLengths[source][target] = path_length;

                            for (int i = (path_length - 1); i > 0; i--)
                            {
                                int out_ind = cluster_nodes[path[i]];
                                int in_ind  = cluster_nodes[path[i - 1]];

                                int adj_out = node_id_to_adj_mat_id[out_ind];
                                int adj_in  = node_id_to_adj_mat_id[in_ind];
                                int edge_id = succesor_node_to_edge_id[adj_out][adj_in];
                                if (edge_id < 0)
                                {
                                    std::cout << "Warning! get_edge_id returned (-1)!" << std::endl;
                                }
                                costs[source][target][objective_id][cross_obj_id] +=
                                    entire_graph[out_ind][edge_id].cost[cross_obj_id];
                            }
                        }
                        else
                        {
                            // Target is not reachable, marking by a (-1)
                            costs[source][target][objective_id][cross_obj_id] = -1;
                        }
                    }
                }
            }
        }
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "Done" << std::endl;
        std::cout << "Finished in " << duration.count() << " [sec]" << std::endl;
    }
    
    // Analysis of pareto-width for edge contraction
    std::cout << "Pareto-set analysis for edge contraction ... ";
    auto start = std::chrono::high_resolution_clock::now();
    int reachable_paths = 0;
    int contractable_paths = 0;
    for (int source = 0; source < boundary_nodes_count; source++)
    {
        int adj_mat_source_id = node_id_to_adj_mat_id[boundary_nodes[source]];
        for (int target = 0; target < boundary_nodes_count; target++)
        {
            // Skipping degenerated paths
            if (source == target)
            {
                continue;
            }

            // Checking that source and target are reachable
            int adj_mat_target_id = node_id_to_adj_mat_id[boundary_nodes[target]];
            bool target_reachable = distances[adj_mat_source_id][adj_mat_target_id] >= 0;
            if (!target_reachable)
            {
                continue;
            }
            reachable_paths++;

            CrossObjectiveCost crossCostsMat = costs[source][target];
            // debug
            if (crossCostsMat[0][0] > crossCostsMat[1][0] || 
                crossCostsMat[1][1] > crossCostsMat[0][1])
            {
                std::cout << "problem\n";
            }
            bool pathContractable = is_path_contractable(crossCostsMat, approx_factor);
            if (pathContractable)
            {
                contractable_paths++;
            }
            //std::cout << "Path Length = " << pathLengths[source][target] << " Is_Contractable = " << pathContractable << std::endl;
        }
    }
    double compression_ratio = (double)contractable_paths / (double)reachable_paths;
    std::cout << "Reachable paths count = " << reachable_paths << std::endl;
    std::cout << "Contractable paths count = " << contractable_paths << " (Compression Ratio = " << compression_ratio << ")" << std::endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Done" << std::endl;
    std::cout << "Finished in " << duration.count() << " [sec]" << std::endl;
    
    return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Normalizing edge costs to [0,1] segment
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool normalize_edge_costs(std::vector<Edge>& edges)
{
    double max_double = std::numeric_limits<double>::max();
    double min_double = -std::numeric_limits<double>::max();

    int n_costs = edges[0].cost.size();
    std::vector<double> min_vals(n_costs, max_double);
    std::vector<double> max_vals(n_costs, min_double);
    
    // Iterating over all edges for min and max values extraction cost-wise
    for (int edge_id = 0; edge_id < edges.size(); edge_id++)
    {
        
        for (int cost_id = 0; cost_id < n_costs; cost_id++)
        {
            if (edges[edge_id].cost[cost_id] < min_vals[cost_id])
            {
                min_vals[cost_id] = edges[edge_id].cost[cost_id];
            }
            if (edges[edge_id].cost[cost_id] > max_vals[cost_id])
            {
                max_vals[cost_id] = edges[edge_id].cost[cost_id];
            }
        }
    }

    // Normalizing each cost according to its corresponding min and max values range
    for (int edge_id = 0; edge_id < edges.size(); edge_id++)
    {
        for (int cost_id = 0; cost_id < n_costs; cost_id++)
        {
            edges[edge_id].cost[cost_id] =
                ((double)edges[edge_id].cost[cost_id] - min_vals[cost_id]) / (max_vals[cost_id] - min_vals[cost_id]);
        }
    }
    return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//All-Pairs Shortest-Paths solver based on Floyd-Warshall algorithm O(V^3)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
bool all_pairs_shortest_paths_obsolete(std::vector<int> cluster_nodes, AdjacencyMatrix entire_graph, 
    AllPairsCostsTensor& costs)
{
    int nodes_count = cluster_nodes.size();
    int objectives_count = entire_graph.get_num_of_objectives();

    // Building a hash table between the original nodes vertices ids and their sequential index among all cluster nodes
    std::unordered_map<int, int> global_to_local_index;
    for (int k = 0; k < nodes_count; k++)
    {
        global_to_local_index[cluster_nodes[k]] = k;
    }

    // Reshaping the costs output argument
    costs.resize(nodes_count); // Resize the matrix to have 2 rows
    for (int i = 0; i < nodes_count; i++)
    {
        costs[i].resize(nodes_count); // Resize each row to have 2 elements
    }

    for (int i = 0; i < nodes_count; i++)
    {
        for (int j = 0; j < nodes_count; j++)
        {
            costs[i][j].resize(objectives_count);
            for (int k = 0; k < objectives_count; k++)
            {
                costs[i][j][k].resize(objectives_count);
            }
        }
    }

    // Computing all-pairs shortest path for every objective function
    for (int objective_id = 0; objective_id < objectives_count; objective_id++)
    {
        // Initialize the needed adjacency lists
        std::vector<std::unordered_map<int, int>> graph(nodes_count);
        std::vector<std::vector<int>> next(nodes_count, std::vector<int>(nodes_count, -1)); // Predecessors array

        // Building the adjacency matrix for the specific objective function
        for (int i = 0; i < nodes_count; i++)
        {
            std::vector<Edge> edges = entire_graph[cluster_nodes[i]];
            for (Edge edge : edges)
            {
                // Check if both target node is part of the cluster
                if (global_to_local_index.find(edge.target) != global_to_local_index.end())
                {
                    // Key found
                    int adj_mat_source_id = i;
                    int adj_mat_target_id = global_to_local_index[edge.target];
                    graph[adj_mat_source_id][adj_mat_target_id] = edge.cost[objective_id];
                }
            }
        }

        // Call the Floyd-Warshall algorithm
        std::cout << "Floyd-Warshall invoked for objective number " << (objective_id + 1) << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        
        floydWarshall(graph, next);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Finished in " << duration.count() << " [sec]" << std::endl;

        // Updating the output argument
        for (int source = 0; source < nodes_count; source++)
        {
            int entire_graph_source_id = cluster_nodes[source];

            for (int target = 0; target < nodes_count; target++)
            {
                // Updating the optimal cost considering the current objective id
                costs[source][target][objective_id][objective_id] = graph[source][target];

                // Updating the rest of the costs (not under the optimization)
                for (int cross_obj_id = 0; cross_obj_id < objectives_count; cross_obj_id++)
                {
                    if (objective_id != cross_obj_id)
                    {
                        costs[source][target][objective_id][cross_obj_id] = 0;
                        std::vector<int> path = reconstructPath_floydWarshall(source, target, next);
                        for (int i = 0; i < path.size() - 1; i++)
                        {
                            int out_ind = cluster_nodes[path[i]];
                            int in_ind = cluster_nodes[path[i + 1]];
                            int edge_id = get_edge_id(out_ind, in_ind, entire_graph);
                            costs[source][target][objective_id][cross_obj_id] += 
                                entire_graph[entire_graph_source_id][edge_id].cost[cross_obj_id];
                        }
                    }
                }
            }
        }
    }

    return true;
}
*/