#include "Utils.h"
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <chrono>
#include <thread>
#include <unordered_map>

using namespace std;

int main(int argc, char** argv)
{
    namespace po = boost::program_options;
    std::vector<string> objective_files;
    std::string cluster_map_file;

    // ==========================================================================
    // Parsing the supported line arguments options
    // ==========================================================================
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("map,m", po::value< std::vector<string> >(&objective_files)->multitoken(), "files for edge weight")
        ("clusters,c", po::value<std::string>()->default_value(""), "clusters mapping file")
        ("eps,e", po::value<double>()->default_value(0), "approximation factor")
        ("output,o", po::value<std::string>()->required(), "Name of the output file")
        ("logging_file,l", po::value<std::string>()->default_value(""), "logging file")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    po::notify(vm);

    // ==========================================================================
    // Loading the objectives files
    // ==========================================================================
    size_t graph_size;
    std::vector<Edge> edges;

    cout << "Loading the following objective graphs:" << endl;
    for (auto file : objective_files) {
        cout << file << std::endl;
    }

    if (load_gr_files(objective_files, edges, graph_size) == false) {
        std::cout << "Failed to load gr files" << std::endl;
        return -1;
    }

    std::cout << "Graph Size: " << graph_size << std::endl;

    // Normalizing costs to [0,1] segment
    //normalize_edge_costs(edges);

    // Constructing the graph's adjacency matrix
    AdjacencyMatrix graph(graph_size, edges);
    AdjacencyMatrix inv_graph(graph_size, edges, true);

    // ==========================================================================
    // Loading the clustering analysis
    // ==========================================================================
    ClusterMapping clusters_map;
    load_clusters_mapping(vm["clusters"].as<std::string>(), clusters_map);

    cout << "Number of clusters is: " << clusters_map.size() << std::endl;

    // ==========================================================================
    // Computing all-pairs shortest-path between all boundary nodes
    // ==========================================================================
    std::vector<double> approx_factor = { 0.1, 0.1 };
    for (int cluster_id = 0; cluster_id < clusters_map.size(); cluster_id++)
    {
        // Creating a lookup table of all nodes inside the cluster
        std::unordered_set<int> lookupTable;
        build_lookup_table(lookupTable, clusters_map[cluster_id]);

        // Collecting only the boundary nodes of the current cluster
        std::vector<int> boundary_nodes;
        get_boundary_nodes(boundary_nodes, graph, lookupTable);

        //boundary_nodes.resize(100);

        // DEBUG EXPORT TO FILE
        /*
        std::ofstream outFile(R"(D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\boundary_nodes.txt)");
        for (int node : boundary_nodes)
        {
            outFile << node << std::endl;
        }
        outFile.close();
        */

        // Computing shortest-paths between all boundary nodes pairs
        //int nodes_count = boundary_nodes.size();
        //AllPairsCostsTensor costs(nodes_count);
        //all_pairs_shortest_paths(boundary_nodes, graph, costs);
        
        std::vector<int> cluster_nodes(lookupTable.begin(), lookupTable.end());

        // debug
        /*
        boundary_nodes.clear();
        for (int node : cluster_nodes)
        {
            boundary_nodes.push_back(node);
        }
        boundary_nodes.resize(5000);
        */

        //int bounary_nodes_count = boundary_nodes.size();
        //AllPairsCostsTensor costs(bounary_nodes_count);

        
        all_pairs_shortest_paths(cluster_nodes, boundary_nodes, graph, approx_factor);
    }
    

    /*
    // Code for displaying progress bar
    int total = 100;
    for (int i = 0; i <= total; ++i) {
        updateProgressBar(i, total);
        // Simulate some work
        this_thread::sleep_for(chrono::milliseconds(50));
    }
    cout << endl;
    */

    return 0;
}