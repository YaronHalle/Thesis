from DIMACS_tools import *
import numpy as np

if __name__ == '__main__':
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Fixing multiple data JSON files
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # directory = r'D:\Thesis\DIMCAS\NY_correlated'
    # json_files = [file for file in os.listdir(directory) if file.endswith('.json')]
    # for filename in json_files:
    #     print(filename)
    #     fix_json_file(directory + '\\' + filename)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plot path length vs contractability histogram
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    stats_file = r"D:\Thesis\A-pex\multiobj\CorrelationClustering\pathLength_vs_Contractability.csv"
    # plot_path_length_vs_contractability(stats_file)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # la = LogAnalysis(r'D:\Thesis\DIMCAS\NY_correlated\log_corr_0_iter_7.json')
    # la.plot_pareto_set()

    gr_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Generating a single correlated instance graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # out_distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY_rho060.gr'
    # out_time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY_rho060.gr'
    # target_corr = 0.6
    out_distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY_rho000.gr'
    out_time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY_rho000.gr'
    target_corr = 0
    # generate_correlated_graph(gr_filename, out_distance_filename, out_time_filename, target_corr)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Testbench for A*pex on different correlations
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output_dir = r'D:\Thesis\DIMCAS\NY_correlated'
    corr_vec = np.hstack([np.arange(-0.99, 0, 0.01), np.arange(0, 1, 0.01)])
    samples_per_corr = 500
    path2ApexExe = r'D:\Thesis\A-pex\multiobj\x64\Release\multiobj.exe'

    # testbench_A(gr_filename, output_dir, corr_vec, samples_per_corr, path2ApexExe)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting testbench A results
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyze_corr_vs_optimality_ratio(output_dir, corr_vec, samples_per_corr)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Testbench for A*pex on original graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output_dir = r'D:\Thesis\DIMCAS\NY_original_MC'
    samples = 1000
    distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    # testbench_B(distance_filename, time_filename, output_dir, samples, path2ApexExe)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Generating a multi-correlated graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\USA-road-d.NY.co'
    # distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    distance_filename = r"D:\Thesis\DIMCAS\NY_correlated\Distance_corr_96.gr"
    # distance_filename = r"D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Rayleigh_Distance_corr_96.gr"

    # time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    time_filename = r"D:\Thesis\DIMCAS\NY_correlated\Time_corr_96.gr"
    # time_filename = r"D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Rayleigh_Time_corr_96.gr"

    clusters_count = 40
    # clusters_count = 0
    new_distance_gr_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Distance_NY_096_corr_clustered.gr'
    new_time_gr_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Time_NY_096_corr_clustered.gr'
    new_coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\NY_096_corr_clustered.co'
    clusters_metafile = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\NY_096_corr_clusters_metafile.txt'
    # min_lat = 4.058e7
    # max_lat = 4.058e7 + (4.08e7 - 4.058e7) * 0.75
    # min_lon = -7.41e7
    # max_lon = -7.41e7 + (-7.36e7 + 7.41e7) * 0.75
    min_lat = 4.058e7
    max_lat =  4.1e7  # 4.1e7  4.08e7
    min_lon = -7.41e7
    max_lon = -7.36e7  # -7.31e7
    # generate_multiple_correlated_graph(distance_filename, time_filename, coords_filename,
    #                                    clusters_count, new_distance_gr_filename,
    #                                    new_time_gr_filename, new_coords_filename, clusters_metafile,
    #                                    min_lon, max_lon, min_lat, max_lat)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting the graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot_graph(gr_filename, coords_filename, clusters_metafile)
    # plot_graph(new_distance_gr_filename, new_coords_filename, clusters_metafile)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Graph Clustering
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eps = 0.01
    # graph_agglomerative_clustering(new_distance_gr_filename, new_time_gr_filename, eps)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Generating enhanced graphs supplemented with contracted edges
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # original_distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    original_distance_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Distance_NY_096_corr_clustered.gr'

    # original_time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    original_time_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Time_NY_096_corr_clustered.gr'

    new_distance_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Distance_NY_096_corr_clustered_with_contraction.gr'
    new_time_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Time_NY_096_corr_clustered_with_contraction.gr'

    contracted_edges_file = r'D:\Thesis\A-pex\multiobj\CorrelationClustering\Contracted_Edges.csv'

    augmentation_ratio = 0.3

    # add_contracted_edges_to_graphs(original_distance_filename, original_time_filename,
    #                                new_distance_filename, new_time_filename, contracted_edges_file,
    #                                augmentation_ratio)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Comparing A*pex performance with/without contraction
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    original_distance_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Distance_NY_096_corr_clustered.gr'
    original_time_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Time_NY_096_corr_clustered.gr'
    coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\NY_096_corr_clustered.co'
    cntrcted_distance_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Distance_NY_096_corr_clustered_with_contraction.gr'
    cntrcted_time_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\Time_NY_096_corr_clustered_with_contraction.gr'
    clusters_metafile = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\NY_096_corr_clusters_metafile.txt'
    # output_dir = r'D:\Thesis\DIMCAS\NY_original_MC'
    output_dir = r"D:\Thesis\DIMCAS\NY_original_check3"
    samples = 500
    epsilon = 0.01
    testbench_D(original_distance_filename, original_time_filename, cntrcted_distance_filename, cntrcted_time_filename,
                coords_filename, clusters_metafile, output_dir, samples, path2ApexExe, epsilon)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting testbench_D results
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyze_apex_performance(output_dir, samples)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Comparing A*pex performance with heuristics inflation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    coords_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.co'
    output_dir = r"D:\Thesis\DIMCAS\NY_original_MC_inflated_heuristics"
    samples = 500
    epsilon = 0.01
    # testbench_E(distance_filename, time_filename, coords_filename, output_dir, samples, path2ApexExe, epsilon)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Analyze testbench_E results
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyze_heuristics_inflation(output_dir, samples, epsilon)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Checking whether heuristics inflation generates eps-bounded solutions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    distance_filename = r"D:\Thesis\DIMCAS\NY_correlated\USA-road-d.NY_rho000.gr"
    # time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    time_filename = r"D:\Thesis\DIMCAS\NY_correlated\USA-road-t.NY_rho000.gr"
    coords_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.co'
    output_dir = r"D:\Thesis\DIMCAS\NY_original_check2"
    samples = 1000
    # testbench_F(distance_filename, time_filename, coords_filename, output_dir, samples, path2ApexExe)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting testbench_E results
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # analyze_testbench_F(output_dir, samples)