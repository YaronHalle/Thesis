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
    testbench_B(distance_filename, time_filename, output_dir, samples, path2ApexExe)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Generating a multi-correlated graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\USA-road-d.NY.co'
    distance_filename = r'D:\Thesis\DIMCAS\NY\USA-road-d.NY.gr'
    time_filename = r'D:\Thesis\DIMCAS\NY\USA-road-t.NY.gr'
    # clusters_count = 40
    clusters_count = 0
    new_distance_gr_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\tUSA-road-d.NY_multiple_clusters.gr'
    new_time_gr_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\tUSA-road-t.NY_multiple_clusters.gr'
    new_coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\tUSA-road-d.NY_multiple_clusters.co'
    clusters_metafile = r'D:\Thesis\DIMCAS\NY_correlated\Multiple_Clusters\tNY_multiple_clusters.txt'
    min_lat = 4.058e7
    max_lat = 4.058e7 + (4.08e7 - 4.058e7) * 0.75
    min_lon = -7.41e7
    max_lon = -7.41e7 + (-7.36e7 + 7.41e7) * 0.75
    # min_lat = 4.058e7
    # max_lat = 4.08e7
    # min_lon = -7.41e7
    # max_lon = -7.36e7
    generate_multiple_correlated_graph(distance_filename, time_filename, coords_filename,
                                       clusters_count, new_distance_gr_filename,
                                       new_time_gr_filename, new_coords_filename, clusters_metafile,
                                       min_lon, max_lon, min_lat, max_lat)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting the graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot_graph(gr_filename, coords_filename, clusters_metafile)
    plot_graph(new_distance_gr_filename, new_coords_filename, clusters_metafile)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Graph Clustering
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eps = 0.01
    # graph_agglomerative_clustering(new_distance_gr_filename, new_time_gr_filename, eps)


