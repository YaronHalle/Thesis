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
    analyze_corr_vs_optimality_ratio(output_dir, corr_vec, samples_per_corr)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plotting the graph
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # coords_filename = r'D:\Thesis\DIMCAS\NY_correlated\USA-road-d.NY.co'
    # plot_graph(gr_filename, coords_filename)


