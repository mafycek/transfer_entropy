import numpy as np
import pandas as pd
import scipy.optimize as optimize
import matplotlib.pyplot as plt

if __name__ == "__main__":
    suffix = "png"
    dpi = 300

    if False:
        files = [
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D3.txt", 3, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5.txt", 5, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10.txt", 10, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15.txt", 15, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20.txt", 20, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25.txt", 25, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30.txt", 30, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D35.txt", 35, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40.txt", 40, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40.txt", 45, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40.txt", 50, 2),
        ]

    if True:
        files = [
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M2_comp.txt", 5, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D7_M2_comp.txt", 7, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2_comp.txt", 10, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D12_M2_comp.txt", 12, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15_M2_comp.txt", 15, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2_comp.txt", 20, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25_M2_comp.txt", 25, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2_comp.txt", 30, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40_M2_comp.txt", 40, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D50_M2_comp.txt", 50, 2),
        ]

    if False:
        files = [
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M2_correlated.txt",
             5, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D7_M2_correlated.txt",
             7, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2_correlated.txt",
             10, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D12_M2_correlated.txt",
             12, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15_M2_correlated.txt",
             15, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D22_M2_correlated.txt",
             22, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25_M2_correlated.txt",
             25, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D27_M2_correlated.txt",
             27, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2_correlated.txt",
             30, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D32_M2_correlated.txt",
             32, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D35_M2_correlated.txt",
             35, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D37_M2_correlated.txt",
             37, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40_M2_correlated.txt",
             40, 2),
        ]

    if False:
        files = [
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M2_full_correlated.txt",
             5, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D7_M2_full_correlated.txt",
             7, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2_full_correlated.txt",
             10, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D12_M2_full_correlated.txt",
             12, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15_M2_full_correlated.txt",
             15, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D22_M2_full_correlated.txt",
             22, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2_full_correlated.txt",
             20, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25_M2_full_correlated.txt",
             25, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D27_M2_full_correlated.txt",
             27, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2_full_correlated.txt",
             30, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D32_M2_full_correlated.txt",
             32, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D35_M2_full_correlated.txt",
             35, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D37_M2_full_correlated.txt",
             37, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40_M2_full_correlated.txt",
             40, 2),
        ]

    limit_of_calculated_values = 10
    filename = "renyi_entropy_LP_normal"
    pandas_table = {}
    for file, dimension, metric in files:
        pandas_table[(dimension, metric)] = pd.read_csv(file, delimiter=" ")

    for key, table in pandas_table.items():
        dimension, metric = key
        xdata = table.iloc[:, 0].values
        plt.title(f"Renyi entropy k=10")
        for k in range(1, limit_of_calculated_values):
            ydata = table.iloc[:, k].values
            plt.plot(xdata, ydata, "-", label=f"k={k}")

        ydata = table.iloc[:, limit_of_calculated_values+1].values
        plt.plot(xdata, ydata, ".-", label="Theoretical value")
        plt.legend(ncol=3)
        plt.savefig(filename + "_{}_{}".format(dimension, metric) + "." + suffix, dpi=dpi, bbox_inches="tight")
        plt.close()

    filename = "relative_renyi_entropy_LP_normal"
    for key, table in pandas_table.items():
        dimension, metric = key
        xdata = table.iloc[:, 0].values
        plt.title(f"Renyi entropy k=10")
        for k in range(limit_of_calculated_values+2, 2*limit_of_calculated_values+2):
            ydata = table.iloc[:, k].values
            plt.plot(xdata, ydata, "-", label=f"k={k}")

        plt.legend(ncol=3)
        plt.savefig(filename + "_{}_{}".format(dimension, metric) + "." + suffix, dpi=dpi, bbox_inches="tight")
        plt.close()

    for key, table in pandas_table.items():
        dimension, metric = key
        ydata = table.iloc[:, limit_of_calculated_values+1].values
        xdata = table.iloc[:, limit_of_calculated_values].values

        plt.plot(xdata, ydata, label=f"D={dimension}")
    plt.legend(ncol=3)
    plt.savefig(filename + "_comparison".format(dimension, metric) + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
