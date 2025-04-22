import numpy as np
import pandas as pd
import scipy.optimize as optimize
import matplotlib.pyplot as plt

if __name__ == "__main__":
    filename = "fig"
    suffix = "png"
    dpi = 300
    if False:
        files = [
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D2.txt", 2, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D2_M1.9.txt", 2,
             1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D2_M2.1.txt", 2,
             2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D3.txt", 3, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D3_M1.9.txt", 3,
             1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D3_M2.1.txt", 3,
             2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5.txt", 5, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M1.9.txt", 5,
             1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M2.1.txt", 5,
             2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D7.txt", 7, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10.txt", 10, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M1.9.txt",
             10, 1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2.1.txt",
             10, 2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2.05.txt",
             10, 2.05),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D12.txt", 12, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15.txt", 15, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D17.txt", 17, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20.txt", 20, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M1.9.txt",
             20, 1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2.1.txt",
             20, 2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2.2.txt",
             20, 2.2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2.15.txt",
             20, 2.15),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D23.txt", 23, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25.txt", 25, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D27.txt", 27, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30.txt", 30, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M1.9.txt",
             30, 1.9),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2.1.txt",
             30, 2.1),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2.2.txt",
             30, 2.2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2.05.txt",
             30, 2.05),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D32.txt", 32, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D35.txt", 35, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D37.txt", 37, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40.txt", 40, 2),
            ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D50.txt", 50, 2),
        ]

        pandas_table = {}
        for file, dimension, metric in files:
            pandas_table[(dimension, metric)] = pd.read_csv(file, delimiter=" ")

        with open("parameter_fitting.txt", "wt") as pf:
            for key, table in pandas_table.items():
                dimension, metric = key
                xdata = table.iloc[:, 0].values
                ydata = table.iloc[:, 21].values

                location = table[table.columns[0]] <= 0.5
                xdata1 = table.loc[location].iloc[:, 0].values
                ydata1 = table.loc[location].iloc[:, 21].values
                function = lambda x, a, b: a * pow(x, b)
                result1, covariance_result = optimize.curve_fit(function, xdata1, ydata1)

                location = table[table.columns[0]] >= 1
                xdata2 = table.loc[location].iloc[:, 0].values
                ydata2 = table.loc[location].iloc[:, 21].values
                result2, covariance_result = optimize.curve_fit(function, xdata2, ydata2)

                plt.plot(xdata1, function(xdata1, *result1), "*", label="fit")
                plt.plot(xdata2, function(xdata2, *result2), "x", label="fit2")
                plt.plot(xdata, ydata, "-", label="Ratio")
                plt.xscale("log")
                plt.yscale("log")
                plt.legend()
                # plt.show()
                plt.savefig(filename + "{}_{}".format(dimension, metric) + "." + suffix, dpi=dpi, bbox_inches="tight")
                print(f"{dimension} {metric} {result1[0]} {result1[1]} {result2[0]} {result2[1]}", file=pf)
                # plt.show()
                plt.close()

    file_set = [
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D5_M2_comp.txt", 5,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D10_M2_comp.txt", 10,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D15_M2_comp.txt", 15,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D20_M2_comp.txt", 20,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D25_M2_comp.txt", 25,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D30_M2_comp.txt", 30,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D40_M2_comp.txt", 40,
         2),
        ("/home/hynek/work/skola/prog/python/transfer_entropy/build_release/src_cpp/renyi_gaussian_D50_M2_comp.txt", 50,
         2),
    ]

    pandas_table = {}
    xdata = []
    ydata = []

    with open("parameter_fitting_absolute.txt", "wt") as pf:
        for file, dimension, metric in file_set:
            pandas_table[(dimension, metric)] = pd.read_csv(file, delimiter=" ")

            for key, table in pandas_table.items():
                dimension, metric = key
                location = table[table.columns[0]] == 0.999
                xdata1 = table.loc[location].iloc[:, 0].values
                ydata1 = table.loc[location].iloc[:, 21].values
                xdata.append(dimension)
                ydata.append(ydata1[0])

                print(f"{dimension} {metric} {xdata1[0]} {ydata1[0]}", file=pf)

        function = lambda x, a, b: a * x + b
        result, covariance_result = optimize.curve_fit(function, xdata, ydata)
        print(f"{result}")
