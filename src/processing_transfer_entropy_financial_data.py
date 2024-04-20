#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import traceback
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from processing_datasets import *
from figure_generation import *

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino"],
    }
)


def load_processed_dataset(
        directory, dataset, dataset_raw, new_columns_base_name="transfer_entropy_"
):
    TE = pd.read_pickle(dataset)
    columns = TE.columns

    TE_raw = pd.read_pickle(dataset_raw)
    TE_statistics = pd.read_pickle(f"{directory}/refined_statistics.bin")

    return (
        TE,
        [item for item in TE.columns.tolist() if "mean" in str(item[1])],
        TE_raw,
        TE_statistics,
    )


swap_column_number = 5
shuffled_column_number = 4

if __name__ == "__main__":
    directories = [
        # "conditional_information_transfer",
        # "conditional_information_transfer_3",
        # "financial_transfer_entropy"
        # , "financial_transfer_entropy_2", "financial_transfer_entropy_1", "financial_transfer_entropy"
        # "indices/short_memory_random=0.1",
        "test"
    ]

    print(f"Available plot styles {plt.style.available}")

    for directory in directories:
        figure_parameters = {
            "dpi": 300,
            "fontsize": 12,
            "output_format": "png",
            "latex_title_size": "\\Huge",
            "latex_label_size": "\\huge",
            "latex_overview_label_size": "\\small",
            "directory": directory,
        }

        name_of_title = "conditional_information_transfer"
        processed_dataset = directory + "/pivot_dataset.bin"
        processed_raw_dataset = directory + "/pivot_dataset_raw.bin"
        files = glob.glob(processed_dataset)
        if len(files) == 0:
            TE, TE_column_names, TE_raw, TE_statistics = process_datasets(
                directory=directory,
                processed_datasets=directory
                                   + "/conditional_information_transfer-*.bin",
                result_dataset=processed_dataset,
                result_raw_dataset=processed_raw_dataset,
                take_k_th_nearest_neighbor=5,
                new_columns_base_name=name_of_title,
                converter_epsilon=lambda x: x.split("-")[1].split(".b")[0],
            )
        else:
            TE, TE_column_names, TE_raw, TE_statistics = load_processed_dataset(
                directory, processed_dataset, processed_raw_dataset
            )

        symbols = generate_epsilon_values(TE, "epsilon", "symbol_1", "symbol_2")
        print(f"{symbols}")

        TE_nonswapped_columns = [
            item for item in TE_column_names if item[swap_column_number] == False
        ]
        TE_cathegories = []
        TE_cathegories_std = []
        for name in names.keys():
            for shuffled in [True, False]:
                TE_cathegory_names = [
                    item
                    for item in TE_column_names
                    if item[0].startswith(name)
                       and item[shuffled_column_number] == shuffled
                       and item[swap_column_number] == False
                ]
                TE_cathegory_names_std = [
                    (item[0], "std", item[2], item[3], item[4], item[5], item[6])
                    for item in TE_column_names
                    if item[0].startswith(name)
                       and item[shuffled_column_number] == shuffled
                       and item[swap_column_number] == False
                ]

                if TE_cathegory_names:
                    TE_cathegories.append(TE_cathegory_names)
                    TE_cathegories_std.append(TE_cathegory_names_std)

        # generate_figures_from_data(TE, TE_cathegories, TE_cathegories_std, TE_nonswapped_columns, symbols, directory, name_of_title, figure_parameters)
        # generate_figures_from_neighborhood_statistics(TE_statistics, figure_parameters)
        generate_overview_figures_from_averages(TE_statistics, figure_parameters)

    print("Finished")
