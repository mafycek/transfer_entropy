#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

from figure_generation import *

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino"],
    }
)

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
        merged_dataset = f"{directory}/merged_{name_of_title}.bin"

        files = glob.glob(merged_dataset)
        if len(files) == 0:
            # create dataset
            pass
        else:
            table = pd.read_pickle(merged_dataset)
            # generate_overview_figures_from_averages(table, figure_parameters)
            generate_figures_from_neighborhood_statistics(table, figure_parameters)
