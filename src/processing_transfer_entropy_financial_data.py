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

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino"],
    }
)


def load_processed_dataset(
        dataset, dataset_raw, new_columns_base_name="transfer_entropy_"
):
    TE = pd.read_pickle(dataset)
    columns = TE.columns

    TE_raw = pd.read_pickle(dataset_raw)

    return TE, [item for item in TE.columns.tolist() if "mean" in str(item[1])], TE_raw


if __name__ == "__main__":
    dpi = 300
    fontsize = 12
    output = "png"
    latex_title_size = "\\Huge"
    latex_label_size = "\\huge"
    latex_overview_label_size = "\\small"

    directories = [
        # "conditional_information_transfer",
        # "conditional_information_transfer_3",
        # "financial_transfer_entropy"
        # , "financial_transfer_entropy_2", "financial_transfer_entropy_1", "financial_transfer_entropy"
        "indices/short_memory_random=0.001",
    ]

    print(plt.style.available)

    for directory in directories:
        name_of_title = "conditional_information_transfer"
        processed_dataset = directory + "/pivot_dataset.bin"
        processed_raw_dataset = directory + "/pivot_dataset_raw.bin"
        files = glob.glob(processed_dataset)
        if len(files) == 0:
            TE, TE_column_names, TE_raw = process_datasets(
                processed_datasets=directory + "/conditional_information_transfer-*.bin",
                result_dataset=processed_dataset,
                result_raw_dataset=processed_raw_dataset,
                take_k_th_nearest_neighbor=5,
                new_columns_base_name=name_of_title,
                converter_epsilon=lambda x: x.split("-")[1].split(".b")[0],
            )
        else:
            TE, TE_column_names, TE_raw = load_processed_dataset(
                processed_dataset, processed_raw_dataset
            )

        epsilons = TE["epsilon"].unique()
        set_symbols = set()
        for epsilon in epsilons:
            splitted_epsilon = epsilon.split("_")
            set_symbols.add(splitted_epsilon[0])
            set_symbols.add(splitted_epsilon[1])

        TE["symbol_1"] = TE["epsilon"].map(lambda value: value.split("_")[0])
        TE["symbol_2"] = TE["epsilon"].map(lambda value: value.split("_")[1])

        symbols = list(set_symbols)
        print(f"{symbols}")

        TE_nonswapped_columns = [item for item in TE_column_names if item[4] == False]
        TE_cathegories = []
        TE_cathegories_std = []
        for name in names.keys():
            for shuffled in [True, False]:
                TE_cathegory_names = [
                    item for item in TE_column_names if
                    item[0].startswith(name) and item[3] == shuffled and item[4] == False
                ]
                TE_cathegory_names_std = [
                    (item[0], "std", item[2], item[3], item[4]) for item in TE_column_names if
                    item[0].startswith(name) and item[3] == shuffled and item[4] == False
                ]

                if TE_cathegory_names:
                    TE_cathegories.append(TE_cathegory_names)
                    TE_cathegories_std.append(TE_cathegory_names_std)

        for symbol in symbols:
            TE_selected = TE.loc[
                (TE["symbol_1"] == symbol) | (TE["symbol_2"] == symbol)
                ]

            for item, item_std in zip(TE_cathegories, TE_cathegories_std):
                complete_column_name = item[0]
                column_name = item[0][0]
                label = generate_label_from_name_of_column_TE(item[0], name_of_title)

                column_name = column_name.replace(
                    "conditional_information_transfer", "transfer_entropy"
                )
                name_of_title = column_name.split("y_")[0] + "y"
                base_filename = column_name.split("y_")[0] + "y"
                balance = "balance" in name_of_title
                if balance:
                    name_of_title = "Balance of" + name_of_title.split("balance")[1]
                pure_title = name_of_title.capitalize().replace("_", " ") + f" of {symbol}"
                latex_title = latex_title_size + f"""{{{pure_title}}}"""
                latex_alpha_label = latex_overview_label_size + r"$\alpha$"
                shuffled_calculation = complete_column_name[3]
                swapped_datasets = complete_column_name[4]

                standard_filename = (
                        directory
                        + f"/{base_filename}_{symbol}"
                        + ("_shuffled" if shuffled_calculation else "")
                        + "_2d"
                )
                std_filename = (
                        directory
                        + f"/{base_filename}_{symbol}"
                        + ("_shuffled" if shuffled_calculation else "")
                        + "_2d_std"
                )
                errorbar_filename = (
                        directory
                        + f"/{base_filename}_{symbol}"
                        + ("_shuffled" if shuffled_calculation else "")
                        + "_2d_bars"
                )

                figures2d_TE_overview_alpha_errorbar(
                    TE_selected,
                    item,
                    item_std,
                    latex_title,
                    latex_alpha_label,
                    (name_of_title, latex_overview_label_size),
                    errorbar_filename,
                    output,
                    dpi=dpi,
                    fontsize=fontsize,
                )

                figures2d_TE_overview_alpha(
                    TE_selected,
                    item,
                    latex_title,
                    latex_alpha_label,
                    (name_of_title, latex_overview_label_size),
                    standard_filename,
                    output,
                    dpi=dpi,
                    fontsize=fontsize,
                )

                figures2d_TE_overview_alpha(
                    TE_selected,
                    item_std,
                    latex_title,
                    latex_alpha_label,
                    (name_of_title, latex_overview_label_size),
                    std_filename,
                    output,
                    dpi=dpi,
                    fontsize=fontsize,
                )

            for item in TE_nonswapped_columns:
                try:
                    label = generate_label_from_name_of_column_TE(item, name_of_title)

                    complete_column_name = list(item)
                    column_name = item[0]
                    column_name = column_name.replace(
                        "conditional_information_transfer", "transfer_entropy"
                    )

                    shuffled_calculation = complete_column_name[3]
                    swapped_datasets = complete_column_name[4]
                    complete_column_name_std = complete_column_name.copy()
                    complete_column_name_std[1] = "std"

                    name_of_title = column_name.split("y_")[0] + "y"
                    balance = "balance" in name_of_title
                    if balance:
                        name_of_title = "Balance of" + name_of_title.split("balance")[1]
                    pure_title = name_of_title.capitalize().replace("_", " ")
                    latex_title = latex_title_size + f"""{{{pure_title}}}"""
                    latex_title_std = (
                            latex_title_size
                            + f"""{{Standard deviation of {pure_title.lower()} }}"""
                    )

                    title_graph = {
                        "transfer_entropy": r"$\Huge\rm{Transfer\ entropy}$",
                        "conditional_information_transfer": r"$\Huge\rm{Conditional\ information\ transfer}$",
                    }
                    filename_direction = {True: "Y->X", False: "X->Y"}

                    label_latex_std = latex_label_size + f"""$\\sigma_{{{label}}}$"""
                    label = latex_label_size + f"${label}$"
                    latex_alpha_label = latex_label_size + r"$\alpha$"
                    print(column_name, label, label_latex_std)

                    errorbar_filename = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                            + "_2d_bars"
                    )
                    standard_filename = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                            + "_2d"
                    )
                    plot_3D_filename = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                    )
                    plot_3D_surf_filename = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled_surf" if shuffled_calculation else "_surf")
                    )
                    plot_2D_filename_implot = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                            + "_implot"
                    )
                    plot_2D_filename_implot_std = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                            + "_implot_std"
                    )
                    std_filename = (
                            directory
                            + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}"
                            + ("_shuffled" if shuffled_calculation else "")
                            + "_2d_std"
                    )

                    figures2d_TE_alpha(
                        TE_selected,
                        item,
                        latex_title,
                        latex_alpha_label,
                        label,
                        standard_filename,
                        output,
                        dpi=dpi,
                        fontsize=fontsize,
                    )
                    figures2d_TE_alpha(
                        TE_selected,
                        tuple(complete_column_name_std),
                        latex_title_std,
                        latex_alpha_label,
                        label_latex_std,
                        std_filename,
                        output,
                        dpi=dpi,
                        fontsize=fontsize,
                    )
                    figures2d_TE_alpha_errorbar(
                        TE_selected,
                        item,
                        tuple(complete_column_name_std),
                        latex_title,
                        latex_alpha_label,
                        label,
                        errorbar_filename,
                        output,
                        dpi=dpi,
                        fontsize=fontsize,
                    )
                except Exception as exc:
                    print(f"Problem {exc} {item}")
                    traceback.print_exc()

        del TE, TE_column_names, TE_raw

    print("Finished")
