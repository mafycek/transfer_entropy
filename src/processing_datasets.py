#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import sys
import math
import datetime
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

title_map = {
    (False, False): r"{\alpha: X\rightarrow Y}",
    (True, False): r"{\alpha: X_{shuffled}\rightarrow Y}",
    (False, True): r"{\alpha: Y\rightarrow X}",
    (True, True): r"{\alpha: Y_{shuffled}\rightarrow X}",
}

names = {
    "balance_effective_transfer_entropy": 5,
    "balance_transfer_entropy": 4,
    "effective_transfer_entropy": 4,
    "transfer_entropy": 3,
}


def max_divisor(number):
    limit = math.ceil(math.sqrt(number))
    for divisor in range(limit, 1, -1):
        if number % divisor == 0:
            return (divisor, number // divisor)
    return (1, number)


def generate_label_from_name_of_column_TE(column_multiindex, name_of_title):
    swap_column_number = 5
    shuffled_column_number = 4

    column_name = column_multiindex[0]
    shift = 0
    for key, value in names.items():
        if key in column_name:
            shift = value
            break

    column_name = column_name.replace(
        "conditional_information_transfer", "transfer_entropy"
    )

    shift = shift - 1
    try:
        future_first_TS = column_name.split("_")[shift + 2]
    except IndexError as err:
        future_first_TS = None
    balance = "balance" in name_of_title

    shuffled_calculation = column_multiindex[shuffled_column_number]
    swapped_datasets = column_multiindex[swap_column_number]
    history_first_TS = column_name.split("_")[shift]
    history_second_TS = column_name.split("_")[shift + 1]
    standard_filename = "figure"

    if future_first_TS is not None:
        if balance:
            label = """T^{}_{} (\\{{{}\\}},\\{{{}\\}},\\{{{}\\}})""".format(
                (
                    "{(R, effective, balance)}"
                    if "effective" in column_name
                    else "{(R, balance)}"
                ),
                title_map[(shuffled_calculation, swapped_datasets)],
                history_first_TS,
                history_second_TS,
                future_first_TS,
            )
            # + "-" + "$T^{}_{} ([{}],[{}],[{}])$".format("{(R, eff)}" if "effective" in column_name else "{(R)}", title_map[(shuffled_calculation, not swapped_datasets)], history_first_TS, history_second_TS, future_first_TS)
        else:
            label = """T^{}_{} (\\{{{}\\}},\\{{{}\\}},\\{{{}\\}})""".format(
                "{(R, eff)}" if "effective" in column_name else "{(R)}",
                title_map[(shuffled_calculation, swapped_datasets)],
                history_first_TS,
                history_second_TS,
                future_first_TS,
            )
    else:
        label = "T^{}_{} ([{}],[{}])".format(
            "{(R, eff)}" if "effective" in column_name else "{(R)}",
            title_map[(shuffled_calculation, swapped_datasets)],
            history_first_TS,
            history_second_TS,
        )
    return label


def translate_to_label_TE(column_name):
    (
        timeseries_name,
        history_first_TS,
        future_first_TS,
        history_second_TS,
        swapped_datasets,
        shuffled_calculation,
    ) = column_name
    balance = "balance" in timeseries_name
    effective = "effective" in timeseries_name
    standard_filename = "figure"

    type_of_TE = None
    if effective and balance:
        type_of_TE = "{(R, effective, balance)}"
    elif not effective and balance:
        type_of_TE = "{(R, balance)}"
    elif effective and not balance:
        type_of_TE = "{(R, effective)}"
    else:
        type_of_TE = "{(R)}"

    return """T^{}_{} (\\{{{}\\}},\\{{{}\\}},\\{{{}\\}})""".format(
        type_of_TE,
        title_map[(shuffled_calculation, swapped_datasets)],
        history_first_TS,
        history_second_TS,
        future_first_TS,
    )


def figures3d_TE(
        dataset,
        selector,
        title,
        xlabel,
        ylabel,
        zlabel,
        filename,
        suffix,
        view=(50, -20),
        dpi=300,
):
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    colors = ["r", "g", "b", "c", "m", "y", "k", "orange", "pink"]
    markers = ["b", "^"]

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    row_size = len(dataset["epsilon"].unique())
    xs = dataset[["alpha"]]
    ys = dataset[["epsilon"]]
    zs = dataset[[selector]]

    try:
        ax.plot_wireframe(
            np.reshape(xs.values, (-1, row_size)),
            np.reshape(ys.values, (-1, row_size)),
            np.reshape(zs.values, (-1, row_size)),
            rstride=1,
            cstride=1,
            color=colors[0],
            linewidth=1,
        )
    except Exception as exc:
        print(f"{exc}: Problem D=")

    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.legend(loc=1)
    ax.view_init(view[0], view[1])

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    # plt.draw()
    # plt.show()
    plt.close()
    del fig


def figures3d_surface_TE(
        dataset,
        selector,
        title,
        xlabel,
        ylabel,
        zlabel,
        filename,
        suffix,
        cmap="magma",
        view=(50, -20),
        dpi=300,
):
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    colors = ["r", "g", "b", "c", "m", "y", "k", "orange", "pink"]
    markers = ["b", "^"]

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    row_size = len(dataset["epsilon"].unique())
    xs = dataset[["alpha"]]
    ys = dataset[["epsilon"]]
    zs = dataset[[selector]]

    try:
        ax.plot_surface(
            np.reshape(xs.values, (-1, row_size)),
            np.reshape(ys.values, (-1, row_size)),
            np.reshape(zs.values, (-1, row_size)),
            rstride=1,
            cstride=1,
            cmap=cmap,
            linewidth=0,
            antialiased=False,
        )
    except Exception as exc:
        print(f"{exc}: Problem D=")

    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.legend(loc=1)
    ax.view_init(view[0], view[1])

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    # plt.draw()
    # plt.show()
    plt.close()
    del fig


def minimal_difference(target):
    epsilon_differences = []
    for item in range(0, len(target) - 1):
        epsilon_differences.append(round(target[item + 1] - target[item], 4))
    return min(epsilon_differences)


def figures2d_imshow(
        dataset, selector, title, xlabel, ylabel, filename, suffix, cmap="magma", dpi=300
):
    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig, ax = plt.subplots(1, 1, figsize=(13, 8))

    ax.set(title=title)
    ax.set(xlabel=xlabel)
    ax.set(ylabel=ylabel)
    ax.grid(True)
    # ax.set_ylim([0.85, 1.5])

    # dataset = dataset.loc[ dataset["alpha"].between(0.85, 1.5) ]

    epsilons = dataset["epsilon"].unique()
    alphas = dataset["alpha"].unique()
    epsilon = dataset[["epsilon"]]
    alpha = dataset[["alpha"]]
    data = dataset[[("epsilon", "", "", "", ""), ("alpha", "", "", "", ""), selector]]
    xs = dataset[["epsilon"]].values.reshape((len(alphas), len(epsilons)))
    ys = dataset[["alpha"]].values.reshape((len(alphas), len(epsilons)))
    zs = dataset[[selector]].values.reshape((len(alphas), len(epsilons)))
    coords = np.array(list(zip(xs.flatten(), ys.flatten())))

    minimal_epsilon_difference = minimal_difference(epsilons)
    changed_epsilons = np.arange(epsilons[0], epsilons[-1], minimal_epsilon_difference)
    flatten_zs = zs.flatten()

    grid = np.dstack(np.meshgrid(changed_epsilons, alphas)).reshape(-1, 2)
    sampled_data = griddata(coords, flatten_zs, grid, method="nearest")

    number_epsilons = len(epsilons)
    number_alphas = len(alphas)
    epsilon_margin = (max(epsilons) - min(epsilons)) / (number_epsilons * 2.0)
    alpha_margin = (max(alphas) - min(alphas)) / (number_alphas * 2.0)
    extent = [
        min(epsilons) - epsilon_margin,
        max(epsilons) + epsilon_margin,
        min(alphas) - alpha_margin,
        max(alphas) + alpha_margin,
    ]
    ims = ax.imshow(
        sampled_data.reshape((len(alphas), len(changed_epsilons))),
        origin="lower",
        interpolation="nearest",
        extent=extent,
        cmap=color_map,
        aspect="auto",
    )

    fig.colorbar(ims)
    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    # plt.show()
    plt.close()
    del fig


def figures2d_TE_alpha(
        dataset,
        selector,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1)

    color_map = matplotlib.colormaps.get_cmap(cmap)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    codes = dataset["epsilon"].unique()
    list_selector = list(selector)
    list_selector[4] = not selector[4]
    selector_not = tuple(list_selector)
    columns = list(dataset.columns.values)
    if selector_not in columns:
        number_of_datasets = float(2 * len(codes)) - 1
    else:
        number_of_datasets = float(len(codes)) - 1

    order_of_dataset = 0
    for code in codes:
        subselection = dataset.loc[dataset["epsilon"] == code]

        for swap in [True, False]:
            list_selector = list(selector)
            list_selector[4] = swap
            selector = tuple(list_selector)
            if selector in columns:
                label = code.replace("_", "-")
                if swap:
                    label_split = label.split("-")
                    if len(label_split) >= 2:
                        label = f"${{ \\scriptstyle {label_split[1]}-{label_split[0]}, Y->X }}$"
                    else:
                        label = f"${label}, Y->X$"
                else:
                    label = f"${{ \\scriptstyle {label}, X->Y }}$"

                ys = subselection[["alpha"]]
                zs = subselection[[selector]]

                try:
                    if number_of_datasets != 0.0:
                        map_position = order_of_dataset / number_of_datasets
                    else:
                        map_position = 0
                    color = color_map(map_position)
                    ax.plot(ys.values, zs.values, linewidth=3, label=label, color=color)
                except Exception as exc:
                    print(f"{exc}: Problem D=")

                order_of_dataset += 1

    plt.legend(loc=0, ncol=2, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_overview_alpha(
        dataset,
        selectors,
        title,
        xlabel,
        y_label_params,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(len(selectors))
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16))

    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig.suptitle(title)
    name_of_title, latex_overview_label_size = y_label_params

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_selector = column + number_columns * row
            selector = selectors[index_of_selector]

            ylabel = generate_label_from_name_of_column_TE(selector, name_of_title)
            ylabel_latex = latex_overview_label_size + f"${ylabel}$"

            if number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel_latex)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel_latex)
            else:
                axs[row, column].set_xlabel(xlabel)
                axs[row, column].set_ylabel(ylabel_latex)

            codes = dataset["epsilon"].unique()
            list_selector = list(selector)
            list_selector[4] = not selector[4]
            selector_not = tuple(list_selector)
            columns = list(dataset.columns.values)
            if selector_not in columns:
                number_of_datasets = float(2 * len(codes)) - 1
            else:
                number_of_datasets = float(len(codes)) - 1

            order_of_dataset = 0
            for code in codes:
                subselection = dataset.loc[dataset["epsilon"] == code]

                for swap in [True, False]:
                    list_selector = list(selector)
                    list_selector[4] = swap
                    selector = tuple(list_selector)
                    if selector in columns:
                        label = code.replace("_", "-")
                        if swap:
                            label_split = label.split("-")
                            if len(label_split) >= 2:
                                label = f"${{ \\scriptstyle {label_split[1]}-{label_split[0]}, Y->X }}$"
                            else:
                                label = f"${label}, Y->X$"
                        else:
                            label = f"${{ \\scriptstyle {label}, X->Y }}$"

                        ys = subselection[["alpha"]]
                        zs = subselection[[selector]]

                        try:
                            if number_of_datasets != 0.0:
                                map_position = order_of_dataset / number_of_datasets
                            else:
                                map_position = 0
                            color = color_map(map_position)
                            if number_rows == 1:
                                axs[column].plot(
                                    ys.values,
                                    zs.values,
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                )
                            elif number_columns == 1:
                                axs[row].plot(
                                    ys.values,
                                    zs.values,
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                )
                            else:
                                axs[row, column].plot(
                                    ys.values,
                                    zs.values,
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                )
                        except Exception as exc:
                            print(f"{exc}: Problem D=")

                        order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_overview_alpha_order_statistics(
        dataset,
        number_selectors,
        selectors,
        title,
        xlabel,
        y_label_params,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(number_selectors)
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16))

    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig.suptitle(title)
    name_of_title, latex_overview_label_size = y_label_params

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_dataset = column + number_columns * row

            ylabel = (
                "AA"  # + generate_label_from_name_of_column_TE(selector, name_of_title)
            )
            ylabel_latex = latex_overview_label_size + f"${ylabel}$"

            if number_rows == 1 and number_columns == 1:
                axs.set_xlabel(xlabel)
                axs.set_ylabel(ylabel_latex)
            elif number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel_latex)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel_latex)
            else:
                axs[row, column].set_xlabel(xlabel)
                if (row == number_rows // 2) and (column == 0):
                    axs[row, column].set_ylabel(ylabel_latex)

            order_of_dataset = 0
            for selector_item in selectors:
                try:
                    # subselection = [dataset[item] for item in selector_item["selector"]]
                    subselection = [
                        dataset.xs(item, level="Statistical value", axis=1)
                        for item in selector_item["selector"]
                    ]

                    label = selector_item["label"]
                    ys = subselection[0].index.to_list()
                    zs = []
                    for index in range(len(subselection)):
                        # for item in ys:
                        #    print(subselection[index].loc[item], subselection[index].loc[item].values[0][index_of_order], len(subselection[index].loc[item].values.shape))
                        if selector_item["aggregation"] == "sample":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index][index_of_dataset].loc[key]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )
                        elif selector_item["aggregation"] == "neighborhood":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index]
                                        .loc[key]
                                        .values[0][index_of_dataset]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )  # [index_of_order]

                    color = selector_item["color"]
                    actual_axes = None
                    if number_rows == 1 and number_columns == 1:
                        actual_axes = axs
                    elif number_rows == 1:
                        actual_axes = axs[column]
                    elif number_columns == 1:
                        actual_axes = axs[row]
                    else:
                        actual_axes = axs[row, column]

                    if (
                            len(subselection) == 2
                            and selector_item["style"] == "quantile band"
                    ):
                        actual_axes.fill_between(
                            ys, zs[0], zs[1], linewidth=1, label=label, color=color
                        )
                    elif len(subselection) == 1 and selector_item["style"] == "line":
                        actual_axes.plot(
                            ys, zs[0], linewidth=1, label=label, color=color
                        )
                    elif (
                            len(subselection) == 2 and selector_item["style"] == "yerrorbar"
                    ):
                        lower_band = zs[0] - zs[1]
                        upper_band = zs[0] + zs[1]
                        actual_axes.fill_between(
                            ys,
                            lower_band,
                            upper_band,
                            linewidth=1,
                            label=label,
                            color=color,
                        )
                    else:
                        print("problem")

                except Exception as exc:
                    print(traceback.format_exc(), selector_item, selectors)

                order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_overview_alpha_timeseries(
        dataset,
        timeseries_groups,
        selectors,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    number_selectors = len(timeseries_groups)
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(number_selectors)
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16))
    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig.suptitle(title)

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_dataset = column + number_columns * row
            subselection_timeseries = dataset.xs(
                timeseries_groups[index_of_dataset], level="Sample", axis=1
            )

            if number_rows == 1 and number_columns == 1:
                axs.set_xlabel(xlabel)
                axs.set_ylabel(ylabel)
            elif number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel)
            else:
                axs[row, column].set_title(f"{timeseries_groups[index_of_dataset]}")
                axs[row, column].set_xlabel(xlabel)
                if column == 0:
                    axs[row, column].set_ylabel(ylabel)

            order_of_dataset = 0
            for selector_item in selectors:
                try:
                    # subselection = [dataset[item] for item in selector_item["selector"]]
                    subselection = [
                        subselection_timeseries.xs(
                            item, level="Statistical value", axis=1
                        )
                        for item in selector_item["selector"]
                    ]

                    label = selector_item["label"]
                    ys = subselection[0].index.to_list()
                    zs = []
                    for index in range(len(subselection)):
                        # for item in ys:
                        #    print(subselection[index].loc[item], subselection[index].loc[item].values[0][index_of_order], len(subselection[index].loc[item].values.shape))
                        if selector_item["aggregation"] == "sample":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index][0].loc[key]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )
                        elif selector_item["aggregation"] == "neighborhood":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index].loc[key].values[0][0]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )  # [index_of_order]

                    # prepare visualization
                    color = selector_item["color"]
                    actual_axes = None
                    if number_rows == 1 and number_columns == 1:
                        actual_axes = axs
                    elif number_rows == 1:
                        actual_axes = axs[column]
                    elif number_columns == 1:
                        actual_axes = axs[row]
                    else:
                        actual_axes = axs[row, column]

                    if (
                            len(subselection) == 2
                            and selector_item["style"] == "quantile band"
                    ):
                        actual_axes.fill_between(
                            ys, zs[0], zs[1], linewidth=1, label=label, color=color
                        )
                    elif len(subselection) == 1 and selector_item["style"] == "line":
                        actual_axes.plot(
                            ys, zs[0], linewidth=1, label=label, color=color
                        )
                    elif (
                            len(subselection) == 2 and selector_item["style"] == "yerrorbar"
                    ):
                        lower_band = zs[0] - zs[1]
                        upper_band = zs[0] + zs[1]
                        actual_axes.fill_between(
                            ys,
                            lower_band,
                            upper_band,
                            linewidth=1,
                            label=label,
                            color=color,
                        )
                    else:
                        print("problem")

                except Exception as exc:
                    print(traceback.format_exc(), selector_item, selectors)

                order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_overview_various_parameters(
        dataset,
        timeseries_groups,
        selectors,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    number_selectors = len(timeseries_groups)
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(number_selectors)
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16), sharex=True, sharey=True)
    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig.suptitle(title)

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_dataset = column + number_columns * row
            subselection_timeseries = dataset.xs(
                timeseries_groups[index_of_dataset], level="Randomness", axis=1
            )  # level="Sample"

            if number_rows == 1 and number_columns == 1:
                axs.set_xlabel(xlabel)
                axs.set_ylabel(ylabel)
            elif number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel)
            else:
                axs[row, column].set_title(f"${timeseries_groups[index_of_dataset]}$")
                axs[row, column].set_xlabel(xlabel)
                if column == 0:
                    axs[row, column].set_ylabel(ylabel)

            order_of_dataset = 0
            for selelectior_order in enumerate(selectors["selector"]):
                try:
                    # print(selectors["selector"], subselection_timeseries, )
                    subselection = subselection_timeseries[selelectior_order[1]]
                    # .xs(
                    #    tuple(selectors["selector"]), level="Sample", axis=1
                    # )

                    label = selectors["label"]
                    # ys = subselection[0].index.to_list()
                    ys = subselection.index.to_list()
                    zs = []

                    #for index in selectors["selector"]:
                        # for item in ys:
                        #    print(subselection[index].loc[item], subselection[index].loc[item].values[0][index_of_order], len(subselection[index].loc[item].values.shape))
                    if selectors["aggregation"] == "sample":
                        zs.append(
                            np.array(
                                [
                                    subselection.loc[key]
                                    for index_key_selection, key in enumerate(ys)
                                ]
                            )
                        )
                    elif selectors["aggregation"] == "neighborhood":
                        zs.append(
                            np.array(
                                [
                                    subselection[index].loc[key].values[0][0]
                                    for index_key_selection, key in enumerate(ys)
                                ]
                            )
                        )  # [index_of_order]

                    # prepare visualization
                    color = selectors["color"][selelectior_order[0]]
                    actual_axes = None
                    if number_rows == 1 and number_columns == 1:
                        actual_axes = axs
                    elif number_rows == 1:
                        actual_axes = axs[column]
                    elif number_columns == 1:
                        actual_axes = axs[row]
                    else:
                        actual_axes = axs[row, column]

                    if len(zs) == 1 and selectors["style"] == "line":
                        actual_axes.plot(
                            ys, zs[0], linewidth=1, label=label, color=color
                        )
                    else:
                        print("problem")

                except Exception as exc:
                    print(traceback.format_exc(), selectors)

                order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_complet_overview_alpha_timeseries(
        dataset,
        timeseries_groups,
        selectors,
        title,
        xlabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=10,
        style="seaborn-v0_8",
):
    number_selectors = len(timeseries_groups)
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(number_selectors)
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16))
    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig.suptitle(title)

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_dataset = column + number_columns * row

            ylabel = timeseries_groups[index_of_dataset][1]
            subselection_timeseries = dataset.xs(
                timeseries_groups[index_of_dataset][0], level="Sample", axis=1
            )

            if number_rows == 1 and number_columns == 1:
                axs.set_xlabel(xlabel)
                axs.set_ylabel(ylabel)
            elif number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel)
            else:
                axs[row, column].set_title(f"{timeseries_groups[index_of_dataset]}")
                axs[row, column].set_xlabel(xlabel)
                if column == 0:
                    axs[row, column].set_ylabel(ylabel)

            order_of_dataset = 0
            for selector_item in selectors:
                try:
                    # subselection = [dataset[item] for item in selector_item["selector"]]
                    subselection = [
                        subselection_timeseries.xs(
                            item, level="Statistical value", axis=1
                        )
                        for item in selector_item["selector"]
                    ]

                    label = selector_item["label"]
                    ys = subselection[0].index.to_list()
                    zs = []
                    for index in range(len(subselection)):
                        # for item in ys:
                        #    print(subselection[index].loc[item], subselection[index].loc[item].values[0][index_of_order], len(subselection[index].loc[item].values.shape))
                        if selector_item["aggregation"] == "sample":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index][0].loc[key]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )
                        elif selector_item["aggregation"] == "neighborhood":
                            zs.append(
                                np.array(
                                    [
                                        subselection[index].loc[key].values[0][0]
                                        for index_key_selection, key in enumerate(ys)
                                    ]
                                )
                            )  # [index_of_order]

                    # prepare visualization
                    color = selector_item["color"]
                    actual_axes = None
                    if number_rows == 1 and number_columns == 1:
                        actual_axes = axs
                    elif number_rows == 1:
                        actual_axes = axs[column]
                    elif number_columns == 1:
                        actual_axes = axs[row]
                    else:
                        actual_axes = axs[row, column]

                    if (
                            len(subselection) == 2
                            and selector_item["style"] == "quantile band"
                    ):
                        actual_axes.fill_between(
                            ys, zs[0], zs[1], linewidth=1, label=label, color=color
                        )
                    elif len(subselection) == 1 and selector_item["style"] == "line":
                        actual_axes.plot(
                            ys, zs[0], linewidth=1, label=label, color=color
                        )
                    elif (
                            len(subselection) == 2 and selector_item["style"] == "yerrorbar"
                    ):
                        lower_band = zs[0] - zs[1]
                        upper_band = zs[0] + zs[1]
                        actual_axes.fill_between(
                            ys,
                            lower_band,
                            upper_band,
                            linewidth=1,
                            label=label,
                            color=color,
                        )
                    else:
                        print("problem")

                except Exception as exc:
                    print(traceback.format_exc(), selector_item, selectors)

                order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_overview_alpha_errorbar(
        dataset,
        selectors,
        error_selectors,
        title,
        xlabel,
        y_label_params,
        filename,
        suffix,
        view=(70, 120),
        cmap="rainbow",
        dpi=300,
        fontsize=15,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)
    number_rows, number_columns = max_divisor(len(selectors))
    fig, axs = plt.subplots(number_rows, number_columns, figsize=(26, 16))

    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig.suptitle(title)
    name_of_title, latex_overview_label_size = y_label_params

    for row in range(number_rows):
        for column in range(number_columns):
            index_of_selector = column + number_columns * row
            selector = selectors[index_of_selector]
            error_selector = error_selectors[index_of_selector]

            ylabel = generate_label_from_name_of_column_TE(selector, name_of_title)
            ylabel_latex = latex_overview_label_size + f"${ylabel}$"

            if number_rows == 1:
                axs[column].set_xlabel(xlabel)
                axs[column].set_ylabel(ylabel_latex)
            elif number_columns == 1:
                axs[row].set_xlabel(xlabel)
                axs[row].set_ylabel(ylabel_latex)
            else:
                axs[row, column].set_xlabel(xlabel)
                axs[row, column].set_ylabel(ylabel_latex)

            codes = dataset["epsilon"].unique()
            list_selector = list(selector)
            list_selector[4] = not selector[4]
            selector_not = tuple(list_selector)
            columns = list(dataset.columns.values)
            if selector_not in columns:
                number_of_datasets = float(2 * len(codes)) - 1
            else:
                number_of_datasets = float(2 * len(codes)) - 1

            order_of_dataset = 0
            for code in codes:
                subselection = dataset.loc[dataset["epsilon"] == code]

                for swap in [True, False]:
                    list_selector = list(selector)
                    list_selector[4] = swap
                    selector = tuple(list_selector)
                    if selector in columns:
                        label = code.replace("_", "-")
                        if swap:
                            label_split = label.split("-")
                            if len(label_split) >= 2:
                                label = f"${label_split[1]}-{label_split[0]}, Y->X$"
                            else:
                                label = f"${label}, Y->X$"
                        else:
                            label = f"${label}, X->Y$"

                        ys = subselection[["alpha"]]
                        zs = subselection[[selector]]
                        error_bar = subselection[[error_selector]].copy()

                        error_selector_negative_std = list(error_selector)
                        error_selector_negative_std[1] = "-std"
                        # error_bar[tuple(error_selector_negative_std)] = error_bar.apply(lambda x: -x, axis=1, raw=True)
                        errors = error_bar.copy().T.to_numpy()

                        try:
                            map_position = order_of_dataset / number_of_datasets
                            color = color_map(map_position)
                            lims = np.array([True] * ys.size, dtype=bool)
                            if number_rows == 1:
                                axs[column].errorbar(
                                    ys.values.flatten(),
                                    zs.values.flatten(),
                                    yerr=errors.flatten(),
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                    ls="dotted",
                                )
                            elif number_columns == 1:
                                axs[row].errorbar(
                                    ys.values.flatten(),
                                    zs.values.flatten(),
                                    yerr=errors.flatten(),
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                    ls="dotted",
                                )
                            else:
                                axs[row, column].errorbar(
                                    ys.values.flatten(),
                                    zs.values.flatten(),
                                    yerr=errors.flatten(),
                                    linewidth=3,
                                    label=label,
                                    color=color,
                                    ls="dotted",
                                )
                        except Exception as exc:
                            print(f"{exc}: {errors.shape}")
                        order_of_dataset += 1

            plt.legend(loc=0, ncol=2, fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_alpha_errorbar(
        dataset,
        selector,
        error_selector,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        view=(70, 120),
        cmap="rainbow",
        dpi=300,
        fontsize=15,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.cm.get_cmap(cmap)

    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1)

    markers = ["b", "^"]

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    codes = dataset["epsilon"].unique()
    list_selector = list(selector)
    list_selector[4] = not selector[4]
    selector_not = tuple(list_selector)
    columns = list(dataset.columns.values)
    if selector_not in columns:
        number_of_datasets = float(2 * len(codes)) - 1
    else:
        number_of_datasets = float(len(codes)) - 1

    order_of_dataset = 0
    for code in codes:
        subselection = dataset.loc[dataset["epsilon"] == code]

        for swap in [True, False]:
            list_selector = list(selector)
            list_selector[4] = swap
            selector = tuple(list_selector)
            if selector in columns:
                label = code.replace("_", "-")
                if swap:
                    label_split = label.split("-")
                    if len(label_split) >= 2:
                        label = f"${label_split[1]}-{label_split[0]}, Y->X$"
                    else:
                        label = f"${label}, Y->X$"
                else:
                    label = f"${label}, X->Y$"

                ys = subselection[["alpha"]]
                zs = subselection[[selector]]
                error_bar = subselection[[error_selector]].copy()

                error_selector_negative_std = list(error_selector)
                error_selector_negative_std[1] = "-std"
                # error_bar[tuple(error_selector_negative_std)] = error_bar.apply(lambda x: -x, axis=1, raw=True)
                errors = error_bar.copy().T.to_numpy()

                try:
                    map_position = order_of_dataset / number_of_datasets
                    color = color_map(map_position)
                    lims = np.array([True] * ys.size, dtype=bool)
                    ax.errorbar(
                        ys.values.flatten(),
                        zs.values.flatten(),
                        yerr=errors.flatten(),
                        linewidth=3,
                        label=label,
                        color=color,
                        ls="dotted",
                    )
                except Exception as exc:
                    print(f"{exc}: {errors.shape}")
                order_of_dataset += 1

    plt.legend(loc=0, ncol=2, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE(
        dataset,
        selector,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=15,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1)

    markers = ["b", "^"]

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.set_ylim([-0.1, 1])

    alphas = dataset["alpha"].unique()
    mean = int(len(alphas) / 2)
    neghborhood = 5
    # subselected_alphas = alphas[mean - neghborhood:  mean + neghborhood]
    subselected_alphas = [
        alpha
        for number, alpha in enumerate(alphas)
        if (0.85 <= alpha <= 1.5 and number % 2 == 0)
    ]

    for alpha in subselected_alphas:
        subselection = dataset.loc[dataset["alpha"] == alpha]
        ys = subselection[["epsilon"]]
        zs = subselection[[selector]]

        trasform = lambda alpha: (alpha - min(subselected_alphas)) / (
                max(subselected_alphas) - min(subselected_alphas)
        )
        color = color_map(trasform(alpha))
        row_size = 100
        try:
            ax.plot(
                ys.values,
                zs.values,
                color=color,
                linewidth=3,
                label=r"$\alpha={}$".format(round(alpha, 3)),
            )
        except Exception as exc:
            print(f"{exc}: Problem D=")

    plt.legend(loc=0, ncol=3, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_TE_errorbar(
        dataset,
        selector,
        error_selector,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        view=(70, 120),
        cmap="rainbow",
        dpi=300,
        fontsize=15,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1)

    markers = ["b", "^"]

    ax.set_title(title)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    alphas = dataset["alpha"].unique()
    mean = int(len(alphas) / 2)
    neghborhood = 5
    subselected_alphas = [
        alpha
        for number, alpha in enumerate(alphas)
        if (0.70 <= alpha <= 2 and number % 2 == 0)
    ]

    for alpha in subselected_alphas:
        subselection = dataset.loc[dataset["alpha"] == alpha]
        ys = subselection[["epsilon"]]
        zs = subselection[[selector]]
        error_bar = subselection[[error_selector]].copy()

        error_selector_negative_std = list(error_selector)
        error_selector_negative_std[1] = "-std"
        # error_bar[tuple(error_selector_negative_std)] = error_bar.apply(lambda x: -x, axis=1, raw=True)
        errors = error_bar.copy().T.to_numpy()

        trasform = lambda alpha: (alpha - min(subselected_alphas)) / (
                max(subselected_alphas) - min(subselected_alphas)
        )
        color = color_map(trasform(alpha))
        row_size = 100
        try:
            label = r"${\tiny" + rf"\alpha={round(alpha, 3)}" + r"}$"
            ax.errorbar(
                ys.values.flatten(),
                zs.values.flatten(),
                yerr=errors.flatten(),
                color=color,
                linewidth=3,
                label=label,
            )
        except Exception as exc:
            print(f"{exc}: {errors.shape}")

    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.legend(loc=0, ncol=3, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def figures2d_samples_TE(
        dataset, selector, title, ylabel, filename, suffix, cmap="rainbow", dpi=300
):
    matplotlib.style.use("seaborn")

    color_map = matplotlib.colormaps.get_cmap(cmap)
    alphas = dataset["alpha"].unique()
    epsilons = dataset["epsilon"].unique()
    subselection = dataset.loc[dataset["alpha"] == alphas[0]]
    subselection = subselection.loc[subselection["epsilon"] == epsilons[0]]

    one_subselection = subselection[[selector]]
    number_of_samples = len(subselection[[selector]].values[0, 0])
    mean = int(len(alphas) / 2)
    neghborhood = 5
    subselected_alphas = alphas[mean - neghborhood: mean + neghborhood]

    for sample in range(number_of_samples):
        fig = plt.figure(figsize=(13, 8))
        ax = fig.add_subplot(1, 1, 1)

        markers = ["b", "^"]

        ax.set_title(title)
        ax.set_xlabel(r"$\varepsilon$")
        ax.set_ylabel(ylabel)

        for alpha in subselected_alphas:
            subselection = dataset.loc[dataset["alpha"] == alpha]
            subselection.sort_values(by=["epsilon"], inplace=True)
            # print(subselection)
            ys = subselection[["epsilon"]]
            zs = subselection[[selector]]

            trasform = lambda alpha: (alpha - min(subselected_alphas)) / (
                    max(subselected_alphas) - min(subselected_alphas)
            )
            color = color_map(trasform(alpha))
            row_size = 100
            try:
                ax.plot(
                    ys.values,
                    [float(item[0][sample]) for item in zs.values],
                    color=color,
                    linewidth=3,
                    label=r"$\alpha={}$".format(round(alpha, 3)),
                )
            except Exception as exc:
                print(f"{exc}: Problem D=")

        plt.legend(loc=4)
        plt.savefig(
            filename.format(sample) + "." + suffix, dpi=dpi, bbox_inches="tight"
        )
        plt.close()
        del fig


def figures2d_fixed_epsilon(
        dataset,
        selector,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=15,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)
    color_map = matplotlib.colormaps.get_cmap(cmap)

    alphas = dataset["alpha"].unique()
    fixed_epsilon = 0.07

    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    subselected_alphas = [
        alpha
        for number, alpha in enumerate(alphas)
        if (0.70 <= alpha <= 2 and number % 2 == 0)
    ]

    constants = []
    for alpha in subselected_alphas:
        subselection = dataset.loc[dataset["alpha"] == alpha]
        constant = (
                subselection.loc[dataset["epsilon"] == fixed_epsilon][[selector]].values[0][
                    0
                ]
                - 1.1
        )
        constants.append(constant)

        ys = subselection[["epsilon"]].values
        zs = [[constant] for item in ys]

        trasform = lambda alpha: (alpha - min(alphas)) / (max(alphas) - min(alphas))
        color = color_map(trasform(alpha))
        try:
            ax.plot(
                ys,
                zs,
                color=color,
                linewidth=3,
                label=r"$\varepsilon={}$".format(round(alpha, 3)),
            )
        except Exception as exc:
            print(f"{exc}: Problem D=")

    ax.set_ylim(
        min(constants), max(constants) + (max(constants) - min(constants)) * 0.3
    )
    plt.legend(ncol=6, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    plt.savefig(filename + "_f" + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def escort_distribution(
        datasets,
        columns,
        title,
        xlabel,
        ylabel,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.colormaps.get_cmap(cmap)

    fig = plt.figure(figsize=(13, 8))
    markers = ["b", "^"]
    fig.suptitle(title)

    for index_dataset, dataset in enumerate(datasets):
        ax = fig.add_subplot(1, len(datasets), index_dataset + 1)
        if index_dataset == 0:
            ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)

        for index_column, column in enumerate(columns):
            subselection_x = dataset[["x"]]
            subselection_y = dataset[[str(column)]]

            try:
                color = color_map(index_column / (len(columns) - 1))
                ax.set_yscale("log")
                if index_dataset == 1:
                    ax.set_xlim(-0.45, 0.2)
                    ax.set_ylim(0.0000000000001, 0.08)
                else:
                    ax.set_xlim(-15, 15)
                    ax.set_ylim(0.00000001, 2)
                ax.plot(
                    subselection_x.values,
                    subselection_y.values,
                    color=color,
                    linewidth=2,
                    label=r"$\alpha={}$".format(round(column, 3)),
                )
            except Exception as exc:
                print(f"{exc}: Problem D=")
        if index_dataset == 0:
            ax.legend(loc=1)

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def granger_function_plot(
        dataset,
        title,
        xlabel,
        ylabel,
        zlabel,
        filename,
        suffix,
        cmap="rainbow",
        view=(50, -20),
        dpi=300,
):
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_zlim(-5, 0)
    ax.set_ylim(1, 10)
    # ax.set_zscale("log")

    row_size = len(dataset["k"].unique())
    xs = dataset["alpha"]
    ys = dataset["k"]
    zs = dataset["granger"]
    Xs = np.reshape(xs.values, (row_size, -1))
    Ys = np.reshape(ys.values, (row_size, -1))
    Zs = np.reshape(zs.values, (row_size, -1))

    try:
        surf = ax.plot_surface(
            Xs, Ys, Zs, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False
        )
        fig.colorbar(surf, shrink=0.5, aspect=10)
    except Exception as exc:
        print(f"{exc}: Problem D=")

    # plt.legend(loc=1)
    ax.view_init(view[0], view[1])

    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()

    del fig


def lyapunov_exponent_plot(
        dataset,
        title,
        xlabel,
        ylabels,
        labels,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=17,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig = plt.figure(figsize=(13, 8))

    x = dataset[[0]].values.flatten().tolist()
    ys = [dataset[[i]].values.flatten().tolist() for i in [1, 5, 4]]
    y0 = [0.0] * len(x)
    ys.insert(1, y0)
    # fig.set_title(title)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabels)
    ax.set_xlim(0, 0.25)
    ax.set_ylim(-0.13, 0.135)

    for index, y in enumerate(ys):
        color = color_map(index / (len(ys) - 1))
        ax.plot(x, y, linewidth=3, color=color, label=labels[index])

    # plt.legend(loc=3)
    ax.set_xticklabels([0.0, 0.05, 0.1, 0.15, 0.2, 0.25], fontsize=fontsize)
    ax.set_yticklabels([-0.1, -0.05, 0, 0.05, 0.15, 0.2], fontsize=fontsize)
    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def qgauss_plot(
        dataset,
        title,
        xlabel,
        ylabels,
        labels,
        filename,
        suffix,
        cmap="rainbow",
        dpi=300,
        fontsize=17,
        style="seaborn-v0_8",
):
    matplotlib.style.use(style)

    color_map = matplotlib.colormaps.get_cmap(cmap)
    fig = plt.figure(figsize=(13, 8))

    selected_datasets = dataset[dataset["X1"] >= 0.001]
    x = selected_datasets["X1"].values.flatten().tolist()
    ys = [
        selected_datasets[[selected_datasets.columns[i]]].values.flatten().tolist()
        for i in [1, 2, 3]
    ]
    fig.suptitle(title)
    ax = fig.add_subplot(1, 2, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabels)

    # ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 3.7)
    ax.set_ylim(1e-8, 1e-3)

    for index, y in enumerate(ys):
        color = color_map(index / (len(ys) - 1))
        ax.plot(x, y, linewidth=3, color=color, label=labels[index])

    plt.legend()

    ax = fig.add_subplot(1, 2, 2)
    ax.set_xlabel(xlabel)
    # ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 3.7)
    ax.set_ylim(1e-8, 1e-3)

    selected_datasets = dataset[dataset["X1"] >= 0.001]
    x = selected_datasets["X1"].values.flatten().tolist()
    ys = [
        selected_datasets[[selected_datasets.columns[i]]].values.flatten().tolist()
        for i in [5, 6, 7]
    ]

    for index, y in enumerate(ys):
        color = color_map(index / (len(ys) - 1))
        ax.plot(x, y, linewidth=3, color=color, label=labels[index + 3])

    plt.legend()
    plt.savefig(filename + "." + suffix, dpi=dpi, bbox_inches="tight")
    plt.close()
    del fig


def fill_column_frame_RTE(
        source_frame,
        processed_column,
        take_k_th_nearest_neighbor,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    if len(unary_operation) == len(destination_multiindex):
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            calculation = source_frame.apply(
                lambda row: actual_unary_operation(
                    row[processed_column][take_k_th_nearest_neighbor:]
                ),
                axis=1,
                raw=False,
            )
            calculated_frame[actual_destination_multiindex] = calculation
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_statistics_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    samples = source_frame[[processed_columns[0]]].iloc[0].values[0].shape[0]
    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        for sample in range(samples):
            sample_actual_destination_multiindex = list(
                actual_destination_multiindex[:]
            )
            sample_actual_destination_multiindex[7] = sample
            calculated_frame[[tuple(sample_actual_destination_multiindex)]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            np.array(row[[processed_column]]) for processed_column in processed_columns
        ]
        merged_columns = np.stack(process_columns, axis=1)
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(merged_columns.tolist())
            for sample in range(samples):
                sample_actual_destination_multiindex = list(
                    actual_destination_multiindex[:]
                )
                sample_actual_destination_multiindex[7] = sample

                calculated_frame.at[
                    index, tuple(sample_actual_destination_multiindex)
                ] = result[0][sample]
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_NN_statistics_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
        take_k_th_nearest_neighbor,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[[actual_destination_multiindex]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            row[[processed_column]]
            .values.tolist()[0]
            .tolist()[take_k_th_nearest_neighbor:]
            for processed_column in processed_columns
        ]

        dataset = np.array(process_columns).flatten()
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(dataset)

            calculated_frame.at[index, actual_destination_multiindex] = result
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_statistics_balanced_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    samples = source_frame[[processed_columns[0]]].iloc[0].values[0].shape[0]
    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        for sample in range(samples):
            sample_actual_destination_multiindex = list(
                actual_destination_multiindex[:]
            )
            sample_actual_destination_multiindex[7] = sample
            calculated_frame[[tuple(sample_actual_destination_multiindex)]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            np.array(row[[processed_column]])
            - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            for processed_column in processed_columns
        ]
        merged_columns = np.stack(process_columns, axis=1)
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(merged_columns.tolist())
            for sample in range(samples):
                sample_actual_destination_multiindex = list(
                    actual_destination_multiindex[:]
                )
                sample_actual_destination_multiindex[7] = sample

                calculated_frame.at[
                    index, tuple(sample_actual_destination_multiindex)
                ] = np.array(result[0])[sample]
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_statistics_effective_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    samples = source_frame[[processed_columns[0]]].iloc[0].values[0].shape[0]
    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        for sample in range(samples):
            sample_actual_destination_multiindex = list(
                actual_destination_multiindex[:]
            )
            sample_actual_destination_multiindex[7] = sample
            calculated_frame[[tuple(sample_actual_destination_multiindex)]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            np.array(row[[processed_column]])
            - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            for processed_column in processed_columns
        ]

        merged_columns = np.stack(process_columns, axis=1)
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(merged_columns.tolist())
            for sample in range(samples):
                sample_actual_destination_multiindex = list(
                    actual_destination_multiindex[:]
                )
                sample_actual_destination_multiindex[7] = sample

                calculated_frame.at[
                    index, tuple(sample_actual_destination_multiindex)
                ] = np.array(result[0])[sample]
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_NN_statistics_effective_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
        take_k_th_nearest_neighbor,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[destination_multiindex] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            (
                    np.array(row[[processed_column]])
                    - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            )
            .tolist()[0]
            .tolist()[take_k_th_nearest_neighbor:]
            for processed_column in processed_columns
        ]
        dataset = np.array(process_columns).flatten()
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(dataset)
            calculated_frame.at[index, actual_destination_multiindex] = result
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_NN_statistics_balanced_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
        take_k_th_nearest_neighbor,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    samples = source_frame[[processed_columns[0]]].iloc[0].values[0].shape[0]
    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[[actual_destination_multiindex]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            (
                    np.array(row[[processed_column]])
                    - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            )
            .tolist()[0]
            .tolist()[take_k_th_nearest_neighbor:]
            for processed_column in processed_columns
        ]

        dataset = np.array(process_columns).flatten()
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(dataset)
            calculated_frame.at[index, actual_destination_multiindex] = result
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_NN_statistics_balanced_effective_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
        take_k_th_nearest_neighbor,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[destination_multiindex] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            (
                    np.array(
                        row[
                            [
                                (
                                    processed_column[0],
                                    processed_column[1],
                                    processed_column[2],
                                    processed_column[3],
                                    not processed_column[4],
                                    processed_column[5],
                                    0,
                                )
                            ]
                        ]
                    )
                    - np.array(row[[processed_column]])
                    - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            not processed_column[5],
                            0,
                        )
                    ]
                ]
            )
                    + np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            processed_column[6],
                        )
                    ]
                ]
            )
            )
            .tolist()[0]
            .tolist()[take_k_th_nearest_neighbor:]
            for processed_column in processed_columns
        ]
        dataset = np.array(process_columns).flatten()
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(dataset)
            calculated_frame.at[index, actual_destination_multiindex] = result
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_sample_statistics_balanced_effective_RTE(
        source_frame,
        processed_columns,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    samples = source_frame[[processed_columns[0]]].iloc[0].values[0].shape[0]
    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        for sample in range(samples):
            sample_actual_destination_multiindex = list(
                actual_destination_multiindex[:]
            )
            sample_actual_destination_multiindex[7] = sample
            calculated_frame[[tuple(sample_actual_destination_multiindex)]] = np.NAN
    for index, row in source_frame.iterrows():
        process_columns = [
            np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            - np.array(row[[processed_column]])
            - np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            not processed_column[5],
                            0,
                        )
                    ]
                ]
            )
            + np.array(
                row[
                    [
                        (
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            processed_column[6],
                        )
                    ]
                ]
            )
            for processed_column in processed_columns
        ]
        merged_columns = np.stack(process_columns, axis=1)
        for actual_unary_operation, actual_destination_multiindex in zip(
                unary_opertions, destination_multiindices
        ):
            result = actual_unary_operation(merged_columns.tolist())
            for sample in range(samples):
                sample_actual_destination_multiindex = list(
                    actual_destination_multiindex[:]
                )
                sample_actual_destination_multiindex[7] = sample

                calculated_frame.at[
                    index, tuple(sample_actual_destination_multiindex)
                ] = np.array(result[0])[sample]
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_effective_RTE(
        source_frame: object,
        processed_column: object,
        take_k_th_nearest_neighbor: object,
        unary_operation: object,
        destination_frame: object,
        destination_multiindex: object,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[actual_destination_multiindex] = source_frame.apply(
            lambda row: float(
                actual_unary_operation(
                    np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            processed_column[5],
                            processed_column[6],
                        ][take_k_th_nearest_neighbor:]
                    )
                    - np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            processed_column[5],
                            0,
                        ][take_k_th_nearest_neighbor:]
                    )
                )
            ),
            axis=1,
            raw=False,
        )
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_balance_RTE(
        source_frame,
        processed_column,
        take_k_th_nearest_neighbor,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[actual_destination_multiindex] = source_frame.apply(
            lambda row: float(
                actual_unary_operation(
                    np.array(row[processed_column][take_k_th_nearest_neighbor:])
                    - np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            processed_column[6],
                        ][take_k_th_nearest_neighbor:]
                    )
                )
            ),
            axis=1,
            raw=False,
        )
    return pd.concat([destination_frame, calculated_frame], axis=1)


def fill_column_frame_balance_effective_RTE(
        source_frame,
        processed_column,
        take_k_th_nearest_neighbor,
        unary_operation,
        destination_frame,
        destination_multiindex,
):
    if isinstance(unary_operation, list or tuple) and isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = unary_operation
        destination_multiindices = destination_multiindex
    elif not isinstance(unary_operation, list or tuple) and not isinstance(
            destination_multiindex, list or tuple
    ):
        unary_opertions = [unary_operation]
        destination_multiindices = [destination_multiindex]
    else:
        raise RuntimeError(
            "Unary operations and destination multi-indices do not have "
        )

    calculated_frame = pd.DataFrame(index=destination_frame.index)
    for actual_unary_operation, actual_destination_multiindex in zip(
            unary_opertions, destination_multiindices
    ):
        calculated_frame[actual_destination_multiindex] = source_frame.apply(
            lambda row: float(
                actual_unary_operation(
                    -np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            processed_column[5],
                            processed_column[6],
                        ][take_k_th_nearest_neighbor:]
                    )
                    + np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            processed_column[5],
                            0,
                        ][take_k_th_nearest_neighbor:]
                    )
                    + np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            processed_column[4],
                            not processed_column[5],
                            processed_column[6],
                        ][take_k_th_nearest_neighbor:]
                    )
                    - np.array(
                        row[
                            processed_column[0],
                            processed_column[1],
                            processed_column[2],
                            processed_column[3],
                            not processed_column[4],
                            not processed_column[5],
                            0,
                        ][take_k_th_nearest_neighbor:]
                    )
                )
            ),
            axis=1,
            raw=False,
        )
    return pd.concat([destination_frame, calculated_frame], axis=1)


def refined_process_dataset(input_dataset, processed_dataset):
    new_columns_base_name = "transfer_entropy"
    effective_new_column_name = f"effective_{new_columns_base_name}"
    balance_new_column_name = f"balance_{new_columns_base_name}"
    balance_effective_new_column_name = f"balance_effective_{new_columns_base_name}"

    take_k_th_nearest_neighbor = 5
    index_shuffle_dataset = 4
    index_swap_dataset = 5
    converter_epsilon = converter_epsilon = lambda x: (
        x.split("-")[2].split(".b")[0] if "1e-" in x else x.split("-")[1].split(".b")[0]
    )

    averaging_over_samples = f"Avaraging over sample for each neighbour k"
    averaging_over_samples_and_neighbors = f"Avaraging over samples and neighbours k"
    averaging_over_neighbors = (
        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}"
    )

    epsilon = converter_epsilon(input_dataset)
    path = Path(input_dataset)
    print(
        f"{datetime.datetime.now().isoformat()} Processing input file {input_dataset}"
    )
    table = pd.read_pickle(path)

    refined_table_converted_multiindices = {}
    for key, value in table.items():
        refined_table_converted_multiindices[key] = {}
        for key_q, items in value.items():
            refined_key = key
            for index, value in enumerate(items):
                refined_table_converted_multiindices[refined_key][
                    (index, key_q)
                ] = value

    frame_refined = pd.DataFrame.from_dict(refined_table_converted_multiindices)
    frame_refined.index.set_names(["k", "q"], inplace=True)

    refined_table_values_arrays = {}
    for key, value in table.items():
        for key_q, items in value.items():
            if key in refined_table_values_arrays:
                refined_table_values_arrays[key][key_q] = np.array(items)
            else:
                refined_table_values_arrays[key] = {key_q: np.array(items)}

    frame = pd.DataFrame.from_dict(
        refined_table_values_arrays
    )  # , columns=refined_column
    frame["epsilon"] = epsilon

    index_dataframe = frame.index
    names = [
        "Name of variable",
        "History first",
        "Future first",
        "History second",
        "Statistical value",
        "Shuffled",
        "Reversed direction",
        "Sample number",
        "Remark",
        "Sample",
    ]
    refined_column = pd.MultiIndex(
        levels=[[]] * len(names), codes=[[]] * len(names), names=names
    )
    refined_frame = pd.DataFrame(index=index_dataframe, columns=refined_column)

    # preparing columns aggregated using samples
    multi_sample_columns = {}
    for item in frame.columns[:-1]:
        (
            variable,
            history_first,
            future_first,
            history_second,
            shuffled_dataset,
            reversed_order,
            sample_number,
        ) = item
        if (
                not (
                            variable,
                            history_first,
                            future_first,
                            history_second,
                            shuffled_dataset,
                            reversed_order,
                    )
                    in multi_sample_columns
        ):
            multi_sample_columns[
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                )
            ] = [sample_number]
        else:
            multi_sample_columns[
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                )
            ].append(sample_number)

    unary_statistical_operations = [
        lambda x: np.mean(x, axis=0),
        lambda x: np.std(x, axis=0),
        lambda x: np.median(x, axis=0),
        lambda x: np.quantile(x, 0.25, axis=0),
        lambda x: np.quantile(x, 0.75, axis=0),
        lambda x: np.quantile(x, 0.1, axis=0),
        lambda x: np.quantile(x, 0.9, axis=0),
    ]
    operations_names = [
        "mean",
        "std",
        "median",
        "quantile 0.25",
        "quantile 0.75",
        "quantile 0.1",
        "quantile 0.9",
    ]
    del item

    # sample statistics
    sample_number = 0
    for key, samples in multi_sample_columns.items():
        (
            variable,
            history_first,
            future_first,
            history_second,
            shuffled_dataset,
            reversed_order,
        ) = key
        columns = [list(key) + [sample] for sample in samples]

        refined_frame = fill_column_frame_sample_statistics_RTE(
            frame,
            columns,
            unary_statistical_operations,
            refined_frame,
            [
                (
                    new_columns_base_name,
                    history_first,
                    future_first,
                    history_second,
                    operation_name,
                    shuffled_dataset,
                    reversed_order,
                    sample_number,
                    averaging_over_samples,
                    epsilon,
                )
                for operation_name in operations_names
            ],
        )

        if shuffled_dataset:
            refined_frame = fill_column_frame_sample_statistics_effective_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_samples,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

        if not reversed_order:
            refined_frame = fill_column_frame_sample_statistics_balanced_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        balance_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_samples,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

        if not reversed_order and shuffled_dataset:
            refined_frame = fill_column_frame_sample_statistics_balanced_effective_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        balance_effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_samples,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

    # nearest neighbour statistics
    columns_to_calculate_RTE = frame.columns
    for key in columns_to_calculate_RTE[:-1]:
        (
            variable,
            history_first,
            future_first,
            history_second,
            shuffled_dataset,
            reversed_order,
            sample_number,
        ) = key

        refined_frame = fill_column_frame_RTE(
            frame,
            key,
            take_k_th_nearest_neighbor,
            unary_statistical_operations,
            refined_frame,
            [
                (
                    new_columns_base_name,
                    history_first,
                    future_first,
                    history_second,
                    operation_name,
                    shuffled_dataset,
                    reversed_order,
                    sample_number,
                    averaging_over_neighbors,
                    epsilon,
                )
                for operation_name in operations_names
            ],
        )

        if shuffled_dataset:
            refined_frame = fill_column_frame_effective_RTE(
                frame,
                key,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

        if not reversed_order:
            refined_frame = fill_column_frame_balance_RTE(
                frame,
                key,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

        if not reversed_order and shuffled_dataset:
            refined_frame = fill_column_frame_balance_effective_RTE(
                frame,
                key,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        balance_effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        averaging_over_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
            )

    # sample and nearest neighbour statistics
    for key, samples in multi_sample_columns.items():
        (
            variable,
            history_first,
            future_first,
            history_second,
            shuffled_dataset,
            reversed_order,
        ) = key
        columns = [list(key) + [sample] for sample in samples]

        if "information" not in str(key[0]):  # more sample present
            refined_frame = fill_column_frame_sample_NN_statistics_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        new_columns_base_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        shuffled_dataset,
                        reversed_order,
                        0,
                        averaging_over_samples_and_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
                take_k_th_nearest_neighbor,
            )

        if (not reversed_order) and "information" not in str(key[0]):
            refined_frame = fill_column_frame_sample_NN_statistics_balanced_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        balance_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        shuffled_dataset,
                        reversed_order,
                        0,
                        averaging_over_samples_and_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
                take_k_th_nearest_neighbor,
            )

        if shuffled_dataset and "information" not in str(key[0]):
            refined_frame = fill_column_frame_sample_NN_statistics_effective_RTE(
                frame,
                columns,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        effective_new_column_name,
                        history_first,
                        future_first,
                        history_second,
                        operation_name,
                        not shuffled_dataset,
                        reversed_order,
                        0,
                        averaging_over_samples_and_neighbors,
                        epsilon,
                    )
                    for operation_name in operations_names
                ],
                take_k_th_nearest_neighbor,
            )

        if not reversed_order and shuffled_dataset and "information" not in str(key[0]):
            refined_frame = (
                fill_column_frame_sample_NN_statistics_balanced_effective_RTE(
                    frame,
                    columns,
                    unary_statistical_operations,
                    refined_frame,
                    [
                        (
                            balance_effective_new_column_name,
                            history_first,
                            future_first,
                            history_second,
                            operation_name,
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            averaging_over_samples_and_neighbors,
                            epsilon,
                        )
                        for operation_name in operations_names
                    ],
                    take_k_th_nearest_neighbor,
                )
            )

    # append frame for processing
    refined_frame.to_pickle(processed_dataset)


def process_datasets(
        directory,
        processed_datasets,
        result_dataset,
        result_raw_dataset,
        new_columns_base_name="transfer_entropy",
        take_k_th_nearest_neighbor=5,
        converter_epsilon=lambda file: float(
            "-" + file.split("--")[1].split(".b")[0]
            if "--" in file
            else file.split("-")[1].split(".b")[0]
        ),
):
    # taking only some nn data to assure that it converge in theory
    files = glob.glob(processed_datasets)
    print(files)
    frames = []
    frames_raw = []
    frames_refined_statistics = []
    index_shuffle_dataset = 4
    index_swap_dataset = 5

    for file in files:
        epsilon = converter_epsilon(file)
        path = Path(file)
        print(f"{datetime.datetime.now().isoformat()} {file}")
        table = pd.read_pickle(path)

        refined_table_converted_multiindices = {}
        for key, value in table.items():
            refined_table_converted_multiindices[key] = {}
            for key_q, items in value.items():
                refined_key = key
                # if refined_key in refined_table:
                #    refined_table[refined_key][key_q] = np.array(items)
                # else:
                #    refined_table[refined_key] = {key_q: np.array(items)}
                for index, value in enumerate(items):
                    refined_table_converted_multiindices[refined_key][
                        (index, key_q)
                    ] = value

        frame_refined = pd.DataFrame.from_dict(refined_table_converted_multiindices)
        frame_refined.index.set_names(["k", "q"], inplace=True)

        refined_table_values_arrays = {}
        for key, value in table.items():
            for key_q, items in value.items():
                if key in refined_table_values_arrays:
                    refined_table_values_arrays[key][key_q] = np.array(items)
                else:
                    refined_table_values_arrays[key] = {key_q: np.array(items)}

        frame = pd.DataFrame.from_dict(
            refined_table_values_arrays
        )  # , columns=refined_column
        frame["epsilon"] = epsilon

        index_dataframe = frame.index
        names = [
            "Name of variable",
            "History first",
            "Future first",
            "History second",
            "Statistical value",
            "Shuffled",
            "Reversed direction",
            "Sample number",
            "Remark",
            "Sample",
        ]
        refined_column = pd.MultiIndex(
            levels=[[]] * len(names), codes=[[]] * len(names), names=names
        )
        refined_frame = pd.DataFrame(index=index_dataframe, columns=refined_column)

        # preparing columns aggregated using samples
        multi_sample_columns = {}
        for item in frame.columns[:-1]:
            (
                variable,
                history_first,
                future_first,
                history_second,
                shuffled_dataset,
                reversed_order,
                sample_number,
            ) = item
            if (
                    not (
                                variable,
                                history_first,
                                future_first,
                                history_second,
                                shuffled_dataset,
                                reversed_order,
                        )
                        in multi_sample_columns
            ):
                multi_sample_columns[
                    (
                        variable,
                        history_first,
                        future_first,
                        history_second,
                        shuffled_dataset,
                        reversed_order,
                    )
                ] = [sample_number]
            else:
                multi_sample_columns[
                    (
                        variable,
                        history_first,
                        future_first,
                        history_second,
                        shuffled_dataset,
                        reversed_order,
                    )
                ].append(sample_number)

        unary_statistical_operations = [
            lambda x: np.mean(x, axis=0),
            lambda x: np.std(x, axis=0),
            lambda x: np.median(x, axis=0),
            lambda x: np.quantile(x, 0.25, axis=0),
            lambda x: np.quantile(x, 0.75, axis=0),
            lambda x: np.quantile(x, 0.1, axis=0),
            lambda x: np.quantile(x, 0.9, axis=0),
        ]

        # drop aggregated elements with only single sample
        for key, samples in multi_sample_columns.items():
            if len(samples) > 1:  # more sample present
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                ) = key
                columns = [list(key) + [sample] for sample in samples]

                refined_frame = fill_column_frame_sample_statistics_RTE(
                    frame,
                    columns,
                    unary_statistical_operations,
                    refined_frame,
                    [
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "mean",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "std",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "median",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.25",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.75",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.1",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.9",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                    ],
                )

            if (
                    (len(samples) > 1)
                    and shuffled_dataset
                    and not ("entropy" in str(key[0]) or "information" in str(key[0]))
            ):
                effective_columns_base_name = f"effective_{new_columns_base_name}"
                refined_frame = fill_column_frame_sample_statistics_RTE(
                    frame,
                    columns,
                    unary_statistical_operations,
                    refined_frame,
                    [
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "mean",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "std",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "median",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.25",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.75",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.1",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            effective_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.9",
                            not shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                    ],
                )

            if (not reversed_order) and "information" not in str(key[0]):
                balance_column_name = f"balance_{new_columns_base_name}"
                refined_frame = fill_column_frame_sample_statistics_balanced_RTE(
                    frame,
                    columns,
                    unary_statistical_operations,
                    refined_frame,
                    [
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "mean",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "std",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "median",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.25",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.75",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.1",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                        (
                            balance_column_name,
                            history_first,
                            future_first,
                            history_second,
                            "quantile 0.9",
                            shuffled_dataset,
                            reversed_order,
                            0,
                            f"Avaraging over sample for each neighbour k",
                            epsilon,
                        ),
                    ],
                )

            if (
                    (len(samples) > 1)
                    and not reversed_order
                    and shuffled_dataset
                    and "information" not in str(item[0])
            ):
                balance_effective_column_name = (
                    f"balance_effective_{new_columns_base_name}"
                )
                refined_frame = (
                    fill_column_frame_sample_statistics_balanced_effective_RTE(
                        frame,
                        columns,
                        unary_statistical_operations,
                        refined_frame,
                        [
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "mean",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "std",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "median",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "quantile 0.25",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "quantile 0.75",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "quantile 0.1",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                            (
                                balance_effective_column_name,
                                history_first,
                                future_first,
                                history_second,
                                "quantile 0.9",
                                not shuffled_dataset,
                                reversed_order,
                                0,
                                f"Avaraging over sample for each neighbour k",
                                epsilon,
                            ),
                        ],
                    )
                )

        columns_to_calculate_RTE = frame.columns
        for item in columns_to_calculate_RTE[:-1]:
            try:
                (
                    variable,
                    history_first,
                    future_first,
                    history_second,
                    shuffled_dataset,
                    reversed_order,
                    sample_number,
                ) = item
                print(f"Processing {item}")

                refined_frame = fill_column_frame_RTE(
                    frame,
                    item,
                    take_k_th_nearest_neighbor,
                    unary_statistical_operations,
                    refined_frame,
                    [
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "mean",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "std",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "median",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "0.25 quantile",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "0.75 quantile",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "0.1 quantile",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                        (
                            new_columns_base_name,
                            history_first,
                            future_first,
                            history_second,
                            "0.9 quantile",
                            shuffled_dataset,
                            reversed_order,
                            sample_number,
                            f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                            epsilon,
                        ),
                    ],
                )

                frame[
                    f"{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                    "mean",
                    "",
                    "",
                    shuffled_dataset,
                    reversed_order,
                    sample_number,
                ] = frame.apply(
                    lambda row: np.mean(row[item][take_k_th_nearest_neighbor:]),
                    axis=1,
                    raw=False,
                )
                frame[
                    f"{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                    "std",
                    "",
                    "",
                    shuffled_dataset,
                    reversed_order,
                    sample_number,
                ] = frame.apply(
                    lambda row: np.std(row[item][take_k_th_nearest_neighbor:]),
                    axis=1,
                    raw=False,
                )

                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.std(x, axis=0), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "std", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))
                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.median(x, axis=0), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "median", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))
                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.quantile(x, 0.25), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "0.25 quantile", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))
                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.quantile(x, 0.75), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "0.75 quantile", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))
                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.quantile(x, 0.1), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "0.1 quantile", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))
                # fill_column_frame(frame, item, take_k_th_nearest_neighbor, lambda x: np.quantile(x, 0.9), refined_frame, (new_columns_base_name, history_first, future_first, history_second, "0.9 quantile", shuffled_dataset, reversed_order, sample_number, f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}", variable))

            except Exception as exc:
                debug = frame[item]
                print(f"{exc} {item} {debug}")
                raise exc

        # effective transfer entropy
        column_to_calculate_effective_transfer_entropy = [
            item
            for item in frame.columns.tolist()
            if item[index_shuffle_dataset]
               and not ("entropy" in str(item[0]) or "information" in str(item[0]))
        ]
        for item in column_to_calculate_effective_transfer_entropy:
            (
                variable,
                history_first,
                future_first,
                history_second,
                shuffled_dataset,
                reversed_order,
                sample_number,
            ) = item

            refined_frame = fill_column_frame_effective_RTE(
                frame,
                item,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "mean",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "std",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "median",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.25 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.75 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.1 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.9 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                ],
            )

            frame[
                f"effective_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "mean",
                "",
                "",
                not shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.mean(
                        np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                        - np.array(row[item][take_k_th_nearest_neighbor:])
                    )
                ),
                axis=1,
                raw=False,
            )
            frame[
                f"effective_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "std",
                "",
                "",
                not shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.std(
                        np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                        - np.array(row[item][take_k_th_nearest_neighbor:])
                    )
                ),
                axis=1,
                raw=False,
            )

        # balance of entropy
        column_to_calculate_balance_transfer_entropy = [
            item
            for item in frame.columns.tolist()
            if not bool(item[index_swap_dataset])
               and "information" not in str(item[0])
               and "epsilon" not in str(item[0])
        ]
        for item in column_to_calculate_balance_transfer_entropy:
            (
                variable,
                history_first,
                future_first,
                history_second,
                shuffled_dataset,
                reversed_order,
                sample_number,
            ) = item

            refined_frame = fill_column_frame_effective_RTE(
                frame,
                item,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "mean",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "std",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "median",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.25 quantile",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.75 quantile",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.1 quantile",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.9 quantile",
                        shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                ],
            )

            frame[
                f"balance_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "mean",
                "",
                "",
                shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.mean(
                        np.array(row[item][take_k_th_nearest_neighbor:])
                        - np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                shuffled_dataset,
                                not reversed_order,
                                sample_number,
                            ][take_k_th_nearest_neighbor:]
                        )
                    )
                ),
                axis=1,
                raw=False,
            )
            frame[
                f"balance_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "std",
                "",
                "",
                shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.std(
                        np.array(row[item][take_k_th_nearest_neighbor:])
                        - np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                shuffled_dataset,
                                not reversed_order,
                                sample_number,
                            ][take_k_th_nearest_neighbor:]
                        )
                    )
                ),
                axis=1,
                raw=False,
            )

        print(f"Number of columns {len(refined_frame.columns)}")

        # balance of effective entropy
        balance_effective_names = [
            item
            for item in frame.columns.tolist()
            if not bool(item[index_swap_dataset])
               and bool(item[index_shuffle_dataset])
               and "information" not in str(item[0])
               and "epsilon" not in str(item[0])
        ]
        for item in balance_effective_names:
            (
                variable,
                history_first,
                future_first,
                history_second,
                shuffled_dataset,
                reversed_order,
                sample_number,
            ) = item

            refined_frame = fill_column_frame_balance_effective_RTE(
                frame,
                item,
                take_k_th_nearest_neighbor,
                unary_statistical_operations,
                refined_frame,
                [
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "mean",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "std",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "median",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.25 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.75 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.1 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                    (
                        f"balance_effective_{new_columns_base_name}",
                        history_first,
                        future_first,
                        history_second,
                        "0.9 quantile",
                        not shuffled_dataset,
                        reversed_order,
                        sample_number,
                        f"Avaraging over nearest neighbours starting from {take_k_th_nearest_neighbor}",
                        epsilon,
                    ),
                ],
            )

            frame[
                f"balance_effective_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "mean",
                "",
                "",
                not shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.mean(
                        -np.array(row[item][take_k_th_nearest_neighbor:])
                        + np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                        + np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                shuffled_dataset,
                                not reversed_order,
                                sample_number,
                            ][take_k_th_nearest_neighbor:]
                        )
                        - np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                not reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                    )
                ),
                axis=1,
                raw=False,
            )
            frame[
                f"balance_effective_{new_columns_base_name}_{history_first}_{future_first}_{history_second}",
                "std",
                "",
                "",
                not shuffled_dataset,
                reversed_order,
                sample_number,
            ] = frame.apply(
                lambda row: float(
                    np.std(
                        -np.array(row[item][take_k_th_nearest_neighbor:])
                        + np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                        + np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                shuffled_dataset,
                                not reversed_order,
                                sample_number,
                            ][take_k_th_nearest_neighbor:]
                        )
                        - np.array(
                            row[
                                item[0],
                                history_first,
                                future_first,
                                history_second,
                                not shuffled_dataset,
                                not reversed_order,
                                0,
                            ][take_k_th_nearest_neighbor:]
                        )
                    )
                ),
                axis=1,
                raw=False,
            )

        # dropping the index
        frame = frame.reset_index()

        column = [
            ("alpha", "", "", "", "", "", "") if "index" == item[0] else item
            for item in frame.columns.tolist()
        ]
        new_columns = pd.MultiIndex.from_tuples(
            [
                ("alpha", "", "", "", "", "", "") if "index" == item[0] else item
                for item in frame.columns
            ]
        )
        frame.columns = new_columns

        # selection of columns
        columns = [
            item
            for item in frame.columns.tolist()
            if "mean" in str(item[1])
               or "std" in str(item[1])
               or "alpha" in str(item[0])
               or "epsilon" in str(item[0])
        ]
        frame_with_processed_results = frame[columns]

        columns = [
            item
            for item in frame.columns.tolist()
            if isinstance(item[0], float)
               or "alpha" in str(item[0])
               or "epsilon" in str(item[0])
        ]
        frame_with_raw_results = frame[columns]
        # print(frame)
        # if item[0] not in ["alpha", "epsilon"] else item[0:3]
        columns = [
            (
                str(item[1]) + "_" + str(item[2]) + "_" + str(item[3])
                if isinstance(item[0], float)
                else item[0]
            )
            for item in frame_with_raw_results.columns.tolist()
        ]
        frame_with_raw_results.columns = columns

        # append frame for processing
        frames.append(frame_with_processed_results)
        frames_raw.append(frame_with_raw_results)
        frames_refined_statistics.append(refined_frame)

    # join the table
    join_table = pd.concat(frames, ignore_index=True)
    join_table_refined_statistics = pd.concat(frames_refined_statistics)
    try:
        join_table_raw = pd.concat(frames_raw, ignore_index=True)
    except:
        print("Problem with columns")
        first_frame = frames_raw[0]
        for file, frame in zip(files, frames_raw):
            comparison = len(frame.columns.tolist()) == len(
                first_frame.columns.tolist()
            )
            if not comparison:
                print(file, frame.columns, comparison)
        sys.exit(1)

    # print(join_table)
    index_alpha = join_table.columns.tolist()
    pivot_table = pd.pivot_table(join_table, index=[index_alpha[0], index_alpha[1]])
    print(pivot_table, join_table.columns.tolist())

    print(join_table_raw)
    index_alpha = join_table_raw.columns.tolist()
    pivot_table_raw = join_table_raw.set_index([index_alpha[0], index_alpha[-1]])
    print(pivot_table_raw)

    TE = pivot_table.reset_index()
    TE_raw = pivot_table_raw.reset_index()

    TE.to_pickle(result_dataset)
    TE_raw.to_pickle(result_raw_dataset)

    join_table_refined_statistics.to_pickle(f"{directory}/refined_statistics.bin")

    return (
        TE,
        [item for item in join_table.columns.tolist() if "mean" in str(item[1])],
        TE_raw,
        join_table_refined_statistics,
    )
