#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import typing
import traceback
from processing_datasets import *


def generate_epsilon_values(TE, epsilon_column, symbol_1, symbol2):
    epsilons = TE[epsilon_column].unique()
    set_symbols = set()
    for epsilon in epsilons:
        splitted_epsilon = epsilon.split("_")
        set_symbols.add(splitted_epsilon[0])
        set_symbols.add(splitted_epsilon[1])

    TE[symbol_1] = TE[epsilon_column].map(lambda value: value.split("_")[0])
    TE[symbol2] = TE[epsilon_column].map(lambda value: value.split("_")[1])

    symbols = list(set_symbols)
    return symbols


def generate_figures_from_data(
        TE,
        TE_cathegories,
        TE_cathegories_std,
        TE_nonswapped_columns,
        symbols,
        directory,
        name_of_title,
        figure_parameters: typing.Dict,
):
    for symbol in symbols:
        TE_selected = TE.loc[(TE["symbol_1"] == symbol) | (TE["symbol_2"] == symbol)]

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
            latex_title = (
                    figure_parameters["latex_title_size"] + f"""{{{pure_title}}}"""
            )
            latex_alpha_label = (
                    figure_parameters["latex_overview_label_size"] + r"$\alpha$"
            )
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
                (name_of_title, figure_parameters["latex_overview_label_size"]),
                errorbar_filename,
                figure_parameters["output_format"],
                dpi=figure_parameters["dpi"],
                fontsize=figure_parameters["fontsize"],
            )

            figures2d_TE_overview_alpha(
                TE_selected,
                item,
                latex_title,
                latex_alpha_label,
                (name_of_title, figure_parameters["latex_overview_label_size"]),
                standard_filename,
                figure_parameters["output_format"],
                dpi=figure_parameters["dpi"],
                fontsize=figure_parameters["fontsize"],
            )

            figures2d_TE_overview_alpha(
                TE_selected,
                item_std,
                latex_title,
                latex_alpha_label,
                (name_of_title, figure_parameters["latex_overview_label_size"]),
                std_filename,
                figure_parameters["output_format"],
                dpi=figure_parameters["dpi"],
                fontsize=figure_parameters["fontsize"],
            )

        for item in []:  # TE_nonswapped_columns:
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
                latex_title = (
                        figure_parameters["latex_title_size"] + f"""{{{pure_title}}}"""
                )
                latex_title_std = (
                        figure_parameters["latex_title_size"]
                        + f"""{{Standard deviation of {pure_title.lower()} }}"""
                )

                title_graph = {
                    "transfer_entropy": r"$\Huge\rm{Transfer\ entropy}$",
                    "conditional_information_transfer": r"$\Huge\rm{Conditional\ information\ transfer}$",
                }
                filename_direction = {True: "Y->X", False: "X->Y"}

                label_latex_std = (
                        figure_parameters["latex_label_size"] + f"""$\\sigma_{{{label}}}$"""
                )
                label = figure_parameters["latex_label_size"] + f"${label}$"
                latex_alpha_label = figure_parameters["latex_label_size"] + r"$\alpha$"
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
                    figure_parameters["output_format"],
                    dpi=figure_parameters["dpi"],
                    fontsize=figure_parameters["fontsize"],
                )
                figures2d_TE_alpha(
                    TE_selected,
                    tuple(complete_column_name_std),
                    latex_title_std,
                    latex_alpha_label,
                    label_latex_std,
                    std_filename,
                    figure_parameters["output_format"],
                    dpi=figure_parameters["dpi"],
                    fontsize=figure_parameters["fontsize"],
                )
                figures2d_TE_alpha_errorbar(
                    TE_selected,
                    item,
                    tuple(complete_column_name_std),
                    latex_title,
                    latex_alpha_label,
                    label,
                    errorbar_filename,
                    figure_parameters["output_format"],
                    dpi=figure_parameters["dpi"],
                    fontsize=figure_parameters["fontsize"],
                )
            except Exception as exc:
                print(f"Problem {exc} {item}")
                traceback.print_exc()


def generate_figures_from_neighborhood_statistics(
        TE_statistics, figure_parameters: typing.Dict, figure_generation
):
    selection = TE_statistics.xs(
        figure_generation["selector"], level=figure_generation["level"], axis=1
    )

    timeserie_names = set()
    RTE_types = set()
    time_sets = set()
    statistics_types = set()
    samples = set()
    for item in selection.columns:
        (
            RTE_type,
            history_first,
            future_first,
            history_second,
            type_of_statistics,
            shuffled_dataset,
            reversed_order,
            sample_number,
            timeserie_name,
        ) = item
        timeserie_names.add(timeserie_name)
        statistics_types.add(type_of_statistics)
        RTE_types.add(RTE_type)
        time_sets.add((history_first, future_first, history_second))
        samples.add(sample_number)

    subList = list(
        filter(lambda elem: "quantile" in elem or "median" in elem, statistics_types)
    )

    for timeserie_name in list(timeserie_names):
        timeserie_selection = selection.xs(timeserie_name, level="Sample", axis=1)

        for RTE_type in RTE_types:
            RTE_type_selection = timeserie_selection.xs(
                RTE_type, level="Name of variable", axis=1
            )
            for time_set in list(time_sets):
                reverse_directions = (
                    [True, False] if "balance" not in RTE_type else [False]
                )
                for actual_reverse_direction in reverse_directions:
                    shuffled_calculations = (
                        [True] if "effective" not in RTE_type else [False]
                    )
                    for actual_shuffled_calculation in shuffled_calculations:
                        history_first, future_first, history_second = time_set
                        time_set_selection_hf = RTE_type_selection.xs(
                            history_first, level="History first", axis=1
                        )
                        time_set_selection_ff = time_set_selection_hf.xs(
                            future_first, level="Future first", axis=1
                        )
                        time_set_selection_hs = time_set_selection_ff.xs(
                            history_second, level="History second", axis=1
                        )
                        time_set_selection_reverse = time_set_selection_hs.xs(
                            actual_reverse_direction, level="Reversed direction", axis=1
                        )
                        time_set_selection_complet = time_set_selection_reverse.xs(
                            actual_shuffled_calculation, level="Shuffled", axis=1
                        )

                        samples_set = set()
                        for item in time_set_selection_complet.columns:
                            samples_set.add(item[1])
                        samples = list(samples_set)

                        name_of_title = RTE_type.split("y_")[0]
                        symbol = timeserie_name
                        base_filename = name_of_title
                        shuffled_calculation = False

                        pure_title = (
                                name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title}}}"""
                        )
                        latex_alpha_label = (
                                figure_parameters["latex_overview_label_size"] + r"$\alpha$"
                        )
                        pure_title_std = (
                                "Standard deviation of "
                                + name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title_std = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title_std}}}"""
                        )

                        for figure_parameter in figure_generation["figure_generation"]:
                            standard_filename = figure_parameter["filename"](
                                base_filename,
                                symbol,
                                history_first,
                                future_first,
                                history_second,
                                actual_reverse_direction,
                                actual_shuffled_calculation,
                            )
                            figure_parameters = figure_parameter["profile"]

                            rows = len(samples)
                            figures2d_TE_overview_alpha_order_statistics(
                                time_set_selection_complet,
                                rows,
                                figure_parameters,
                                latex_title,
                                latex_alpha_label,
                                (
                                    name_of_title,
                                    figure_parameters["latex_overview_label_size"],
                                ),
                                standard_filename,
                                figure_parameters["output_format"],
                                dpi=figure_parameters["dpi"],
                                fontsize=figure_parameters["fontsize"],
                            )


def generate_figures_from_complete_statistics(
        TE_statistics, figure_parameters: typing.Dict, figure_generation
):
    selection = TE_statistics.xs(
        figure_generation["selector"], level=figure_generation["level"], axis=1
    )

    timeserie_names = set()
    RTE_types = set()
    time_sets = set()
    statistics_types = set()
    samples = set()
    for item in selection.columns:
        (
            RTE_type,
            history_first,
            future_first,
            history_second,
            type_of_statistics,
            shuffled_dataset,
            reversed_order,
            sample_number,
            timeserie_name,
        ) = item
        timeserie_names.add(timeserie_name)
        statistics_types.add(type_of_statistics)
        RTE_types.add(RTE_type)
        time_sets.add((history_first, future_first, history_second))
        samples.add(sample_number)

    subList = list(
        filter(lambda elem: "quantile" in elem or "median" in elem, statistics_types)
    )

    for timeserie_name in list(timeserie_names):
        timeserie_selection = selection.xs(timeserie_name, level="Sample", axis=1)

        for RTE_type in RTE_types:
            RTE_type_selection = timeserie_selection.xs(
                RTE_type, level="Name of variable", axis=1
            )
            for time_set in list(time_sets):
                reverse_directions = (
                    [True, False] if "balance" not in RTE_type else [False]
                )
                for actual_reverse_direction in reverse_directions:
                    shuffled_calculations = (
                        [True] if "effective" not in RTE_type else [False]
                    )
                    for actual_shuffled_calculation in shuffled_calculations:
                        history_first, future_first, history_second = time_set
                        time_set_selection_hf = RTE_type_selection.xs(
                            history_first, level="History first", axis=1
                        )
                        time_set_selection_ff = time_set_selection_hf.xs(
                            future_first, level="Future first", axis=1
                        )
                        time_set_selection_hs = time_set_selection_ff.xs(
                            history_second, level="History second", axis=1
                        )
                        time_set_selection_reverse = time_set_selection_hs.xs(
                            actual_reverse_direction, level="Reversed direction", axis=1
                        )
                        time_set_selection_complet = time_set_selection_reverse.xs(
                            actual_shuffled_calculation, level="Shuffled", axis=1
                        )

                        samples_set = set()
                        for item in time_set_selection_complet.columns:
                            samples_set.add(item[1])
                        samples = list(samples_set)

                        name_of_title = RTE_type.split("y_")[0]
                        symbol = timeserie_name
                        base_filename = name_of_title
                        shuffled_calculation = False

                        pure_title = (
                                name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title}}}"""
                        )
                        latex_alpha_label = (
                                figure_parameters["latex_overview_label_size"] + r"$\alpha$"
                        )
                        pure_title_std = (
                                "Standard deviation of "
                                + name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title_std = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title_std}}}"""
                        )

                        for figure_parameter in figure_generation["figure_generation"]:
                            standard_filename = figure_parameter["filename"](
                                base_filename,
                                symbol,
                                history_first,
                                future_first,
                                history_second,
                                actual_reverse_direction,
                                actual_shuffled_calculation,
                            )
                            figure_parameters = figure_parameter["profile"]

                            rows = 1

                            figures2d_TE_overview_alpha_order_statistics(
                                time_set_selection_complet,
                                rows,
                                figure_parameters,
                                latex_title,
                                latex_alpha_label,
                                (
                                    name_of_title,
                                    figure_parameters["latex_overview_label_size"],
                                ),
                                standard_filename,
                                figure_parameters["output_format"],
                                dpi=figure_parameters["dpi"],
                                fontsize=figure_parameters["fontsize"],
                            )

    print(selection)


def generate_overview_figures_from_averages(
        TE_statistics, figure_parameters, figure_generation
):
    selection = TE_statistics.xs(
        figure_generation["selector"], level=figure_generation["level"], axis=1
    )

    timeserie_names = set()
    RTE_types = set()
    time_sets = set()
    statistics_types = set()
    samples = set()

    for item in selection.columns:
        (
            RTE_type,
            history_first,
            future_first,
            history_second,
            type_of_statistics,
            shuffled_dataset,
            reversed_order,
            sample_number,
            timeserie_name,
        ) = item
        timeserie_names.add(timeserie_name)
        statistics_types.add(type_of_statistics)
        RTE_types.add(RTE_type)
        time_sets.add((history_first, future_first, history_second))

    print(
        list(timeserie_names), list(statistics_types), list(RTE_types), list(time_sets)
    )

    for timeserie_name in list(timeserie_names):
        timeserie_selection = selection.xs(timeserie_name, level="Sample", axis=1)

        for RTE_type in RTE_types:
            RTE_type_selection = timeserie_selection.xs(
                RTE_type, level="Name of variable", axis=1
            )
            for time_set in list(time_sets):
                reverse_directions = (
                    [True, False] if "balance" not in RTE_type else [False]
                )
                for actual_reverse_direction in reverse_directions:
                    shuffled_calculations = (
                        [False, True] if "effective" not in RTE_type else [False]
                    )
                    for actual_shuffled_calculation in shuffled_calculations:
                        history_first, future_first, history_second = time_set
                        time_set_selection_hf = RTE_type_selection.xs(
                            history_first, level="History first", axis=1
                        )
                        time_set_selection_ff = time_set_selection_hf.xs(
                            future_first, level="Future first", axis=1
                        )
                        time_set_selection_hs = time_set_selection_ff.xs(
                            history_second, level="History second", axis=1
                        )
                        time_set_selection_reverse = time_set_selection_hs.xs(
                            actual_reverse_direction, level="Reversed direction", axis=1
                        )
                        time_set_selection_complet = time_set_selection_reverse.xs(
                            actual_shuffled_calculation, level="Shuffled", axis=1
                        )

                        samples_set = set()
                        for item in time_set_selection_complet.columns:
                            samples_set.add(item[1])
                        samples = list(samples_set)

                        name_of_title = RTE_type.split("y_")[0]
                        symbol = timeserie_name
                        base_filename = name_of_title

                        # figure_parameters["directory"] + f"/{column_name}_{symbol}_{filename_direction[swapped_datasets]}" + ("_shuffled" if shuffled_calculation else "")
                        pure_title = (
                                name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title}}}"""
                        )
                        latex_alpha_label = (
                                figure_parameters["latex_overview_label_size"] + r"$\alpha$"
                        )
                        pure_title_std = (
                                "Standard deviation of "
                                + name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title_std = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title_std}}}"""
                        )

                        rows = (
                            (len(samples) if actual_shuffled_calculation else 1)
                            if isinstance(
                                time_set_selection_complet.iloc[0].iloc[0], float
                            )
                            else time_set_selection_complet.iloc[0].iloc[0].shape[0]
                        )

                        for figure_parameter in figure_generation["figure_generation"]:
                            standard_filename = figure_parameter["filename"](
                                base_filename,
                                symbol,
                                history_first,
                                future_first,
                                history_second,
                                actual_reverse_direction,
                                actual_shuffled_calculation,
                            )
                            figure_parameters = figure_parameter["profile"]

                            figures2d_TE_overview_alpha_order_statistics(
                                time_set_selection_complet,
                                rows,
                                figure_parameters,
                                latex_title,
                                latex_alpha_label,
                                (
                                    name_of_title,
                                    figure_parameters["latex_overview_label_size"],
                                ),
                                standard_filename,
                                figure_parameters["output_format"],
                                dpi=figure_parameters["dpi"],
                                fontsize=figure_parameters["fontsize"],
                            )


def generate_figures_aggregated_figures(TE_statistics, figure_parameters, figure_generation
                                        ):
    selection = TE_statistics.xs(
        figure_generation["selector"], level=figure_generation["level"], axis=1
    )

    timeserie_names = set()
    RTE_types = set()
    time_sets = set()
    statistics_types = set()
    samples = set()
    timeserie_names_categories = set()
    timeserie_names_groups = {}

    for item in selection.columns:
        (
            RTE_type,
            history_first,
            future_first,
            history_second,
            type_of_statistics,
            shuffled_dataset,
            reversed_order,
            sample_number,
            timeserie_name,
        ) = item
        timeserie_names.add(timeserie_name)
        timeserie_names_categories.add(timeserie_name.split("_")[0])
        timeserie_base_name = timeserie_name.split("_")[0]
        if timeserie_base_name not in timeserie_names_groups:
            timeserie_names_groups[timeserie_base_name] = set([timeserie_name])
        else:
            timeserie_names_groups[timeserie_base_name].add(timeserie_name)

        statistics_types.add(type_of_statistics)
        RTE_types.add(RTE_type)
        time_sets.add((history_first, future_first, history_second))

    print(
        list(timeserie_names), list(statistics_types), list(RTE_types), list(time_sets)
    )

    for timeseries_name, timeseries_groups in timeserie_names_groups.items():
        selected_columns = [item for item in selection.columns if item[8] in timeseries_groups]
        timeserie_selection = selection[[*selected_columns]]
        symbol = timeseries_name

        for RTE_type in RTE_types:
            RTE_type_selection = timeserie_selection.xs(
                RTE_type, level="Name of variable", axis=1
            )
            for time_set in list(time_sets):
                reverse_directions = (
                    [True, False] if "balance" not in RTE_type else [False]
                )
                for actual_reverse_direction in reverse_directions:
                    shuffled_calculations = (
                        [True, False] if "effective" not in RTE_type else [False]
                    )
                    for actual_shuffled_calculation in shuffled_calculations:
                        try:
                            history_first, future_first, history_second = time_set
                            time_set_selection_hf = RTE_type_selection.xs(
                                history_first, level="History first", axis=1
                            )
                            time_set_selection_ff = time_set_selection_hf.xs(
                                future_first, level="Future first", axis=1
                            )
                            time_set_selection_hs = time_set_selection_ff.xs(
                                history_second, level="History second", axis=1
                            )
                            time_set_selection_reverse = time_set_selection_hs.xs(
                                actual_reverse_direction, level="Reversed direction", axis=1
                            )
                            time_set_selection_complet = time_set_selection_reverse.xs(
                                actual_shuffled_calculation, level="Shuffled", axis=1
                            )
                        except KeyError as exc:
                            print(time_set_selection_reverse, actual_shuffled_calculation, RTE_type)

                        samples_set = set()
                        for item in time_set_selection_complet.columns:
                            samples_set.add(item[1])
                        samples = list(samples_set)

                        name_of_title = RTE_type.split("y_")[0]
                        base_filename = name_of_title

                        pure_title = (
                                name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        latex_title = (
                                figure_parameters["latex_title_size"]
                                + f"""{{{pure_title}}}"""
                        )
                        latex_alpha_label = (
                                figure_parameters["latex_overview_label_size"] + r"$\alpha$"
                        )
                        ylabel = figure_parameters["latex_overview_label_size"] + r"$" + translate_to_label_TE((
                            RTE_type,
                            history_first,
                            future_first,
                            history_second,
                            actual_reverse_direction,
                            actual_shuffled_calculation)) + r"$"
                        pure_title_std = (
                                "Standard deviation of "
                                + name_of_title.capitalize().replace("_", " ")
                                + f" of {symbol} for nearest neighbor"
                        )
                        rows = len(timeseries_groups)

                        for figure_parameter in figure_generation["figure_generation"]:
                            standard_filename = figure_parameter["filename"](
                                base_filename,
                                symbol,
                                history_first,
                                future_first,
                                history_second,
                                actual_reverse_direction,
                                actual_shuffled_calculation,
                            )
                            figure_generation_configuration = figure_parameter["profile"]

                            figures2d_TE_overview_alpha_timeseries(
                                time_set_selection_complet,
                                sorted(tuple(timeseries_groups)),
                                figure_generation_configuration,
                                latex_title,
                                latex_alpha_label,
                                ylabel,
                                standard_filename,
                                figure_parameters["output_format"],
                                dpi=figure_parameters["dpi"],
                                fontsize=figure_parameters["fontsize"],
                            )
