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

            figure_generation_NN = {
                "selector": "Avaraging over nearest neighbours starting from 5",
                "level": "Remark",
                "figure_generation": [
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: (
                                figure_parameters["directory"]
                                + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                + ("_reversed" if actual_reverse_direction else "")
                                + ("_shuffled" if actual_shuffled_calculation else "")
                                + "_quantiles"
                                + "_NN"
                        ),
                        "profile": [
                            {
                                "selector": ["0.1 quantile", "0.9 quantile"],
                                "color": "moccasin",
                                "style": "quantile band",
                                "label": "0.1-0.9 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["0.25 quantile", "0.75 quantile"],
                                "color": "lightcoral",
                                "style": "quantile band",
                                "label": "0.25-0.75 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["median"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "median",
                                "aggregation": "sample",
                            },
                        ],
                    },
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_std"
                                                                                                  + "_NN",
                        "profile": [
                            {
                                "selector": ["mean", "std"],
                                "color": "moccasin",
                                "style": "yerrorbar",
                                "label": "errorbar",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["mean"],
                                "color": "green",
                                "style": "line",
                                "label": "mean",
                                "aggregation": "sample",
                            },
                        ],
                    },
                ],
            }
            # generate_overview_figures_from_averages(
            #    table, figure_parameters, figure_generation_NN
            # )

            figure_generation_SA = {
                "selector": "Avaraging over nearest neighbours starting from 5",
                "level": "Remark",
                "figure_generation": [
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_quantiles_SA",
                        "profile": [
                            {
                                "selector": ["quantile 0.1", "quantile 0.9"],
                                "color": "moccasin",
                                "style": "quantile band",
                                "label": "0.1-0.9 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["quantile 0.25", "quantile 0.75"],
                                "color": "lightcoral",
                                "style": "quantile band",
                                "label": "0.25-0.75 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["median"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "median",
                                "aggregation": "sample",
                            },
                        ],
                    },
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_std_SA",
                        "profile": [
                            {
                                "selector": ["mean", "std"],
                                "color": "moccasin",
                                "style": "yerrorbar",
                                "label": "errorbar",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["mean"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "mean",
                                "aggregation": "sample",
                            },
                        ],
                    },
                ],
            }

            # generate_figures_from_neighborhood_statistics(
            #    table, figure_parameters, figure_generation_SA
            # )

            figure_generation_NNSA = {
                "selector": f"Avaraging over samples and neighbours k",
                "level": "Remark",
                "figure_generation": [
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_quantiles_NNSA",
                        "profile": [
                            {
                                "selector": ["quantile 0.1", "quantile 0.9"],
                                "color": "moccasin",
                                "style": "quantile band",
                                "label": "0.1-0.9 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["quantile 0.25", "quantile 0.75"],
                                "color": "lightcoral",
                                "style": "quantile band",
                                "label": "0.25-0.75 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["median"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "median",
                                "aggregation": "sample",
                            },
                        ],
                    },
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_std_NNSA",
                        "profile": [
                            {
                                "selector": ["mean", "std"],
                                "color": "moccasin",
                                "style": "yerrorbar",
                                "label": "errorbar",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["mean"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "mean",
                                "aggregation": "sample",
                            },
                        ],
                    },
                ],
            }
            # generate_figures_from_complete_statistics(
            #    table, figure_parameters, figure_generation_NNSA
            # )

            figure_generation_aggregation = {
                "selector": f"Avaraging over samples and neighbours k",
                "level": "Remark",
                "figure_generation": [
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_quantiles_aggregation",
                        "profile": [
                            {
                                "selector": ["quantile 0.1", "quantile 0.9"],
                                "color": "moccasin",
                                "style": "quantile band",
                                "label": "0.1-0.9 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["quantile 0.25", "quantile 0.75"],
                                "color": "lightcoral",
                                "style": "quantile band",
                                "label": "0.25-0.75 quantile band",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["median"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "median",
                                "aggregation": "sample",
                            },
                        ],
                    },
                    {
                        "filename": lambda base_filename, symbol, history_first, future_first, history_second,
                                           actual_reverse_direction, actual_shuffled_calculation: figure_parameters[
                                                                                                      "directory"
                                                                                                  ]
                                                                                                  + f"/{base_filename}_{symbol}_{history_first}_{future_first}_{history_second}"
                                                                                                  + (
                                                                                                      "_reversed" if actual_reverse_direction else "")
                                                                                                  + (
                                                                                                      "_shuffled" if actual_shuffled_calculation else "")
                                                                                                  + "_std_aggregation",
                        "profile": [
                            {
                                "selector": ["mean", "std"],
                                "color": "moccasin",
                                "style": "yerrorbar",
                                "label": "errorbar",
                                "aggregation": "sample",
                            },
                            {
                                "selector": ["mean"],
                                "color": "greenyellow",
                                "style": "line",
                                "label": "mean",
                                "aggregation": "sample",
                            },
                        ],
                    },
                ],
            }
            generate_figures_aggregated_figures(table, figure_parameters, figure_generation_aggregation)
