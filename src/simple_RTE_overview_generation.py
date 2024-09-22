#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from figure_generation import *

if __name__ == "__main__":
    directory = "test/random"
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
    configuration_figures = {
        "merged_dataset": f"merged_{name_of_title}.bin",
        "collect_dataset": f"collect_{name_of_title}.bin",
        "title": name_of_title,
        "folders": [
            {"name": "test/random/r=0.0", "title": "r=0"},
            {"name": "test/random/r=0.001", "title": "r=10^{-3}"},
            {"name": "test/random/r=0.000001", "title": "r=10^{-6}"},
            {"name": "test/random/r=0.0000001", "title": "r=10^{-7}"},
            {"name": "test/random/r=0.00000001", "title": "r=10^{-8}"},
            {"name": "test/random/r=0.000000001", "title": "r=10^{-9}"},
            {"name": "test/random/r=1e-12", "title": "r=10^{-12}"},
            {"name": "test/random/r=1e-15", "title": "r=10^{-15}"},
            {"name": "test/random/r=1e-20", "title": "r=10^{-20}"},
            {"name": "test/random/r=1e-23", "title": "r=10^{-23}"},
        ],
        "profile": [
            {
                "selector": ["median"],
                "color": "greenyellow",
                "style": "line",
                "label": "median",
                "aggregation": "sample",
            },
        ],
    }

    generate_overview_figures_within_category(
        figure_parameters, configuration_figures
    )
