#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import sys
from pathlib import Path

if __name__ == "__main__":
    directories = [
        # "conditional_information_transfer",
        # "conditional_information_transfer_3",
        # "financial_transfer_entropy"
        # , "financial_transfer_entropy_2", "financial_transfer_entropy_1", "financial_transfer_entropy"
        # "indices/short_memory_random=0.1",
        "test"
    ]

    for directory in directories:
        name_of_title = "conditional_information_transfer"
        processed_dataset = f"{directory}/refined_{name_of_title}-*.bin"
        merged_dataset = f"{directory}/merged_{name_of_title}.bin"
        merged_files = glob.glob(merged_dataset)
        if len(merged_files) == 0:
            files = glob.glob(processed_dataset)
            print(f"Number of inputs: {len(files)}, files: {files}")

            files_to_merge = []
            for file in files:
                table = pd.read_pickle(Path(file))
                table.columns.names = names = [
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
                files_to_merge.append(table)

            refined_frame = pd.concat(files_to_merge, axis=1)
            refined_frame.to_pickle(merged_dataset)
        else:
            print(f"Merged file exists")
