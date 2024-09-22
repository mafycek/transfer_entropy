#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import argparse

import processing_datasets

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates conditional information transfer for financial datasets."
    )
    parser.add_argument(
        "--symbol",
        metavar="XXX",
        type=str,
        nargs="+",
        help="Name to process",
    )

    args = parser.parse_args()

    directories = [
        # "conditional_information_transfer",
        # "conditional_information_transfer_3",
        # "financial_transfer_entropy"
        # , "financial_transfer_entropy_2", "financial_transfer_entropy_1", "financial_transfer_entropy"
        # "indices/short_memory_random=0.1",
        # "test/random/r=0.001",
        # "test/random/r=0.000001",
        # "test/random/r=0.0",
        # "test/random/r=0.000000001",
        # "test/random/r=0.00000001",
        # "test/random/r=0.0000001",
        # "test"
        # "test/random/r=1e-12",
        # "test/random/r=1e-15",
        # "test/random/r=1e-20",
        "test/random/r=1e-23",
    ]

    for directory in directories:
        for symbol in args.symbol:
            name_of_title = "conditional_information_transfer"
            processed_dataset = f"{directory}/refined_{name_of_title}-{symbol}.bin"
            input_dataset = f"{directory}/{name_of_title}-{symbol}.bin"
            files = glob.glob(processed_dataset)

            if len(files) == 0:
                processing_datasets.refined_process_dataset(
                    input_dataset, processed_dataset
                )
