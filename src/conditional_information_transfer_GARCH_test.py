#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datetime
import os
import pickle
import time
from pathlib import Path

import numpy as np

from cli_helpers import process_CLI_arguments
from src.data_plugin.GARCH_data_plugin import GARCHDataPlugin
from transfer_entropy import renyi_conditional_information_transfer
from src.data_plugin.sample_generator import preparation_dataset_for_transfer_entropy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates conditional information transfer for GARCH dataset."
    )
    parser.add_argument(
        "--directory",
        type=str,
        default="conditional_information_transfer",
        help="Folder to export results",
    )
    parser.add_argument(
        "--skip",
        metavar="XXX",
        type=int,
        default=2000,
        help="Skipped results of integration",
    )
    parser.add_argument(
        "--blockwise",
        metavar="XXX",
        type=int,
        default=0,
        help="Blockwise calculation of distances to prevent excessive memory usage",
    )
    parser.add_argument(
        "--history_first",
        metavar="XXX",
        type=str,
        nargs="+",
        help="History to take into account",
    )
    parser.add_argument(
        "--future_first",
        metavar="XXX",
        type=str,
        nargs="+",
        help="Future to take into account",
    )
    parser.add_argument(
        "--history_second",
        metavar="XXX",
        type=str,
        nargs="+",
        help="History to take into account",
    )
    parser.add_argument(
        "--maximal_neighborhood",
        metavar="XXX",
        type=int,
        default=2,
        help="Maximal neighborhood",
    )
    parser.add_argument(
        "--arbitrary_precision",
        action="store_true",
        help="Calculates the main part in arbitrary precision",
    )
    parser.add_argument(
        "--arbitrary_precision_decimal_numbers",
        metavar="XXX",
        type=int,
        default=100,
        help="Sets number saved in arbitrary precision arithmetic",
    )
    parser.add_argument(
        "--dataset", action="store_true", help="Use dataset", default=False
    )
    parser.add_argument(
        "--dataset_range", metavar="XXX-YYY", type=str, help="Dataset with range"
    )
    parser.add_argument(
        "--alpha_params",
        metavar="XXX",
        type=str,
        nargs="+",
        help="Calculation for alphas, min, max, number between",
    )

    args = parser.parse_args()
    # print(args.epsilon, flush=True)

    if args.history_first:
        histories_firsts = process_CLI_arguments(args.history_first)
    else:
        histories_firsts = range(2, 25)

    if args.future_first:
        future_firsts = process_CLI_arguments(args.future_first)
    else:
        future_firsts = range(2, 25)

    if args.history_second:
        histories_seconds = process_CLI_arguments(args.history_second)
    else:
        histories_seconds = range(2, 25)

    if args.alpha_params:
        alpha_params = (
            float(args.alpha_params[0]),
            float(args.alpha_params[1]),
            int(args.alpha_params[2]),
        )
    else:
        alpha_params = (0.1, 1.9, 37)

    dataset_handler = GARCHDataPlugin(os.getcwd() + "/data")

    # create alphas that are been calculated
    alphas = np.round(
        np.linspace(alpha_params[0], alpha_params[1], alpha_params[2], endpoint=True), 3
    )

    print(
        f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} History parameters: {histories_firsts} {histories_seconds} {future_firsts}"
    )
    print(f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Alphas: {alphas}")

    # create alphas that are been calculated
    alphas = np.round(np.linspace(0.1, 1.9, 37, endpoint=True), 3)

    arbitrary_precision = args.arbitrary_precision
    arbitrary_precision_decimal_numbers = args.arbitrary_precision_decimal_numbers

    # load static dataset
    if args.dataset:
        datasets, epsilons = dataset_handler.load_static_dataset(args)
    else:
        datasets = None

    print(
        f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Calculated epsilons: {epsilons}"
    )

    # loop over different realizations for various epsilon
    for index_epsilon, epsilon in enumerate(epsilons):
        # create structure for results
        results = {}

        for swap_datasets in [False, True]:
            # loop over shuffling
            for shuffle_dataset in [True, False]:
                configuration_of_integration = {
                    "cache": True,
                    "epsilon": epsilon,
                    "arbitrary_precision": arbitrary_precision,
                    "arbitrary_precision_decimal_numbers": arbitrary_precision_decimal_numbers,
                }

                # prepare dataset that is been processed
                (
                    marginal_solution_1,
                    marginal_solution_2,
                ) = GARCHDataPlugin.prepare_dataset(
                    args,
                    index_epsilon=index_epsilon,
                    datasets=datasets,
                    swap_datasets=swap_datasets,
                    shuffle_dataset=shuffle_dataset,
                    configuration_of_integration=configuration_of_integration,
                )

                for future_first in future_firsts:
                    for histories_first in histories_firsts:
                        for histories_second in histories_seconds:
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} History first: {histories_first}, future first: {future_first}, history second: {histories_second} and epsilon: {epsilon} is processed",
                                flush=True,
                            )

                            # preparation of the configuration dictionary
                            # additional +1 is there for separation
                            configuration = {
                                "transpose": True,
                                "blockwise": args.blockwise,
                                "history_index_x": histories_first,
                                "history_index_y": histories_second,
                                "future_index_x": future_first,
                            }

                            # prepare samples to be used to calculate transfer entropy
                            t0 = time.process_time()
                            (
                                y_fut,
                                y_hist,
                                z_hist,
                            ) = preparation_dataset_for_transfer_entropy(
                                marginal_solution_2,
                                marginal_solution_1,
                                **configuration,
                            )
                            t1 = time.process_time()
                            duration = t1 - t0
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Preparation of datasets [s]: {duration}",
                                flush=True,
                            )

                            # create range of indices that will be used for calculation
                            indices_to_use = list(
                                range(1, args.maximal_neighborhood + 1)
                            )
                            configuration = {
                                "transpose": True,
                                "axis_to_join": 0,
                                "method": "LeonenkoProzanto",
                                "alphas": alphas,
                                "enhanced_calculation": True,
                                "indices_to_use": indices_to_use,
                                "arbitrary_precision": arbitrary_precision,
                                "arbitrary_precision_decimal_numbers": arbitrary_precision_decimal_numbers,
                            }

                            # calculation of transfer entropy
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Transfer entropy for history first: {histories_first}, future first: {future_first}, history second: {histories_second} and epsilon: {epsilon} shuffling: {shuffle_dataset} and swapped dataset {swap_datasets} is calculated",
                                flush=True,
                            )
                            t0 = time.process_time()
                            transfer_entropy = renyi_conditional_information_transfer(
                                y_fut, y_hist, z_hist, **configuration
                            )
                            t1 = time.process_time()
                            duration = t1 - t0
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Duration of calculation of transfer entropy [s]: {duration}",
                                flush=True,
                            )
                            # print(f" * Transfer Renyi entropy with {history} {epsilon}: {transfer_entropy}", flush=True)

                            # store transfer entropy to the result structure
                            string_histories_first = ",".join(
                                str(x) for x in histories_first
                            )
                            string_future_first = ",".join(str(x) for x in future_first)
                            string_histories_second = ",".join(
                                str(x) for x in histories_second
                            )
                            results[
                                (
                                    epsilon,
                                    f"{string_histories_first}_{string_future_first}",
                                    string_histories_second,
                                    shuffle_dataset,
                                    swap_datasets,
                                )
                            ] = transfer_entropy
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Transfer entropy calculation for history first: {histories_first}, future first: {future_first}, history second: {histories_second} and epsilon: {epsilon}, shuffling: {shuffle_dataset} and swapped dataset {swap_datasets} is finished",
                                flush=True,
                            )

        # save result structure to the file
        path = Path(f"{args.directory}/Conditional_information_transfer-{epsilon}.bin")
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Save to file {path}",
            flush=True,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fb:
            pickle.dump(results, fb)
