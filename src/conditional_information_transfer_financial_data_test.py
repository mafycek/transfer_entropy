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
from src.data_plugin.generic_data_plugin import GenericDataPlugin
from src.data_plugin.finance_data_plugin import FinanceDataPlugin
from src.data_plugin.finance_database_plugin import (
    FinanceDatabasePlugin,
    CalculationStatusType,
)
from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin
from transfer_entropy import renyi_conditional_information_transfer
from src.data_plugin.sample_generator import preparation_dataset_for_transfer_entropy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates conditional information transfer for financial datasets."
    )
    parser.add_argument(
        "--database",
        type=bool,
        default=False,
        help="Switch to use db",
    )

    # inserts support for SQL database arguments
    FinanceDatabasePlugin.CLISupport(parser)
    # inserts support for NoSQL database arguments
    FinanceMongoDatabasePlugin.CLISupport(parser)

    parser.add_argument(
        "--directory",
        type=str,
        default="conditional_information_transfer",
        help="Folder to export results",
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
        default=[1],
        help="History to take into account",
    )
    parser.add_argument(
        "--future_first",
        metavar="XXX",
        type=str,
        nargs="+",
        default=[1],
        help="Future to take into account",
    )
    parser.add_argument(
        "--history_second",
        metavar="XXX",
        type=str,
        nargs="+",
        default=[1],
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
        default=False,
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
        "--dataset_1_code", help="First dataset to load", default="sp500"
    )
    parser.add_argument(
        "--dataset_2_code", help="Second dataset to load", default="visa"
    )
    parser.add_argument(
        "--dataset_1_selector",
        metavar="XXX",
        type=int,
        nargs="+",
        help="Selector of columns for first dataset",
        default=[0],
    )
    parser.add_argument(
        "--dataset_2_selector",
        metavar="XXX",
        type=int,
        nargs="+",
        help="Selector of columns for second dataset",
        default=[0],
    )
    parser.add_argument(
        "--postselection_X_future",
        metavar="XXX",
        type=int,
        default=None,
        nargs="+",
        help="Postselector of future X",
    )
    parser.add_argument(
        "--postselection_X_history",
        metavar="XXX",
        type=int,
        default=None,
        nargs="+",
        help="Postselector of history X",
    )
    parser.add_argument(
        "--postselection_Y_history",
        metavar="XXX",
        type=int,
        default=None,
        nargs="+",
        help="Postselector of history Y",
    )
    parser.add_argument(
        "--random_source",
        metavar="XXX",
        type=float,
        default=None,
        nargs="+",
        help="Random sources added to X",
    )
    parser.add_argument(
        "--alpha_params",
        metavar="XXX",
        type=str,
        nargs="+",
        help="Calculation for alphas, min, max, number between",
    )
    parser.add_argument(
        "--start_date",
        metavar="XXX",
        type=str,
        help="Start date of joint timeseries",
        default=None,
    )
    parser.add_argument(
        "--end_date",
        metavar="XXX",
        type=str,
        help="End date of joint timeseries",
        default=None,
    )

    parser.add_argument(
        "--comment",
        type=str,
        help="Comment that is added to database",
        default=None,
    )

    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help="Performs test run without calculation",
    )

    args = parser.parse_args()

    postselection_X_future = args.postselection_X_future
    postselection_X_history = args.postselection_X_history
    postselection_Y_history = args.postselection_Y_history

    args.start_datetime = (
        datetime.datetime.fromisoformat(args.start_date) if args.start_date else None
    )
    args.end_datetime = (
        datetime.datetime.fromisoformat(args.end_date) if args.end_date else None
    )

    if args.test:
        test_run = True
    else:
        test_run = False


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

    future_history_joint = zip(future_firsts, histories_firsts, histories_seconds)

    if args.alpha_params:
        alpha_params = (
            float(args.alpha_params[0]),
            float(args.alpha_params[1]),
            int(args.alpha_params[2]),
        )
    else:
        alpha_params = (0.1, 1.9, 37)

    # create alphas that are been calculated
    alphas = np.round(
        np.linspace(alpha_params[0], alpha_params[1], alpha_params[2], endpoint=True), 3
    )
    print(
        f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} History parameters: {histories_firsts} {histories_seconds} {future_firsts}"
    )
    print(f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Alphas: {alphas}")

    arbitrary_precision = args.arbitrary_precision
    arbitrary_precision_decimal_numbers = args.arbitrary_precision_decimal_numbers
    dataset_1_selector = args.dataset_1_selector
    dataset_2_selector = args.dataset_2_selector
    random_source = args.random_source

    print(f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Load datasets")
    if args.database:
        # dataset_handler = FinanceDatabasePlugin(
        #    args.database_sql_url,
        #    args.database_sql_table,
        #    args.database_sql_username,
        #    args.database_sql_password,
        # )
        dataset_handler = FinanceMongoDatabasePlugin(
            args.database_nosql_url,
            args.database_nosql_table,
            args.database_nosql_username,
            args.database_nosql_password,
        )

    else:
        dataset_handler = FinanceDataPlugin(os.getcwd() + "/../data/1Q23_sp_stocks")

    # load static dataset

    dataset_handler.load_datasets()

    # selection of the dataset 1
    # dataset_handler
    dataset1, metadata1 = dataset_handler.select_dataset_with_code(
        args.dataset_1_code, args.start_datetime, args.end_datetime
    )
    # dataset_handler
    dataset2, metadata2 = dataset_handler.select_dataset_with_code(
        args.dataset_2_code, args.start_datetime, args.end_datetime
    )
    dataset_handler.disconnect()

    if not dataset1 or not dataset2:
        raise RuntimeError(
            f"Any of input datasets are empty. The selection code {args.dataset_1_code} or {args.dataset_2_code} is wrong or time period {args.start_datetime}--{args.end_datetime} is incorrect."
        )

    joint_dataset = FinanceDataPlugin.time_join_dataset(
        dataset1, dataset2, dataset_1_selector, dataset_2_selector
    )

    print(
        f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Datasets {metadata1['code']} and {metadata2['code']} loaded and merged together"
    )
    symbol = f"{metadata1['code']}_{metadata2['code']}"

    # preparation of
    parameters = {
        "histories_firsts": histories_firsts,
        "future_firsts": future_firsts,
        "histories_seconds": histories_seconds,
        "arbitrary_precision": arbitrary_precision,
        "arbitrary_precision_decimal_numbers": arbitrary_precision_decimal_numbers,
        "dataset_1_selector": dataset_1_selector,
        "dataset_2_selector": dataset_2_selector,
        "dataset_1_code": args.dataset_1_code,
        "dataset_2_code": args.dataset_2_code,
        "maximal_neighborhood": args.maximal_neighborhood,
        "alphas": alphas.tolist(),
        "symbol": symbol,
        "collection": FinanceMongoDatabasePlugin.conditional_information_transfer_name,
        "start_timestamp": datetime.datetime.now(),
        "comment": args.comment,
        "postselection_X_history": postselection_X_history,
        "postselection_Y_history": postselection_Y_history,
        "postselection_X_future": postselection_X_future,
        "random_source": random_source,
    }

    if args.database and "id" in metadata1 and "id" in metadata2:
        # insert record to database
        calculation_id = dataset_handler.add_start_of_calculation(
            metadata1["id"], metadata2["id"], parameters
        )
        parameters["calculation_id"] = calculation_id

    # create structure for results
    results = {}

    for future_firsts, histories_firsts, histories_seconds in future_history_joint:
        for swap_datasets in [False, True]:
            # loop over shuffling
            for shuffle_dataset in [False, True]:
                # prepare dataset that is being processed

                (
                    marginal_solution_1,
                    marginal_solution_2,
                ) = GenericDataPlugin.prepare_dataset(
                    datasets=joint_dataset,
                    swap_datasets=swap_datasets,
                    shuffle_dataset=shuffle_dataset,
                    selection1=len(dataset_1_selector),
                    selection2=len(dataset_2_selector),
                )

                for future_first in future_firsts:
                    for histories_first in histories_firsts:
                        for histories_second in histories_seconds:
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} History first: {histories_first}, future first: {future_first}, history second: {histories_second} and symbol: {symbol} is processed",
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
                                "postselection_y_fut": postselection_X_future,
                                "postselection_z_hist": postselection_Y_history,
                                "postselection_y_hist": postselection_X_history,
                                "random_source": random_source,
                            }

                            # prepare samples to be used to calculate transfer entropy
                            t0 = time.process_time()
                            (
                                y_fut,
                                y_hist,
                                z_hist,
                            ) = preparation_dataset_for_transfer_entropy(
                                marginal_solution_1, marginal_solution_2, **configuration
                            )
                            t1 = time.process_time()
                            duration = t1 - t0
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Preparation of datasets [s]: {duration}",
                                flush=True,
                            )

                            # create range of indices that will be used for calculation
                            indices_to_use = list(range(1, args.maximal_neighborhood + 1))
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
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Transfer entropy for history first: {histories_first}, future first: {future_first}, history second: {histories_second} and symbol: {symbol} shuffling: {shuffle_dataset} and swapped dataset {swap_datasets} is calculated",
                                flush=True,
                            )
                            t0 = time.process_time()
                            if not test_run:
                                transfer_entropy = renyi_conditional_information_transfer(
                                    y_fut, y_hist, z_hist, **configuration
                                )
                            else:
                                transfer_entropy = 0

                            t1 = time.process_time()
                            duration = t1 - t0
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Duration of calculation of transfer entropy [s]: {duration}",
                                flush=True,
                            )

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
                                    symbol,
                                    f"{string_histories_first}_{string_future_first}",
                                    string_histories_second,
                                    shuffle_dataset,
                                    swap_datasets,
                                )
                            ] = transfer_entropy
                            print(
                                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} * Transfer entropy calculation for history first: {histories_first}, future first: {future_first}, history second: {histories_second} and symbol: {symbol}, shuffling: {shuffle_dataset} and swapped dataset {swap_datasets} is finished",
                                flush=True,
                            )

    if args.database and not test_run:
        dataset_handler.reconnect()

        # save to database
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Save to database",
            flush=True,
        )

        pickled_result = pickle.dumps(results, -1)
        document_id = dataset_handler.upload_to_gridfs(
            f"Conditional_information_transfer-{symbol}.bin", pickled_result
        )
        parameters["document_id"] = document_id
        parameters["format"] = "pickle"
        parameters["end_timestamp"] = datetime.datetime.now()
        mongo_id = dataset_handler.insert_document(
            FinanceMongoDatabasePlugin.conditional_information_transfer_name, parameters
        )
        parameters["document_id"] = str(mongo_id)

        # insert record to database

        if "calculation_id" in parameters:
            del parameters["_id"]
            dataset_handler.update_state_of_calculation(
                parameters["calculation_id"], CalculationStatusType.FINISHED, parameters
            )
        dataset_handler.disconnect()

    else:
        # save result structure to the file
        path = Path(f"{args.directory}/Conditional_information_transfer-{symbol}.bin")
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Save to file {path}",
            flush=True,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fb:
            pickle.dump(results, fb)

    print(
        f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Calculation finished",
        flush=True,
    )
