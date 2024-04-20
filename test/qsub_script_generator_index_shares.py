#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pprint
from dotenv import load_dotenv

from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin

if __name__ == "__main__":
    load_dotenv()
    url = os.getenv("URL")
    table = os.getenv("NOSQL_TABLE")
    username = os.getenv("NOSQL_USERNAME")
    password = os.getenv("NOSQL_PASSWORD")
    sql_username = os.getenv("SQL_USERNAME")
    sql_password = os.getenv("SQL_PASSWORD")
    sql_database = os.getenv("SQL_DATABASE")

    dbHandler = FinanceMongoDatabasePlugin(
        url,
        table,
        username,
        password,
    )

    collection_handler = dbHandler.database.get_collection(
        FinanceMongoDatabasePlugin.dataset_collection_name
    )
    dataset_cursor = collection_handler.find({})
    maximal_neighborhood = 50

    update_dicts = [
        {
            "dataset_X_selector": "6",
            "filename_template": "script_{}.sh",
            "postselection_X_future": None,
            "random_source": None,
            "comment": "Stock without external source",
        },
        {
            "dataset_X_selector": "6 8",
            "filename_template": "script_{}_i.sh",
            "postselection_X_future": "0",
            "random_source": None,
            "comment": "Stock with index excluding stock",
        },
        {
            "dataset_X_selector": "6",
            "filename_template": "script_{}_r.sh",
            "postselection_X_future": None,
            "random_source": 0.1,
            "comment": "Stock with random source",
        },
        {
            "dataset_X_selector": "6 8",
            "filename_template": "script_{}_ri.sh",
            "postselection_X_future": "0",
            "random_source": 0.1,
            "comment": "Stock with index excluding stock and randomness",
        },
    ]

    selected_codes = [
        "united_health_group",
        "visa",
        "amazon",
        "johnson",
        "tesla",
        "microsoft",
        "apple",
        "exxon",
    ]

    codes = set()
    dataset_metadata = {}
    for dataset in dataset_cursor:
        if dataset["type"] == "stock_with_index":  # aggregated, tick
            codes.add(dataset["code"])
            dataset_metadata[dataset["code"]] = dataset

    codes = codes.intersection(set(selected_codes))
    list_of_codes = sorted(list(codes))
    pprint.pprint(f"{list_of_codes}")

    for update_dict in update_dicts:
        number_of_samples = 30
        history_first = "0, 0 1 , 0 1 2 3"
        history_first = "0 1 2 , 0 1 2 3 , 0 1 2 3 4"
        history_second = "0 , 0 1 , 0 1 2 3"
        history_second = "0 1 2 , 0 1 2 3 , 0 1 2 3 4"
        alpha_params = "0.1 3.0 117"
        dataset_X_selector = update_dict["dataset_X_selector"]
        dataset_Y_selector = "7"
        postselection_X_future = update_dict["postselection_X_future"]
        postselection_X_history = None
        postselection_Y_history = None
        random_source = update_dict["random_source"]
        comment = update_dict["comment"]

        directory = os.getcwd() + "/"
        for count, base_code in enumerate(list_of_codes):
            #        for code in list_of_codes[count + 1: len(list_of_codes)]:
            with open(directory + update_dict["filename_template"].format(base_code), "wt") as fh:
                print(
                    f'./calculations/transfer_entropy.sif --database=True --database_nosql_url={url} --database_sql_url={url} --database_sql_username={sql_username} --database_sql_password={sql_password} --database_sql_table={sql_database} --database_nosql_username={username} --database_nosql_password={password} --database_nosql_table={table} --dataset_1_code={base_code} --dataset_2_code={base_code} --dataset_1_selector {dataset_X_selector} --dataset_2_selector {dataset_Y_selector} --alpha_params {alpha_params} --history_first {history_first} --history_second {history_second} --number_of_samples {number_of_samples} --maximal_neighborhood {maximal_neighborhood} --comment "{comment}"' + (
                        "" if not random_source else f" --random_source {random_source}") + (
                        "" if not postselection_X_future else f" --postselection_X_future {postselection_X_future}") + (
                        "" if not postselection_X_history else f" --postselection_X_history {postselection_X_history}") + (
                        "" if not postselection_Y_history else f" --postselection_Y_history {postselection_Y_history}"),
                    file=fh,
                )
