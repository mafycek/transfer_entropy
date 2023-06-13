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

    dataset_collection = dbHandler.database.get_collection(
        FinanceMongoDatabasePlugin.dataset_collection_name
    )
    dataset_cursor = dataset_collection.find({})
    maximal_neighborhood = 100

    codes = set()
    dataset_metadata = {}
    for dataset in dataset_cursor:
        if dataset["type"] == "aggregated":  # aggregated, tick
            codes.add(dataset["code"])
            dataset_metadata[dataset["code"]] = dataset

    list_of_codes = list(codes)
    pprint.pprint(f"{list_of_codes}")

    base_code = list_of_codes[0]
    directory = os.getcwd() + "/../data/scripts/"
    for code in list_of_codes[1:500]:
        with open(directory + f"script_{base_code}_{code}.sh", "wt") as fh:
            print(
                f"./calculations/transfer_entropy.sif --database=True --database_nosql_url={url} --database_sql_url={url} --database_sql_username={sql_username} --database_sql_password={sql_password} --database_sql_table={sql_database} --database_nosql_username={username} --database_nosql_password={password} --database_nosql_table=t{table} --dataset_1_code={base_code} --dataset_2_code={code} --dataset_1_selector 4 5 6 7 --dataset_2_selector 4 5 6 7 --alpha_params 0.1 3.0 117 --history_first 0 , 0 1 , 0 1 2 3  --history_second 0 , 0 1 , 0 1 2 3 --maximal_neighborhood {maximal_neighborhood}",
                file=fh,
            )
