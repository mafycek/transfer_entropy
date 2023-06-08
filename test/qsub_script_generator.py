#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pprint

from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin

if __name__ == "__main__":
    url = "ashley.fjfi.cvut.cz"
    table = "test"
    username = "admin"
    password = "admin"

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
    for code in list_of_codes[1:200]:
        with open(directory + f"script_{base_code}_{code}.sh", "wt") as fh:
            print(
                f"./calculations/transfer_entropy.sif --database=True --database_nosql_url=ashley.fjfi.cvut.cz --database_sql_url=ashley.fjfi.cvut.cz --database_sql_username=postgresrenyi --database_sql_password=postgresrenyi --database_sql_table=postgres --database_nosql_username=admin --database_nosql_password=admin --database_nosql_table=test --dataset_1_code={base_code} --dataset_2_code={code} --dataset_1_selector 4 5 6 7 --dataset_2_selector 4 5 6 7 --alpha_params 0.1 3.0 117 --history_first 0 , 0 1 , 0 1 2 3  --history_second 0 , 0 1 , 0 1 2 3 --maximal_neighborhood {maximal_neighborhood}",
                file=fh,
            )
