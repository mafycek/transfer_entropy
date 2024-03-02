#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pprint
import os
from datetime import datetime
from dotenv import load_dotenv

from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin

if __name__ == "__main__":
    load_dotenv()
    url = os.getenv("URL")
    table = os.getenv("NOSQL_TABLE")
    username = os.getenv("NOSQL_USERNAME")
    password = os.getenv("NOSQL_PASSWORD")

    dbHandler = FinanceMongoDatabasePlugin(
        url,
        table,
        username,
        password,
    )

    result = dbHandler.database.command("dbstats")

    collection_handler = dbHandler.database[
        FinanceMongoDatabasePlugin.dataset_collection_name
    ]
    cursor_find = collection_handler.find(
        {"type": "tick"}
    )

    code_set = set()
    dateset_dates = set()
    datasets = {}
    for cursor in cursor_find:
        document_id = cursor["dataset_id"]
        code = cursor["code"]
        start_date_time = cursor["start_date_time"]
        dateset_date = f"{start_date_time.year}-{start_date_time.month}"
        code_set.add(code)
        dateset_dates.add(dateset_date)
        if code not in datasets:
            datasets[code] = {}

        datasets[code][dateset_date] = cursor

    fx_datesets = list(code_set)
    print(list(code_set))
    print(list(dateset_dates))

    for code, record in datasets.items():
        print(f"{code} ", end="")
        for dataset_date in record.keys():
            print(f"{dataset_date} ", end="")
        print("")

    for dataset_date in list(dateset_dates):
        print(f"{dataset_date} ", end="")
        for code, record in datasets.items():
            if dataset_date in record:
                print(f"{code} ", end="")
        print("")
