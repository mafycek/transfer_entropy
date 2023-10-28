#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import datetime
import os
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

    history_first = "0 , 0 1 , 0 1 2 3"
    history_first = "0 1"
    history_second = "0 , 0 1 , 0 1 2 3"
    history_second = "0 1"
    alpha_params = "0.1 3.0 117"
    dataset_1_selector = "2 3"
    dataset_2_selector = "2 3"

    dbHandler = FinanceMongoDatabasePlugin(
        url,
        table,
        username,
        password,
    )

    collection_handler = dbHandler.database.get_collection(
        FinanceMongoDatabasePlugin.dataset_collection_name
    )
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

    date_code = "2022-3"
    list_of_codes = []
    for code, record in datasets.items():
        if date_code in record:
            list_of_codes.append(code)

    maximal_neighborhood = 20
    start_date = datetime.datetime(int(date_code.split("-")[0]), int(date_code.split("-")[1]), 1, 0, 0, 0)
    end_date = datetime.datetime(int(date_code.split("-")[0]), int(date_code.split("-")[1]), 7, 0, 0, 0)

    base_code = list_of_codes[0]
    directory = os.getcwd() + "/"
    for count, base_code in enumerate(list_of_codes):
        for code in list_of_codes[count + 1: len(list_of_codes)]:
            with open(directory + f"script_{base_code}_{code}.sh", "wt") as fh:
                print(
                    f"./calculations/transfer_entropy.sif --database=True --database_nosql_url={url} --database_sql_url={url} --database_sql_username={sql_username} --database_sql_password={sql_password} --database_sql_table={sql_database} --database_nosql_username={username} --database_nosql_password={password} --database_nosql_table={table} --dataset_1_code={base_code} --dataset_2_code={code} --dataset_1_selector {dataset_1_selector} --dataset_2_selector {dataset_2_selector} --alpha_params {alpha_params} --history_first {history_first} --history_second {history_second} --maximal_neighborhood {maximal_neighborhood} --start_date {start_date.isoformat()} --end_date {end_date.isoformat()} --comment {date_code}",
                    file=fh,
                )
