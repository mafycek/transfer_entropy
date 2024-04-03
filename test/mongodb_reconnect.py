#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import os
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

    dbHandler.disconnect()
    dbHandler.reconnect()
    collection_handler = dbHandler.database[FinanceMongoDatabasePlugin.conditional_information_transfer_name]
    code = "BATRK"
    cursor_find = collection_handler.find({"dataset_1_code": code})

    for item in cursor_find:
        print(f"{item}")

    dbHandler.disconnect()
