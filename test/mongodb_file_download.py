#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
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

    # results = dbHandler.get_all_documents("conditional_information_transfer", {"symbol": "APD_BAC"})

    # for result in results:
    #    decoded_result = dbHandler.download_from_gridfs(result["document_id"])
    #    variable = pickle.loads(decoded_result, encoding="base64")
    #    print(variable)

    result = dbHandler.database.command("dbstats")

    collection_handler = dbHandler.database[
        FinanceMongoDatabasePlugin.conditional_information_transfer_name
    ]
    code = "BATRK"
    cursor_find = collection_handler.find(
        {"end_timestamp": {"$gt": datetime.fromisoformat("2023-07-30T00:00:01.879Z")}}
    )

    for cursor in cursor_find:
        document_id = cursor["document_id"]
        symbol = cursor["symbol"]
        document_pickle = dbHandler.download_from_gridfs(document_id)
        document = pickle.loads(document_pickle)

        with open(f"Conditional_information_transfer-{symbol}.bin", "wb") as fh:
            fh.write(document_pickle)

    if False:
        dbHandler.DATABASE.get_collection(
            FinanceMongoDatabasePlugin.conditional_information_transfer_name
        ).delete_many({})
        results = dbHandler.get_all_documents("dataset", {})
        id = dbHandler.gridfs_client.find_one({})
        dbHandler.gridfs_client.delete(id._id)
        for result in results:
            if dbHandler.gridfs_client.exists(result["dataset_id"]):
                dbHandler.gridfs_client.delete(result["dataset_id"])

        results = dbHandler.gridfs_client.find({})
        for result in results:
            if dbHandler.gridfs_client.exists(result["_id"]):
                dbHandler.gridfs_client.delete(result["_id"])
