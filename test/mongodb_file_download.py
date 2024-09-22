#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import os
import glob
from datetime import datetime
from dotenv import load_dotenv

from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin

if __name__ == "__main__":
    load_dotenv()
    url = os.getenv("URL")
    table = os.getenv("NOSQL_TABLE")
    username = os.getenv("NOSQL_USERNAME")
    password = os.getenv("NOSQL_PASSWORD")
    download_folder = "test/random"

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

    symbols = [
        "tesla_tesla",
        "amazon_amazon",
        "microsoft_microsoft",
        "exxon_exxon",
        "apple_apple",
        "johnson_johnson",
        "visa_visa",
        "united_united",
    ]

    for symbol in symbols:
        result = dbHandler.database.command("dbstats")
        collection = dbHandler.database[
            FinanceMongoDatabasePlugin.conditional_information_transfer_name
        ]
        document_cursor = collection.find(
            {
                "symbol": symbol,
                "end_timestamp": {
                    "$gt": datetime.fromisoformat("2024-05-30T00:00:01.879Z")
                },
            }
        )
        for metadata_transfer_entropy in document_cursor:
            random_source = (
                False if metadata_transfer_entropy["random_source"] is None else True
            )
            index = (
                False
                if metadata_transfer_entropy["postselection_X_future"] is None
                else True
            )
            if random_source or index:
                appendix = (
                        "_"
                        + ("" if index is False else "i")
                        + ("" if random_source is False else "r")
                )
            else:
                appendix = ""

            document_id = metadata_transfer_entropy["document_id"]
            symbol_in_file = symbol
            random_source_strength = "" if len(metadata_transfer_entropy['random_source']) == 0 else str(
                metadata_transfer_entropy['random_source'][0])
            filename = f"{download_folder}/conditional_information_transfer-{symbol_in_file}{appendix}-{random_source_strength}.bin"
            files = glob.glob(filename)
            if len(files) == 0:
                dataset_pickled = dbHandler.download_from_gridfs(document_id)
                dataset = pickle.loads(dataset_pickled)

                with open(filename, "wb") as fh:
                    fh.write(dataset_pickled)

    collection_handler = dbHandler.database[
        FinanceMongoDatabasePlugin.conditional_information_transfer_name
    ]
    code = "BATRK"
    cursor_find = collection_handler.find(
        {"end_timestamp": {"$gt": datetime.fromisoformat("2023-11-08T12:00:01.879Z")}}
    )

    for cursor in cursor_find:
        document_id = cursor["document_id"]
        symbol = cursor["symbol"]
        document_pickle = dbHandler.download_from_gridfs(document_id)
        document = pickle.loads(document_pickle)

        with open(f"Conditional_information_transfer-{symbol}.bin", "wb") as fh:
            fh.write(document_pickle)

    if False:
        # remove data from database
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
