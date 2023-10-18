#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
from src.data_plugin.mongodb_data_plugin import FinanceMongoDatabasePlugin

if __name__ == "__main__":
    url = "ashley.fjfi.cvut.cz"
    table = "test"
    username = "mongo"
    password = "mongo"
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

    symbol = "BATRK_XONE"
    result = dbHandler.database.command("dbstats")
    collection = dbHandler.database[FinanceMongoDatabasePlugin.conditional_information_transfer_name]
    document_cursor = collection.find({"symbol": symbol})
    for metadata_transfer_entropy in document_cursor:
        document_id = metadata_transfer_entropy["document_id"]

        dataset_pickled = dbHandler.download_from_gridfs(document_id)
        dataset = pickle.loads(dataset_pickled)

        with open(f"Conditional_information_transfer-{symbol}.bin", "wb") as fh:
            fh.write(dataset_pickled)

    if False:
        dbHandler.DATABASE.get_collection(FinanceMongoDatabasePlugin.conditional_information_transfer_name).delete_many(
            {})
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
