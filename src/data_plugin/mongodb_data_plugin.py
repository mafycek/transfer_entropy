#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import base64
import pymongo
import pickle
import datetime
import pprint
from urllib.parse import quote_plus
from gridfs import GridFS
from pymongo.errors import ConnectionFailure
from bson.objectid import ObjectId

from data_plugin.generic_database_plugin import GenericDatabasePlugin
from data_plugin.finance_data_plugin import FinanceDataPlugin


class FinanceMongoDatabasePlugin(GenericDatabasePlugin):
    dataset_collection_name = "dataset"
    conditional_information_transfer_name = "conditional_information_transfer"

    def __init__(self, database_url, database, username, password):
        username = "admin"
        password = "admin"
        mongo_uri = (
            f"mongodb://{quote_plus(username)}:{quote_plus(password)}@{database_url}/"
        )

        super().__init__(mongo_uri, database, username, password)

        print(f"Connecting to {mongo_uri}")
        self.client = pymongo.MongoClient(mongo_uri)
        self.database = self.client[database]
        self.gridfs_client = GridFS(self.database)

        try:
            self.client.admin.command('ping')
        except ConnectionFailure:
            print("Server not available")

    def upload_to_gridfs(self, filename, data):
        return self.gridfs_client.put(data, filename=filename)

    def download_from_gridfs(self, document_id):
        return self.gridfs_client.get(document_id).read()

    def __del__(self):
        pass

    @staticmethod
    def CLISupport(parser):
        parser.add_argument(
            "--database_nosql_url",
            type=str,
            default="",
            help="NoSQL database URL",
        )
        parser.add_argument(
            "--database_nosql_username",
            type=str,
            default="mongo",
            help="NoSQL database username",
        )
        parser.add_argument(
            "--database_nosql_password",
            type=str,
            default="mongo",
            help="NoSQL database password",
        )
        parser.add_argument(
            "--database_nosql_table",
            type=str,
            default="test",
            help="NoSQL database table",
        )

    def select_dataset_with_code(self, code):
        pass

    def insert_document(self, collection: str, document):
        collection = self.database[collection]
        post_id = collection.insert_one(document).inserted_id

        return post_id

    def get_all_documents(self, collection: str, querry):
        collection = self.database[collection]
        return collection.find(querry)


def setup_database(database_plugin):
    list_of_collections = database_plugin.database.list_collection_names()
    if FinanceMongoDatabasePlugin.dataset_collection_name not in list_of_collections:
        database_plugin.database.create_collection(FinanceMongoDatabasePlugin.dataset_collection_name)


if __name__ == "__main__":
    url = "ashley.fjfi.cvut.cz:27017"
    database = "test"
    username = "mongo"
    password = "mongo"
    mongo_uri = f"mongodb://{quote_plus(username)}:{quote_plus(password)}@ashley.fjfi.cvut.cz:27017/"

    mongo_engine = FinanceMongoDatabasePlugin(url
                                              , database, username, password
                                              )

    results = mongo_engine.get_all_documents("conditional_information_transfer", {"symbol": "APD_BAC"})

    for result in results:
        decoded_result = mongo_engine.download_from_gridfs(result["document_id"])
        variable = pickle.loads(decoded_result, encoding="base64")
        print(variable)

    # post_structure = posts.find_one({"_id": ObjectId(str(post_id))})
    # print(post_structure)

    dataPlugin = FinanceDataPlugin(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data"
    )
    old_data = False
    if old_data:
        with open(dataPlugin.file_pickled, "rb") as fh:
            dataset, metadata = pickle.load(fh)
    else:
        datasets_generator = dataPlugin.new_load_datasets()

    setup_database(mongo_engine)
    # updated_row = mongo_engine.add_start_of_calculation(1, 2, {"ABC": "CDE"})
    # mongo_engine.update_finish_of_calculation(updated_row)
    # data = mongo_engine.select_dataset_with_code("ADUS")

    for dataset in datasets_generator:
        metadata = dataset["metadata"]
        raw_dataset = dataset["dataset"]
        pickled_raw_dataset = pickle.dumps(raw_dataset)

        gridfs_raw_dataset_id = mongo_engine.upload_to_gridfs(str(metadata["full_filename"]), pickled_raw_dataset)
        metadata["dataset_id"] = gridfs_raw_dataset_id
        dataset_keys = raw_dataset.keys()
        metadata["start_date_time"] = min(dataset_keys)
        metadata["end_date_time"] = max(dataset_keys)
        mongo_engine.insert_document(FinanceMongoDatabasePlugin.dataset_collection_name, metadata)
