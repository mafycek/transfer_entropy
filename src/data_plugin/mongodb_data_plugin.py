#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Handler of MongoDB for financial datasets"""

import pickle
import os
import datetime
from urllib.parse import quote_plus
import pymongo
from dotenv import load_dotenv
from gridfs import GridFS
from pymongo.errors import ConnectionFailure

from data_plugin.generic_database_plugin import GenericDatabasePlugin
from data_plugin.finance_data_plugin import FinanceDataPlugin


class FinanceMongoDatabasePlugin(GenericDatabasePlugin):
    """
    Handler of finance dataset stored in MongoDB
    """

    dataset_collection_name = "dataset"
    conditional_information_transfer_name = "conditional_information_transfer"

    def __init__(
            self, database_url, database_table_name, username="admin", password="admin"
    ):
        mongo_uri = (
            f"mongodb://{quote_plus(username)}:{quote_plus(password)}@{database_url}/"
        )
        super().__init__(mongo_uri, database_table_name, username, password)
        self.client = None
        self.database = None
        self.gridfs_client = None

        self.reconnect()

        try:
            self.client.admin.command("ping")
        except ConnectionFailure:
            print("Server not available")

    def upload_to_gridfs(self, filename, data):
        """
        Uploads dataset in MongoDB with a filename
        :param filename:  of dataset
        :param data: Raw data
        :return: ObjectId of stared dataset
        """
        return self.gridfs_client.put(data, filename=filename)

    def download_from_gridfs(self, document_id):
        """
        Downloads dataset stored in MongoDB
        :param document_id: ObjectId of stored dataset
        :return: Raw data
        """
        return self.gridfs_client.get(document_id).read()

    def disconnect(self):
        if self.client:
            print(f"{self.client}")
            self.client.close()
            self.client = None
            self.database = None
            self.gridfs_client = None

            pid = os.getpid()
            timestamp = datetime.datetime.now().isoformat()
            print(f"PID:{pid} {timestamp} Disconnecting from Mongo DB")

    def reconnect(
            self, connect_timeout_ms=600000, socket_timeout_ms=300000, timeout_ms=150000
    ):
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Connecting to {self.database_url}"
        )
        self.client = pymongo.MongoClient(
            self.database_url,
            timeoutMS=timeout_ms,
            socketTimeoutMS=socket_timeout_ms,
            connectTimeoutMS=connect_timeout_ms,
        )
        self.database = self.client[self.database_table_name]
        self.gridfs_client = GridFS(self.database)

    def add_start_of_calculation(
            self, dataset_1_fk: int, dataset_2_fk: int, json_data={}
    ):
        # unimplemented method
        pass

    def __del__(self):
        # unimplemented method
        pass

    @staticmethod
    def CLISupport(parser):
        """
        Additional parameters for CLI interface
        :param parser: parser object used to handle  CLI
        """
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

    def insert_document(self, collection: str, document):
        collection = self.database[collection]
        post_id = collection.insert_one(document).inserted_id

        return post_id

    def get_all_documents(self, collection: str, querry):
        collection = self.database[collection]
        return collection.find(querry)

    def get_dataset_with_filename(self, file):
        dataset_cursor = self.get_all_documents(
            FinanceMongoDatabasePlugin.dataset_collection_name, {"file": file}
        )
        output_documents = []
        for dataset_document in dataset_cursor:
            output_documents.append(dataset_document)
        return output_documents

    def select_dataset_with_code(self, code, start_date=None, end_date=None):
        dataset_documents = self.get_all_documents(
            FinanceMongoDatabasePlugin.dataset_collection_name, {"code": code}
        )
        dataset_metadata_complete = {}
        dataset_raw_complete = {}
        for dataset_document in dataset_documents:
            if not (
                    (end_date and end_date < dataset_document["start_date_time"])
                    or (start_date and dataset_document["end_date_time"] < start_date)
            ):
                pickled_dataset_raw = self.download_from_gridfs(
                    dataset_document["dataset_id"]
                )
                dataset_raw = pickle.loads(pickled_dataset_raw)
                select_data_from_dataset(
                    dataset_raw, start_date, end_date, dataset_raw_complete
                )
                dataset_metadata_complete.update(dataset_document)

        return dataset_raw_complete, dataset_metadata_complete


def select_data_from_dataset(dataset_raw, start_date, end_date, dataset_raw_complete):
    for row_data in dataset_raw:
        if (
                (start_date and end_date and start_date <= row_data[0] <= end_date)
                or (start_date and not end_date and start_date <= row_data[0])
                or (not start_date and end_date and row_data[0] <= end_date)
                or (not start_date and not end_date)
        ):
            dataset_raw_complete[row_data[0]] = row_data[1]


def setup_database(database_plugin):
    list_of_collections = database_plugin.database.list_collection_names()
    if FinanceMongoDatabasePlugin.dataset_collection_name not in list_of_collections:
        database_plugin.database.create_collection(
            FinanceMongoDatabasePlugin.dataset_collection_name
        )


if __name__ == "__main__":
    load_dotenv()
    URL = os.getenv("URL")
    DATABASE = os.getenv("NOSQL_TABLE")
    USERNAME = os.getenv("NOSQL_USERNAME")
    PASSWORD = os.getenv("NOSQL_PASSWORD")
    mongo_uri = f"mongodb://{quote_plus(USERNAME)}:{quote_plus(PASSWORD)}@ashley.fjfi.cvut.cz:27017/"

    mongo_engine = FinanceMongoDatabasePlugin(URL, DATABASE, USERNAME, PASSWORD)

    # dataset = mongo_engine.select_dataset_with_code("APD")
    # results = mongo_engine.get_all_documents(
    #    "conditional_information_transfer", {"symbol": "APD_BAC"}
    # )

    # for result in results:
    #    decoded_result = mongo_engine.download_from_gridfs(result["document_id"])
    #    variable = pickle.loads(decoded_result, encoding="base64")
    #    print(variable)

    # post_structure = posts.find_one({"_id": ObjectId(str(post_id))})
    # print(post_structure)

    dataPlugin = FinanceDataPlugin(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data/1Q23_sp_stocks/"
        #"/home/hynek/work/skola/prog/python/transfer_entropy/data"
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

    # upload dataset to mongodb
    for dataset in datasets_generator:
        metadata = dataset["metadata"]
        raw_dataset = dataset["dataset"]
        pickled_raw_dataset = pickle.dumps(raw_dataset)

        metadata_in_database = mongo_engine.get_dataset_with_filename(metadata["file"])
        if len(metadata_in_database) == 0:
            gridfs_raw_dataset_id = mongo_engine.upload_to_gridfs(
                str(metadata["full_filename"]), pickled_raw_dataset
            )
            metadata["dataset_id"] = gridfs_raw_dataset_id
            metadata["start_date_time"] = raw_dataset[0][0]
            metadata["end_date_time"] = raw_dataset[-1][0]
            mongo_engine.insert_document(
                FinanceMongoDatabasePlugin.dataset_collection_name, metadata
            )
        else:
            print(f"Attepmpt to insert duplicated dataset with medatada: {metadata}")
