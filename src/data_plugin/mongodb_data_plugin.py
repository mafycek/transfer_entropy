import pymongo
import datetime
import pprint
from urllib.parse import quote_plus
from bson.objectid import ObjectId

from data_plugin.generic_database_plugin import GenericDatabasePlugin


class FinanceMongoDatabasePlugin(GenericDatabasePlugin):
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

    def __del__(self):
        self.client.close()
        self.database = None

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


if __name__ == "__main__":
    username = "mongo"
    password = "mongo"
    mongo_uri = f"mongodb://{quote_plus(username)}:{quote_plus(password)}@ashley.fjfi.cvut.cz:27017/"

    mongo_engine = FinanceMongoDatabasePlugin(
        "ashley.fjfi.cvut.cz:27017", "test", username, password
    )
    mongo_engine.insert_document(
        "test-calculation",
        {"metadata": "meta", "blob": b"asdasdceaaer3452345jk;jasdjhfljasgfadf"},
    )

    post = {
        "author": "Mike",
        "text": "My first blog post!",
        "tags": ["mongodb", "python", "pymongo"],
        "date": datetime.datetime.utcnow(),
        "document": "ABC",
        "document2": "DCE",
    }
    db = mongo_engine.client["test"]
    posts = db.posts
    post_id = posts.insert_one(post).inserted_id
    print(str(post_id))

    post_structure = posts.find_one({"_id": ObjectId(str(post_id))})
    print(post_structure)
    post_structure = posts.find_one({"author": "Mike"})
    print(post_structure)
    pprint.pprint(posts.find_one())
