import os
import pprint

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

    dataset_collection = dbHandler.database.get_collection(
        FinanceMongoDatabasePlugin.dataset_collection_name
    )
    dataset_cursor = dataset_collection.find({})

    codes = set()
    for dataset in dataset_cursor:
        if dataset["type"] == "aggregated":
            codes.add(dataset["code"])

    list_of_codes = list(codes)
    pprint.pprint(f"{list_of_codes}")

    base_code = list_of_codes[0]
    directory = os.getcwd() + "/../data/scripts/"
    for code in list_of_codes[1:100]:
        with open(directory + f"script_{base_code}_{code}.sh", "wt") as fh:
            print(
                f"./calculations/transfer_entropy.sif --database=True --database_nosql_url=ashley.fjfi.cvut.cz --database_sql_url=ashley.fjfi.cvut.cz --database_sql_username=postgresrenyi --database_sql_password=postgresrenyi --database_sql_table=postgres --database_nosql_table=test --dataset_1_code={base_code} --dataset_2_code={code} --dataset_1_selector 4 5 6 7 --dataset_2_selector 4 5 6 7 --alpha_params 0.1 3.0 80 --history_first 0 , 0 1 , 0 1 2 3  --history_second 0 , 0 1 , 0 1 2 3 --maximal_neighborhood 80",
                file=fh,
            )
