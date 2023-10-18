#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import enum
import json
import pickle
import pprint
from sqlalchemy.orm import sessionmaker, Session, declarative_base
from sqlalchemy import (
    Column,
    Integer,
    String,
    DateTime,
    Float,
    ForeignKey,
    select,
    schema,
    Enum,
    types,
    MetaData,
    PrimaryKeyConstraint,
    __version__,
    create_engine,
    inspect,
    insert,
    update,
)
from sqlalchemy.dialects.postgresql import JSONB

from data_plugin.generic_database_plugin import GenericDatabasePlugin
from data_plugin.finance_data_plugin import FinanceDataPlugin


db_schema = "TEST_SCHEMA"
metadata_obj = MetaData(schema=db_schema)

Base = declarative_base(metadata=metadata_obj)


class DatasetType(enum.Enum):
    TIMESERIES_ASK_BID = 1  # tick currency
    TIMESERIES_PRICE = 2  # price
    TIMESERIES_OHLC = 3  # aggregated

    @staticmethod
    def to_data_storage(dataset):
        table_lookup = {
            DatasetType.TIMESERIES_PRICE: TimeseriesPrice,
            DatasetType.TIMESERIES_ASK_BID: TimeseriesBidAsk,
            DatasetType.TIMESERIES_OHLC: TimeseriesOHLC,
        }
        return table_lookup[dataset]

    @staticmethod
    def from_str(label):
        if label == "aggregated":
            return DatasetType.TIMESERIES_OHLC
        elif label == "tick":
            return DatasetType.TIMESERIES_ASK_BID
        elif label == "shortened":
            return DatasetType.TIMESERIES_PRICE
        else:
            return None


class CalculationStatusType(enum.Enum):
    STARTED = 1  # tick currency
    FINISHED = 2  # price
    FAILED = 3  # aggregated


class FinanceDatabasePlugin(GenericDatabasePlugin):
    def __init__(self, database_url, database, username, password):
        super().__init__(database_url, database, username, password)
        print(f"SqlAlchemy version: {__version__}")

        print(
            f"Connecting to postgresql://{self.username}:{self.password}@{self.database_url}/{self.database_table_name}"
        )
        self.engine = create_engine(
            f"postgresql+psycopg2://{self.username}:{self.password}@{self.database_url}/{self.database_table_name}",
            echo=True,
        )
        self.Session = sessionmaker(self.engine)

    def __del__(self):
        self.engine.dispose()

    @staticmethod
    def CLISupport(parser):
        parser.add_argument(
            "--database_sql_url",
            type=str,
            default="",
            help="SQL database URL",
        )
        parser.add_argument(
            "--database_sql_username",
            type=str,
            default="postgres",
            help="SQL database username",
        )
        parser.add_argument(
            "--database_sql_password",
            type=str,
            default="postgres",
            help="SQL database password",
        )
        parser.add_argument(
            "--database_sql_table",
            type=str,
            default="postgres",
            help="SQL database table",
        )

    def select_dataset_with_code(self, code, start_date, end_date):
        # implement method
        sql_statement = select(Dataset).where(Dataset.code == code)
        with Session(self.engine) as session:
            dataset_infos = session.execute(sql_statement).all()
            pprint.pprint(dataset_infos)

        if len(dataset_infos) > 1:
            print(f"More than one record found. Taking the first one.")
        dataset_info = dataset_infos[0][0]

        dataset_metadata = dataset_info.__dict__

        with Session(self.engine) as session:
            timeseries_storage = DatasetType.to_data_storage(dataset_metadata["type"])
            sql_statement = (
                select(timeseries_storage)
                .where(timeseries_storage.dataset_fk == dataset_metadata["id"])
                .order_by(timeseries_storage.date)
            )
            dataset_data_rows = session.execute(sql_statement).all()

        dataset_data = {}
        price_multiplicator = dataset_metadata["price_multiplicator"]
        for row in dataset_data_rows:
            row_dict = row[0].__dict__
            if (
                    (start_date and end_date and start_date <= row_dict["date"] <= end_date)
                    or (start_date and not end_date and start_date <= row_dict["date"])
                    or (not start_date and end_date and row_dict["date"] <= end_date)
                    or (not start_date and not end_date)
            ):
                if dataset_metadata["type"] == DatasetType.TIMESERIES_OHLC:
                    dataset_data[row_dict["date"]] = (
                        row_dict["open"] / price_multiplicator,
                        row_dict["high"] / price_multiplicator,
                        row_dict["low"] / price_multiplicator,
                        row_dict["close"] / price_multiplicator,
                    )
                elif dataset_metadata["type"] == DatasetType.TIMESERIES_PRICE:
                    dataset_data[row_dict["date"]] = row_dict["value"] / price_multiplicator
                elif dataset_metadata["type"] == DatasetType.TIMESERIES_ASK_BID:
                    dataset_data[row_dict["date"]] = (
                        row_dict["bid"] / price_multiplicator,
                        row_dict["ask"] / price_multiplicator,
                    )

        return dataset_data, dataset_metadata

    def add_start_of_calculation(self, dataset_1_fk: int, dataset_2_fk: int, json_data={}):
        sql_statement = insert(RenyiEntropyCalculation).values(
            dataset_1_fk=dataset_1_fk,
            dataset_2_fk=dataset_2_fk,
            start_date=datetime.datetime.now(),
            status=CalculationStatusType.STARTED,
            parameters=json_data,
        )
        with Session(self.engine) as session:
            dataset_infos = session.execute(sql_statement)
            session.commit()
        return dataset_infos.inserted_primary_key[0]

    def update_state_of_calculation(
            self, id: int, status=CalculationStatusType.FINISHED, json_data={}
    ):
        sql_statement = (
            update(RenyiEntropyCalculation)
            .where(RenyiEntropyCalculation.id == id)
            .values(
                end_date=datetime.datetime.now(),
                status=status,
                parameters=json_data,
            )
        )
        with Session(self.engine) as session:
            dataset_infos = session.execute(sql_statement)
            session.commit()


class Dataset(Base):
    __tablename__ = "DATASET"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    description = Column(String)
    type = Column(Enum(DatasetType))
    header = Column(String)
    code = Column(String)
    file = Column(String)
    full_filename = Column(String)
    price_multiplicator = Column(Float)

    def __repr__(self):
        return f"<Dataset(id={self.id}, name={self.name}, description={self.description}, type={self.type}, header={self.header}, code={self.code}, file={self.file}, full_filename={self.full_filename}, price_multiplicator={self.price_multiplicator})>"


class Parameters(Base):
    __tablename__ = "PARAMETER"
    __table_args__ = {"schema": db_schema}

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    name = Column(String)
    value = Column(String)

    def __repr__(self):
        return f"<Dataset(id={self.id}, dataset_fk={self.dataset_fk}, name={self.name}, value={self.value})>"


class Timeseries(Base):
    __tablename__ = "TIMESERIES"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    value = Column(Integer)
    date = Column(DateTime, primary_key=True)

    def __repr__(self):
        return f"<Timeseries(id={self.id}, dataset_fk={self.dataset_fk}, value={self.value}, date={self.date})>"


class TimeseriesBidAsk(Base):
    __tablename__ = "TIMESERIES_BID_ASK"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    bid = Column(Integer)
    ask = Column(Integer)
    date = Column(DateTime, primary_key=True)

    def __repr__(self):
        return f"<TimeseriesBidAsk(id={self.id}, dataset_fk={self.dataset_fk}, bid={self.bid}, ask={self.ask}, date={self.date})>"


class TimeseriesPrice(Base):
    __tablename__ = "TIMESERIES_PRICE"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    value = Column(Integer)
    date = Column(DateTime, primary_key=True)

    def __repr__(self):
        return f"<TimeseriesPrice(id={self.id}, dataset_fk={self.dataset_fk}, value={self.value}, date={self.date})>"


class TimeseriesOHLC(Base):
    __tablename__ = "TIMESERIES_OHLCV"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    open = Column(Integer)
    high = Column(Integer)
    low = Column(Integer)
    close = Column(Integer)
    date = Column(DateTime, primary_key=True)

    def __repr__(self):
        return f"<TimeseriesOHLCV(id={self.id}, dataset_fk={self.dataset_fk}, open={self.open}, high={self.high}, low={self.low}, close={self.close}, date={self.date})>"


class RenyiEntropyCalculation(Base):
    __tablename__ = "RENYI_ENTROPY_CALCULATION"
    __table_args__ = (PrimaryKeyConstraint("id"),)
    id = Column(Integer, primary_key=True, autoincrement=True)
    start_date = Column(DateTime)
    end_date = Column(DateTime)
    status = Column(Enum(CalculationStatusType))
    dataset_1_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    dataset_2_fk = Column(Integer, ForeignKey(f"{db_schema}.DATASET.id"))
    parameters = Column(JSONB)

    def __repr__(self):
        return f"<RenyiEntropyCalculation(id={self.id}, start_date={self.start_date}, end_date={self.end_date}, status={self.status}, dataset_1_fk={self.dataset_1_fk}, dataset_2_fk={self.dataset_2_fk}, parameters={self.parameters})>"


def setup_database(financial_plugin):
    # creation of schema
    with financial_plugin.engine.connect() as conn:
        if not inspect(financial_plugin.engine).has_schema(db_schema):
            financial_plugin.engine.execute(schema.CreateSchema(db_schema))

        # create tables
        if not inspect(financial_plugin.engine).has_table(Dataset.__tablename__):
            Dataset.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(Parameters.__tablename__):
            Parameters.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(
                TimeseriesBidAsk.__tablename__
        ):
            TimeseriesBidAsk.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(
                TimeseriesPrice.__tablename__
        ):
            TimeseriesPrice.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(TimeseriesOHLC.__tablename__):
            TimeseriesOHLC.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(Timeseries.__tablename__):
            Timeseries.metadata.create_all(financial_plugin.engine)
        if not inspect(financial_plugin.engine).has_table(
                RenyiEntropyCalculation.__tablename__
        ):
            RenyiEntropyCalculation.metadata.create_all(financial_plugin.engine)


if __name__ == "__main__":
    dataPlugin = FinanceDataPlugin(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data"
    )
    old_data = False
    if old_data:
        with open(dataPlugin.file_pickled, "rb") as fh:
            dataset, metadata = pickle.load(fh)
    else:
        datasets_generator = dataPlugin.new_load_datasets()

    financial_plugin = FinanceDatabasePlugin(
        "ashley.fjfi.cvut.cz", "postgres", "postgresrenyi", "postgresrenyi"
    )
    setup_database(financial_plugin)
    updated_row = financial_plugin.add_start_of_calculation(1, 2, {"ABC": "CDE"})
    financial_plugin.update_finish_of_calculation(updated_row)
    data = financial_plugin.select_dataset_with_code("ADUS")

    for dataset in datasets_generator:
        metadata = dataset["metadata"]
        raw_dataset = dataset["dataset"]
        dataset_id = None
        with financial_plugin.Session.begin() as session:
            record = session.scalars(
                select(Dataset).where(Dataset.code == metadata["code"])
            ).all()
            data_to_add = []
            if not record:
                # add new record to the table of datasets
                data_to_add.append(
                    Dataset(
                        code=metadata["code"],
                        description=metadata["code"],
                        file=metadata["file"],
                        full_filename=str(metadata["full_filename"]),
                        header=" ".join(metadata["header"])
                        if "header" in metadata
                        else metadata["format"],
                        name="",
                        price_multiplicator=metadata["price_multiplicator"],
                        type=DatasetType.from_str(metadata["type"]),
                    )
                )
                session.add_all(data_to_add)
                session.commit()
            else:
                session.query(Dataset).filter(Dataset.id == record[0].id).update(
                    {"file": Dataset.file + f" {metadata['file']}"}
                )
                session.commit()

        with financial_plugin.Session.begin() as session:
            record = session.scalars(
                select(Dataset).where(Dataset.code == metadata["code"])
            ).all()
            dataset_id = record[0].id

        size_of_chunk = 7000000
        chunks = (len(raw_dataset) // size_of_chunk) + 1
        for chunk_number in range(chunks):
            with financial_plugin.Session.begin() as session:
                new_records = []
                start_of_chunk = size_of_chunk * chunk_number
                for date, record in raw_dataset[
                    start_of_chunk : start_of_chunk + size_of_chunk
                ]:
                    if metadata["type"] == "aggregated":
                        new_records.append(
                            TimeseriesOHLC(
                                dataset_fk=dataset_id,
                                open=record[0],
                                high=record[1],
                                low=record[2],
                                close=record[3],
                                date=date,
                            )
                        )
                    elif metadata["type"] == "shortened":
                        new_records.append(
                            TimeseriesPrice(
                                dataset_fk=dataset_id, value=record[0], date=date
                            )
                        )
                    elif metadata["type"] == "tick":
                        new_records.append(
                            TimeseriesBidAsk(
                                dataset_fk=dataset_id,
                                bid=record[0],
                                ask=record[1],
                                date=date,
                            )
                        )
                session.add_all(new_records)
                session.commit()
                del new_records
        del raw_dataset
        del metadata
        del dataset
