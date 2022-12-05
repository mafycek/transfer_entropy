#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import enum
import pickle
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
)

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


class FinanceDatabasePlugin(GenericDatabasePlugin):
    def __init__(self, database_url, database, username, password):
        super().__init__(database_url, database, username, password)
        print(f"SqlAlchemy version: {__version__ }")

        print(
            f"Connecting to postgresql://{self.username}:{self.password}@{self.database_url}/{self.database}"
        )
        self.engine = create_engine(
            f"postgresql+psycopg2://{self.username}:{self.password}@{self.database_url}/{self.database}",
            echo=True,
        )
        self.Session = sessionmaker(self.engine)

    def __del__(self):
        pass


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
        return f"<TimeseriesOHLCV(id={self.id}, timeseries_fk={self.timeseries_fk}, open={self.open}, high={self.high}, low={self.low}, close={self.close}, date={self.date})>"


def setup_database(financial_plugin):
    # creation of schema
    if not financial_plugin.engine.dialect.has_schema(
        financial_plugin.engine, db_schema
    ):
        financial_plugin.engine.execute(schema.CreateSchema(db_schema))

    # create tables
    Dataset.metadata.create_all(financial_plugin.engine)
    Parameters.metadata.create_all(financial_plugin.engine)
    TimeseriesBidAsk.metadata.create_all(financial_plugin.engine)
    TimeseriesPrice.metadata.create_all(financial_plugin.engine)
    TimeseriesOHLC.metadata.create_all(financial_plugin.engine)
    Timeseries.metadata.create_all(financial_plugin.engine)


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
        "localhost", "financial_dataset", "postgres", "qwerty"
    )

    setup_database(financial_plugin)

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

        size_of_chunk = 100000
        chunks = (len(raw_dataset)//size_of_chunk)+1
        for chunk_number in range(chunks):
            with financial_plugin.Session.begin() as session:
                new_records = []
                start_of_chunk = size_of_chunk*chunk_number
                for date, record in raw_dataset[start_of_chunk: start_of_chunk+size_of_chunk]:
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
