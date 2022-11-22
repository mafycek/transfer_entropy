#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
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
    Enum,
    PrimaryKeyConstraint,
    __version__,
    create_engine,
)

from data_plugin.generic_database_plugin import GenericDatabasePlugin

Base = declarative_base()


class Dataset(Base):
    __tablename__ = "dataset"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    description = Column(String)
    type = Column(String)
    header = Column(String)
    code = Column(String)
    file = Column(String)
    full_filename = Column(String)
    price_multiplicator = Column(Float)

    def __repr__(self):
        return f"<Dataset(id={self.id}, name={self.name}, description={self.description}, type={self.type}, header={self.header}, code={self.code}, file={self.file}, full_filename={self.full_filename}, price_multiplicator={self.price_multiplicator})>"


class Parameters(Base):
    __tablename__ = "parameter"

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey("dataset.id"))
    name = Column(String)
    value = Column(String)

    def __repr__(self):
        return f"<Dataset(id={self.id}, dataset_fk={self.dataset_fk}, name={self.name}, value={self.value})>"


class Timeseries(Base):
    __tablename__ = "timeseries"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey("dataset.id"))
    value = Column(Integer)
    date = Column(DateTime)

    def __repr__(self):
        return f"<Timeseries(id={self.id}, dataset_fk={self.dataset_fk}, value={self.value}, date={self.date})>"


class TimeseriesBidAsk(Base):
    __tablename__ = "timeseries_bid_ask"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey("dataset.id"))
    bid = Column(Integer)
    ask = Column(Integer)
    date = Column(DateTime)

    def __repr__(self):
        return f"<TimeseriesBidAsk(id={self.id}, dataset_fk={self.dataset_fk}, bid={self.bid}, ask={self.ask}, date={self.date})>"


class TimeseriesPrice(Base):
    __tablename__ = "timeseries_price"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey("dataset.id"))
    value = Column(Integer)
    date = Column(DateTime)

    def __repr__(self):
        return f"<TimeseriesPrice(id={self.id}, dataset_fk={self.dataset_fk}, value={self.value}, date={self.date})>"


class TimeseriesOHLC(Base):
    __tablename__ = "timeseries_OHLCV"
    __table_args__ = (PrimaryKeyConstraint("id", "date"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataset_fk = Column(Integer, ForeignKey("dataset.id"))
    open = Column(Integer)
    high = Column(Integer)
    low = Column(Integer)
    close = Column(Integer)
    date = Column(DateTime)

    def __repr__(self):
        return f"<TimeseriesOHLCV(id={self.id}, timeseries_fk={self.timeseries_fk}, open={self.open}, high={self.high}, low={self.low}, close={self.close}, date={self.date})>"


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


def setup_database():
    financial_plugin = FinanceDatabasePlugin(
        "localhost", "financial_dataset", "postgres", "qwerty"
    )
    Dataset.metadata.create_all(financial_plugin.engine)
    Parameters.metadata.create_all(financial_plugin.engine)
    TimeseriesBidAsk.metadata.create_all(financial_plugin.engine)
    TimeseriesPrice.metadata.create_all(financial_plugin.engine)
    TimeseriesOHLC.metadata.create_all(financial_plugin.engine)
    Timeseries.metadata.create_all(financial_plugin.engine)

    from data_plugin.finance_data_plugin import FinanceDataPlugin

    dataPlugin = FinanceDataPlugin(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data"
    )
    with open(dataPlugin.file_pickled, "rb") as fh:
        dataset, metadata = pickle.load(fh)

    assets = []
    for data, info in zip(dataset, metadata):
        assets.append(
            Dataset(
                code=info["code"],
                description=info["code"],
                file=info["file"],
                full_filename=str(info["full_filename"]),
                header=" ".join(info["header"]) if "header" in info else info["format"],
                name="",
                price_multiplicator=info["price_multiplicator"],
                type=info["type"],
            )
        )

    with financial_plugin.Session.begin() as session:
        session.add_all(assets)
        session.commit()


if __name__ == "__main__":
    setup_database()

    financial_plugin = FinanceDatabasePlugin(
        "localhost", "financial_dataset", "postgres", "qwerty"
    )

    with financial_plugin.Session.begin() as session:
        records = session.scalars(select(Dataset)).all()
        print(records)

    from data_plugin.finance_data_plugin import FinanceDataPlugin

    dataPlugin = FinanceDataPlugin(
        "/home/hynek/work/skola/prog/python/transfer_entropy/data"
    )
    with open(dataPlugin.file_pickled, "rb") as fh:
        datasets, metadata_set = pickle.load(fh)

    print(metadata_set)
    for dataset, metadata in zip(datasets, metadata_set):
        record = session.scalars(
            select(Dataset).where(Dataset.file == metadata["file"])
        ).all()
        if record:
            dataset_id = record[0].id
            with financial_plugin.Session.begin() as session:
                new_records = []
                for date, record in dataset.items():
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
