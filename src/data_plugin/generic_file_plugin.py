#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import abstractmethod
from pathlib import Path
import pickle

from data_plugin.generic_data_plugin import GenericDataPlugin


class GenericFilePlugin(GenericDataPlugin):
    def __init__(self, base_directory=Path(__file__).parents[0]):
        super().__init__()
        self.base_directory = Path(base_directory)
        self.directory = self.base_directory / "finance_data_reference"
        self.file_pickled = self.directory / "dataset.bin"

    @abstractmethod
    def load_datasets(self):
        pass

    def pickle_dataset(self, dataset):
        with open(self.file_pickled, "wb") as file_handler:
            pickle.dump(dataset, file_handler)
