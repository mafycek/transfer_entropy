#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import logging
import os
import pickle
from pathlib import Path
import numpy as np

from sample_generator import shuffle_sample
from data_plugin.generic_file_plugin import GenericFilePlugin


class GARCHDataPlugin(GenericFilePlugin):
    def __init__(self, base_directory=Path(__file__).parents[0]):
        super().__init__(base_directory)
        self.directory = self.base_directory / "GARCH_data_reference"
        self.file_pickled = self.base_directory / "GARCH_data_reference" / "dataset.bin"

    def read_header(self, fh):
        return None

    def read_dataset(self, fh, parameter):
        dataset = []
        for line in fh:
            splitted_line = line.split("  ")
            first = float(splitted_line[0])
            second = float(splitted_line[1])
            dataset.append([first, second])

        frame = np.array(dataset)
        return frame

    def load_datasets(self):
        dataset = []
        datafiles = os.listdir(self.directory)
        datafiles = [
            file
            for file in datafiles
            if "bin" not in file and os.path.isfile(self.directory / file)
        ]
        for datafile in datafiles:
            try:
                parameters = {"epsilon": float(datafile)}
                filename = self.directory / datafile
                with open(filename, "rt") as fh:
                    self.read_header(fh)
                    frame = self.read_dataset(fh, parameters)
                    if frame is not None:
                        dataset.append([parameters, frame])
                    else:
                        break
            except EOFError as exc:
                pass

        return dataset

    def prepare_dataset(
            self,
            args,
            index_epsilon,
            datasets=None,
            swap_datasets=False,
            shuffle_dataset=False,
            configuration_of_integration=None,
    ):
        if not args.dataset:
            pass
        else:
            filtrated_solution = datasets[index_epsilon][1].T

            print(
                f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Shape of solution: {filtrated_solution.shape}",
                flush=True,
            )
            marginal_solution_1 = filtrated_solution[0:1, :].T
            marginal_solution_2 = filtrated_solution[1:2, :].T

        if swap_datasets:
            marginal_solution_1, marginal_solution_2 = (
                marginal_solution_2,
                marginal_solution_1,
            )

        if shuffle_dataset:
            marginal_solution_1 = shuffle_sample(marginal_solution_1)

        return marginal_solution_1, marginal_solution_2

    def load_static_dataset(self, args):
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Load dataset",
            flush=True,
        )
        datasets = self.load_datasets()

        if args.dataset_range:
            dataset_start = int(args.dataset_range.split("-")[0])
            dataset_end = int(args.dataset_range.split("-")[1])
            datasets = datasets[dataset_start:dataset_end]

        epsilons = []
        for dataset in datasets:
            epsilons.append(dataset[0]["epsilon"])
        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Epsilons: {epsilons}",
            flush=True,
        )

        return datasets, epsilons


if __name__ == "__main__":
    data_plugin = GARCHDataPlugin()
    dataset = data_plugin.load_datasets()
    print(f"We aggregated {len(dataset)} records")
    data_plugin.pickle_dataset(dataset)
