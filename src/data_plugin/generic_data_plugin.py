#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import abc
from abc import ABC, abstractmethod
import os
import datetime
from data_plugin.sample_generator import shuffle_sample


class GenericDataPlugin(ABC):
    def __init__(self):
        pass

    def load_datasets(self):
        pass

    @abc.abstractmethod
    def select_dataset_with_code(self, code, start_date=None, end_date=None):
        pass

    @staticmethod
    def prepare_dataset(
            datasets=None,
            swap_datasets=False,
            shuffle_dataset=False,
            selection1=1,
            selection2=1,
    ):
        filtrated_solution = datasets

        print(
            f"PID:{os.getpid()} {datetime.datetime.now().isoformat()} Shape of solution: {filtrated_solution.shape}",
            flush=True,
        )
        marginal_solution_1 = filtrated_solution[:, 0:selection1]
        marginal_solution_2 = filtrated_solution[
                              :, selection1: selection1 + selection2
                              ]

        if swap_datasets:
            marginal_solution_1, marginal_solution_2 = (
                marginal_solution_2,
                marginal_solution_1,
            )

        if shuffle_dataset:
            marginal_solution_1 = shuffle_sample(marginal_solution_1)

        return marginal_solution_1, marginal_solution_2

    def disconnect(self):
        pass
