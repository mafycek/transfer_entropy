#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from data_plugin.generic_data_plugin import GenericDataPlugin


class GenericDatabasePlugin(GenericDataPlugin):
    def __init__(self, database_url, database_table_name, username, password):
        super().__init__()
        self.database_url = database_url
        self.database_table_name = database_table_name
        self.username = username
        self.password = password

    def reconnect(self):
        pass

    @abc.abstractmethod
    def update_state_of_calculation(self):
        pass
