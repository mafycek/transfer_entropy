#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from data_plugin.generic_data_plugin import GenericDataPlugin


class GenericDatabasePlugin(GenericDataPlugin):
    def __init__(self, database_url, database, username, password):
        super().__init__()
        self.database_url = database_url
        self.database = database
        self.username = username
        self.password = password
