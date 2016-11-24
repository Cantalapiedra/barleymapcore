#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatabasesConfig.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

# Fields in references.conf file
REF_NAME = 0
REF_ID = 1
REF_TYPE = 2

# REF_TYPE values
REF_TYPE_BIG = "big"
REF_TYPE_STD = "std"

class DatabasesConfig(object):
    
    _config_file = ""
    _verbose = False
    _config_dict = {}
    # [ref_id] = {"ref_name", "ref_type"}
    
    def __init__(self, config_file, verbose = True):
        self._config_file = config_file
        self._load_config(config_file)
        self._verbose = verbose
    
    def _load_config(self, config_file):
        self._config_dict = {}
        conf_rows = load_conf(config_file, self._verbose) # data_utils.load_conf
        
        for conf_row in conf_rows:
            config_data = conf_row.strip().split(" ")
            ref_id = config_data[REF_ID]
            ref_name = config_data[REF_NAME]
            ref_type = config_data[REF_TYPE]
            
            self._config_dict[ref_id] = {REF_NAME:ref_name, REF_TYPE:ref_type}
        
    
    def get_databases(self):
        return self._config_dict
    
    def get_database(self, database):
        return self._config_dict[database]
    
    def get_databases_ids(self, databases_names):
        databases_ids = []
        
        # Doing this in a loop to conserve order
        for database_name in databases_names:
            found = False
            for database in self._config_dict:
                if self._config_dict[database][REF_NAME] == database_name:
                    databases_ids.append(database)
                    found = True
                    break
            
            if not found:
                sys.stderr.write("DatabasesConfig: database name "+database_name+" not found in config.\n")
        
        return databases_ids
    
    def get_databases_names(self, databases_ids):
        databases_names = []
        
        for database in databases_ids:
            found = False
            if database in self._config_dict:
                databases_names.append(self._config_dict[database][REF_NAME])
                found = True
            
            if not found:
                sys.stderr.write("DatabasesConfig: database ID "+database+" not found in config.\n")
                databases_names.append(database)
        
        return databases_names
    
    
    
    # Obtain type (big or std size) of fasta DB
    def get_ref_type(self, database_config):
        return database_config[REF_TYPE]

## END