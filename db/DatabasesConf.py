#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatabasesConf.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class DatabasesConf(object):
    
    # Fields in references.conf file
    REF_NAME = 0
    REF_ID = 1
    REF_TYPE = 2
    
    _config_file = ""
    _verbose = False
    _config_dict = {}
    # [ref_id] = {"ref_name", "ref_type"}
    
    def __init__(self, config_file, verbose = True):
        self._config_file = config_file
        self._load_config()
        self._verbose = verbose
    
    def _load_config(self):
        for config_line in open(self._config_file, 'r'):
            config_data = config_line.strip().split(" ")
            ref_id = config_data[REF_ID]
            ref_name = config_data[REF_NAME]
            ref_type = config_data[REF_TYPE]
            # ref_type:
            # ----big (AlignmentFacade._REF_TYPE_BIG)
            # ----std (AlignmentFacade._REF_TYPE_NORMAL)
            
            self._config_dict[ref_id] = {REF_NAME:ref_name, REF_TYPE:ref_type}
        
    
    def get_databases(self):
        return self._config_dict
    
    def get_database(self, database):
        return self._config_dict[database]
    
    # Obtain type (big or std size) of fasta DB
    def get_ref_type(self, database_config):
        return database_config[REF_TYPE]

## END