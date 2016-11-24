#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatasetsConfig.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

from barleymapcore.utils.data_utils import load_conf

class DatasetsConfig(object):
    
    DATASET_NAME = 0
    DATASET_ID = 1
    DATASET_TYPE = 2
    FASTA_PATH = 3
    
    # DATASET_TYPE values
    DATASET_TYPE_GENETIC_MARKER = "genetic_marker"
    DATASET_TYPE_RNA = "RNA"
    
    _config_file = ""
    _verbose = False
    _config_dict = {} # dict with data from configuration file (default: conf/datasets.conf)
    
    def __init__(self, config_file, verbose = True):
        self._config_file = config_file
        self._load_config()
        self._verbose = verbose
    
    def _load_config(self):
        self._config_dict = {}
        conf_rows = load_conf(self._config_file, self._verbose) # data_utils.load_conf
        
        #self._config_dict = load_maps(self._config_file, self._verbose) # data_utils.load_maps
        for conf_row in conf_rows:
            
            dataset_name = conf_row[DATASET_NAME]
            dataset_id = conf_row[DATASET_ID]
            dataset_type = conf_row[DATASET_TYPE]
            dataset_fasta_path = conf_row[FASTA_PATH]
            
            dataset_dict = {DATASET_NAME:dataset_name,
                        DATASET_ID:dataset_id,
                        DATASET_TYPE:dataset_type,
                        DATASET_FASTA_PATH:dataset_fasta_path}
            
            self._config_dict[dataset_id] = dataset_dict
    
    def get_datasets(self):
        return self._config_dict
    
    def get_dataset(self, dataset):
        return self._config_dict[dataset]
    
    def get_dataset_name(self, dataset_config):
        return dataset_config[DATASET_NAME]
    
    def get_dataset_type(self, dataset_config):
        return dataset_config[DATASET_TYPE]

## END    