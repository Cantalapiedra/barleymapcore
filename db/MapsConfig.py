#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapsConfig.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

from barleymapcore.utils.data_utils import load_conf

class MapsConfig(object):
    
    MAP_NAME = 0
    MAP_ID = 1
    HAS_CM_POS = 2
    HAS_BP_POS = 3
    AS_PHYSICAL = 4
    IS_HIERARCHICAL = 5
    DB_LIST = 6
    
    # HAS_CM_POS values
    #HAS_CM_POS_FALSE = "cm_false"
    HAS_CM_POS_TRUE = "cm_true"
    
    # HAS_BP_POS values
    #HAS_BP_POS_FALSE = "bp_false"
    HAS_BP_POS_TRUE = "bp_true"
    
    # AS_PHYSICAL values
    AS_PHYSICAL_TRUE = "physical"
    #AS_PHYSICAL_FALSE = "genetic"
    
    # IS_HIERARCHICAL values
    IS_HIERARCHICAL_TRUE = "hierarchical"
    #IS_HIERARCHICAL_FALSE = "greedy"
    
    _config_file = ""
    _verbose = False
    _config_dict = {} # dict with data from maps configuration file (default: conf/maps.conf)
    
    def __init__(self, config_file, verbose = True):
        self._config_file = config_file
        self._load_config(config_file)
        self._verbose = verbose
    
    def _load_config(self, config_file):
        self._config_dict = {}
        conf_rows = load_conf(config_file, self._verbose) # data_utils.load_conf
        
        #self._config_dict = load_maps(self._config_file, self._verbose) # data_utils.load_maps
        for conf_row in conf_rows:
            
            map_name = conf_row[MAP_NAME]
            map_id = conf_row[MAP_ID]
            
            if conf_row[HAS_CM_POS] == HAS_CM_POS_TRUE: map_cm = True
            else: map_cm = False
            
            if conf_row[HAS_BP_POS] == HAS_BP_POS_TRUE: map_bp = True
            else: map_bp = False
            
            if conf_row[AS_PHYSICAL] == AS_PHYSICAL_TRUE: map_physical = True
            else: map_physical = False
            
            if conf_row[IS_HIERARCHICAL] == IS_HIERARCHICAL_TRUE: map_hierarchical = True
            else: map_hierarchical = False
            
            map_db_list = conf_row[DB_LIST].split(",")
            
            map_dict = {MAP_NAME:map_name,
                        MAP_ID:map_id,
                        HAS_CM_POS:map_cm,
                        HAS_BP_POS:map_bp,
                        AS_PHYSICAL:map_physical,
                        IS_HIERARCHICAL:map_hierarchical,
                        DB_LIST:map_db_list}
            
            self._config_dict[map_id] = map_dict
    
    def get_maps(self):
        return self._config_dict
    
    def get_map(self, genetic_map):        
        if genetic_map in self._config_dict:#self._config_dict:
            map_config = self._config_dict[genetic_map]
        else:
            raise Exception("Genetic map "+genetic_map+" is not in config file.")
        
        return map_config
    
    # These are wrappers to use the config_dict fields just within MapsConfig class
    def get_map_db_list(self, map_config):
        return map_config[DB_LIST]
    
    def get_map_has_cm_pos(self, map_config):
        return map_config[HAS_CM_POS]
    
    def get_map_has_bp_pos(self, map_config):
        return map_config[HAS_BP_POS]
    
    def get_map_is_hierarchical(self, map_config):
        return map_config[IS_HIERARCHICAL]
    
## END