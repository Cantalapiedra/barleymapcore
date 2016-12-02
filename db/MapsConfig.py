#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapsConfig.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from barleymapcore.utils.data_utils import load_conf

MAP_NAME = 0
MAP_ID = 1
HAS_CM_POS = 2
HAS_BP_POS = 3
AS_PHYSICAL = 4
IS_HIERARCHICAL = 5
BEST_SCORE = 6
DB_LIST = 7

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

# BEST_SCORE values
BEST_SCORE_YES = "yes"
#BEST_SCORE_NO = "no"

class MapsConfig(object):
    
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
            
            if conf_row[BEST_SCORE] == BEST_SCORE_YES: best_score = True
            else: best_score = False
            
            map_db_list = conf_row[DB_LIST].split(",")
            
            map_dict = {MAP_NAME:map_name,
                        MAP_ID:map_id,
                        HAS_CM_POS:map_cm,
                        HAS_BP_POS:map_bp,
                        AS_PHYSICAL:map_physical,
                        IS_HIERARCHICAL:map_hierarchical,
                        BEST_SCORE:best_score,
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
    def get_map_name(self, map_config):
        return map_config[MAP_NAME]
    
    def get_map_id(self, map_config):
        return map_config[MAP_ID]
    
    def get_map_has_cm_pos(self, map_config):
        return map_config[HAS_CM_POS]
    
    def get_map_has_bp_pos(self, map_config):
        return map_config[HAS_BP_POS]
    
    def get_map_as_physical(self, map_config):
        return map_config[AS_PHYSICAL]
    
    def get_map_is_hierarchical(self, map_config):
        return map_config[IS_HIERARCHICAL]
    
    def get_map_is_best_score(self, map_config):
        return map_config[BEST_SCORE]
    
    def get_map_db_list(self, map_config):
        return map_config[DB_LIST]
    
    def get_maps_names(self, maps_ids):
        maps_names = []
        
        for map_ids in maps_ids:
            found = False
            if map_ids in self._config_dict:
                maps_names.append(self.get_map_name(self.get_map(map_id)))
                #maps_names.append(self._config_dict[map_ids][MAP_NAME])
                found = True
            
            if not found:
                sys.stderr.write("MapsConfig: map ID "+database+" not found in config.\n")
                maps_names.append(map_ids)
        
        
        return maps_names
    
    def get_maps_ids(self, maps_names):
        maps_ids = []
        
        map_names_set = set([
                            (self.get_map_name(self.get_map(map_id)),map_id)
                            for map_id in self.get_maps()
                            ])
        
        # Doing this in a loop to conserve order
        for map_name in maps_names:
            found = False
            if map_name in map_names_set:
                map_id = map_names_set[map_name]
            #for map_id in self._config_dict:
                #if self.get_map_name(self.get_map(map_id)) == map_name:
                #if self._config_dict[map_id][MAP_NAME] == map_name:
                maps_ids.append(map_id)
                found = True
                break
            
            if not found:
                sys.stderr.write("MapsConfig: map name "+map_name+" not found in config.\n")
        
        return maps_ids
    
## END