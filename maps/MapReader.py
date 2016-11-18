#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapReader.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from MapsBase import MapFile, MapTypes

class MapReader(object):
    
    MAP_CONFIG_MAP_NAME = 0
    MAP_CONFIG_MAP_ID = 1
    MAP_CONFIG_HAS_CM_POS = 2
    MAP_CONFIG_HAS_BP_POS = 3
    MAP_CONFIG_AS_PHYSICAL = 4
    MAP_CONFIG_DB_LIST = 5
    
    _maps_path = ""
    _config_file = ""
    _verbose = False
    
    _config_dict = {}
    # [map_id] = {"map_name", "has_cm_pos", "has_bp_pos", "db_list":set([db_id])}
    
    def __init__(self, maps_path, config_file, verbose = True):
        self._maps_path = maps_path
        self._config_file = config_file
        self._load_config()
        self._verbose = verbose
    
    def _load_config(self):
        for config_line in open(self._config_file, 'r'):
            config_data = config_line.strip().split(" ")
            if config_data[self.MAP_CONFIG_HAS_CM_POS] == "cm_true":
                has_cm_pos = True
            else:
                has_cm_pos = False
                
            if config_data[self.MAP_CONFIG_HAS_BP_POS] == "bp_true":
                has_bp_pos = True
            else:
                has_bp_pos = False
            
            if config_data[self.MAP_CONFIG_AS_PHYSICAL] == "physical":
                as_physical = True
            else:
                as_physical = False
                
            self._config_dict[config_data[self.MAP_CONFIG_MAP_ID]] = {"map_name":config_data[self.MAP_CONFIG_MAP_NAME], \
                                                                      MapTypes.MAP_HAS_CM_POS:has_cm_pos, MapTypes.MAP_HAS_BP_POS:has_bp_pos, \
                                                                      MapTypes.MAP_AS_PHYSICAL:as_physical, \
                                                                      "db_list":set(config_data[self.MAP_CONFIG_DB_LIST].split(","))}
        
    
    def __map_config(self, genetic_map):
        
        if genetic_map in self._config_dict:
            map_config = self._config_dict[genetic_map]
        else:
            raise Exception("Genetic map "+genetic_map+" is not in config file.")
        
        return map_config
    
    def __get_position_indexed_data(self, genetic_map, dbs_list, datasets_contig_index, sort_param):
        
        positions_index = {}
        
        map_config = self.__map_config(genetic_map)
        has_cm_pos = map_config[MapTypes.MAP_HAS_CM_POS]
        has_bp_pos = map_config[MapTypes.MAP_HAS_BP_POS] # "has_bp_pos"
        
        if sort_param == "cm":
            if has_cm_pos:
                position_offset = MapFile.MAP_FILE_CM
            else:
                position_offset = MapFile.MAP_FILE_BP
        elif sort_param == "bp":
            if has_bp_pos:
                position_offset = MapFile.MAP_FILE_BP
            else:
                position_offset = MapFile.MAP_FILE_CM
        
        # For each database configured for this map
        for db in dbs_list:
            if db in map_config["db_list"]:
                if self._verbose: sys.stderr.write("\tMapReader: DB: "+db+"\n")
                
                map_path = self._maps_path+genetic_map+"/"+genetic_map+"."+db
                if self._verbose: sys.stderr.write("\tMapReader: \t map file --> "+map_path+"\n")
                
                # For each contig position, find its chromosome index and position data
                for map_line in open(map_path, 'r'):
                    map_data = map_line.strip().split("\t")
                    contig_id = map_data[MapFile.MAP_FILE_MARKER]
                    
                    # If the contig has not aligned marker associated, it wont yield markers to append
                    if contig_id not in datasets_contig_index: continue
                    
                    contig_chr = int(map_data[MapFile.MAP_FILE_CHR])
                    
                    if has_cm_pos:
                        contig_cm = float(map_data[MapFile.MAP_FILE_CM])
                    else:
                        contig_cm = -1.0
                    
                    if has_bp_pos:
                        contig_bp = long(map_data[MapFile.MAP_FILE_BP])
                    else:
                        contig_bp = -1
                    
                    contig_data = [contig_id, contig_chr, contig_cm, contig_bp]
                    contig_pos = contig_data[position_offset]
                    
                    if contig_chr in positions_index:
                        chr_dict = positions_index[contig_chr]
                        if contig_pos in chr_dict:
                            chr_dict[contig_pos].append(contig_data)
                        else:
                            chr_dict[contig_pos] = [contig_data]
                    else:
                        positions_index[contig_chr] = {}
                        positions_index[contig_chr][contig_pos] = [contig_data]
        
        return positions_index
    
    def get_positions_index(self, genetic_map, dbs_list, datasets_contig_index, sort_param):
        positions_index = {}
        
        # Creates an index (by position) of the contigs of databases list (dbs_list) for this genetic map
        positions_index = self.__get_position_indexed_data(genetic_map, dbs_list, datasets_contig_index, sort_param)
        
        # positions_index[chromsome][poisition] --> list of [contig_id, contig_chr, contig_pos]
        
        return positions_index
    
    def obtain_positions(self, contig_set, genetic_map, dbs_list, filter_results = True):
        positions_dict = {}
        # [contig_id] = {"chr", "cm_pos", "bp_pos"}
        
        map_config = self.__map_config(genetic_map)
        
        # Preload contigs to search for
        positions_dict = dict((contig, {"chr":-1, "cm_pos":-1, "bp_pos":-1}) for contig in contig_set)
        test_set = set(contig_set) # A clone of contig_set. Used to shorten the search of contigs
        
        # For this genetic_map, read the info related to each database of contigs
        for db in dbs_list:
            db_records_read = 0
            if db in map_config["db_list"]:
                if self._verbose: sys.stderr.write("\tMapReader: DB: "+db+"\n")
                
                map_path = self._maps_path+genetic_map+"/"+genetic_map+"."+db
                if self._verbose: sys.stderr.write("\tMapReader: map file --> "+map_path+"\n")
                
                # Map data for this database
                for map_line in open(map_path, 'r'):
                    db_records_read += 1
                    map_data = map_line.strip().split("\t")
                    contig_id = map_data[0]
                    
                    if contig_id in test_set:
                        positions_dict[contig_id]["chr"] = int(map_data[MapFile.MAP_FILE_CHR])
                        
                        if map_config[MapTypes.MAP_HAS_CM_POS]:
                            positions_dict[contig_id]["cm_pos"] = float(map_data[MapFile.MAP_FILE_CM])
                        else:
                            positions_dict[contig_id]["cm_pos"] = -1.0
                        
                        if map_config[MapTypes.MAP_HAS_BP_POS]: # "has_bp_pos"
                            positions_dict[contig_id]["bp_pos"] = long(map_data[MapFile.MAP_FILE_BP])
                        else:
                            positions_dict[contig_id]["bp_pos"] = -1
                            
                        test_set.remove(contig_id)
                        
                        if len(test_set) == 0:
                            if self._verbose: sys.stderr.write("\t\t all sequences found -->")
                            break
                
                if self._verbose: sys.stderr.write("\t\t records read: "+str(db_records_read)+"\n")
                
            else:
                if self._verbose: sys.stderr.write("\tMapReader: warning: DB "+db+" is not in config file.\n")
        
        if filter_results:
            positions_dict = dict((contig, positions_dict[contig]) for contig in positions_dict if positions_dict[contig]["chr"] != -1)
        
        return positions_dict
    
    def get_maps_data(self):
        return self._config_dict
    
## END