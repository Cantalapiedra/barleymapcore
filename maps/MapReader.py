#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapReader.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from MapsBase import MapFile, MapTypes
from barleymapcore.db.MapsConfig import MapsConfig

class MapReader(object):
    
    _maps_path = ""
    _maps_config = None
    _verbose = False
    
    def __init__(self, maps_path, maps_config, verbose = True):
        self._maps_path = maps_path
        self._maps_config = maps_config
        self._verbose = verbose
    
    def __get_position_indexed_data(self, genetic_map, dbs_list, datasets_contig_index, sort_param):
        
        positions_index = {}
        
        map_config = self._maps_config.get_map(genetic_map)#__map_config(genetic_map)
        has_cm_pos = self._maps_config.get_map_has_cm_pos(map_config) #map_config[MapTypes.MAP_HAS_CM_POS]
        has_bp_pos = self._maps_config.get_map_has_bp_pos(map_config) #map_config[MapTypes.MAP_HAS_BP_POS] # "has_bp_pos"
        
        if sort_param == MapTypes.MAP_SORT_PARAM_CM:
            if has_cm_pos:
                position_offset = MapFile.MAP_FILE_CM
            else:
                position_offset = MapFile.MAP_FILE_BP
        elif sort_param == MapTypes.MAP_SORT_PARAM_BP:
            if has_bp_pos:
                position_offset = MapFile.MAP_FILE_BP
            else:
                position_offset = MapFile.MAP_FILE_CM
        
        # For each database configured for this map
        for db in dbs_list:
            if db in self._maps_config.get_map_db_list(map_config):
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
        
        map_config = self._maps_config.get_map(genetic_map)
        
        # Preload contigs to search for
        positions_dict = dict((contig, {"chr":-1, "cm_pos":-1, "bp_pos":-1}) for contig in contig_set)
        test_set = set(contig_set) # A clone of contig_set. Used to shorten the search of contigs
        
        # For this genetic_map, read the info related to each database of contigs
        for db in dbs_list:
            db_records_read = 0
            if db in self._maps_config.get_map_db_list(map_config):
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
                        
                        map_has_cm_pos = self._maps_config.get_map_has_cm_pos(map_config)
                        if map_has_cm_pos:
                            positions_dict[contig_id]["cm_pos"] = float(map_data[MapFile.MAP_FILE_CM])
                        else:
                            positions_dict[contig_id]["cm_pos"] = -1.0
                        
                        map_has_bp_pos = self._maps_config.get_map_has_bp_pos(map_config)
                        if map_has_bp_pos: # "has_bp_pos"
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
        return self._maps_config
    
## END