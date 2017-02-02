#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapReader.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from MapFiles import MapFile, ChromosomesFile

from barleymapcore.maps.MapsBase import MapTypes

from barleymapcore.db.MapsConfig import MapsConfig

class MapReader(object):
    
    _maps_path = ""
    _map_config = None
    _verbose = False
    _chrom_dict = None
    
    def __init__(self, maps_path, map_config, verbose = True):
        self._maps_path = maps_path
        self._map_config = map_config
        self._verbose = verbose
        self._chrom_dict = self._load_chrom_dict()
    
    def get_map_config(self):
        return self._map_config
    
    def get_chrom_dict(self, ):
        return self._chrom_dict
    
    def _load_chrom_dict(self):#, filter_results = True):
        #
        chrom_dict = {}
        # [chrom_name] = {"chrom_order"}
        
        map_config = self.get_map_config()
        map_id = map_config.get_id()
        
        # File with map-DB positions
        map_path = self._maps_path+map_id+"/"+map_id+ChromosomesFile.FILE_EXT
        if self._verbose: sys.stderr.write("\tMapReader: reading chromosome order from "+map_path+"\n")
        
        # Map data for this database
        for map_line in open(map_path, 'r'):
            map_data = map_line.strip().split("\t")
            
            chrom_name = map_data[ChromosomesFile.CHROM_NAME]
            chrom_order = map_data[ChromosomesFile.CHROM_ORDER]
            
            if chrom_name in chrom_dict:
                raise m2pException("Duplicated chromosome name "+chrom_name+" in "+map_path+".")
            else:
                chrom_dict[chrom_name] = chrom_order
        
        return chrom_dict
    
    def obtain_map_positions(self, contig_set):#, filter_results = True):
        #
        positions_dict = {}
        # [contig_id] = {"chr", "cm_pos", "bp_pos"}
        
        map_config = self.get_map_config()
        map_id = map_config.get_id()
        map_db_list = map_config.get_db_list()
        
        #contig_set = set(contig_list) # A clone of contig_list. Used to shorten the search of contigs
        
        # For this genetic_map, read the info related to each database of contigs
        for db in map_db_list:
            db_records_read = 0
            
            # File with map-DB positions
            map_path = self._maps_path+map_id+"/"+map_id+"."+db
            if self._verbose: sys.stderr.write("\tMapReader: map file --> "+map_path+"\n")
            
            # Map data for this database
            for map_line in open(map_path, 'r'):
                db_records_read += 1
                map_data = map_line.strip().split("\t")
                
                contig_id = map_data[MapFile.MAP_FILE_CONTIG]
                
                # Create positions for this contig
                if contig_id in contig_set:
                    
                    map_pos_chr = map_data[MapFile.MAP_FILE_CHR]
                
                    #if map_pos_chr == -1:
                    #    continue
                    
                    if not contig_id in positions_dict:
                        positions_dict[contig_id] = {}
                    
                    positions_dict[contig_id]["chr"] = map_pos_chr
                    
                    map_has_cm_pos = map_config.has_cm_pos()
                    if map_has_cm_pos:
                        positions_dict[contig_id]["cm_pos"] = map_data[MapFile.MAP_FILE_CM]#float(map_data[MapFile.MAP_FILE_CM])
                    else:
                        positions_dict[contig_id]["cm_pos"] = -1.0
                    
                    map_has_bp_pos = map_config.has_bp_pos()
                    if map_has_bp_pos: # "has_bp_pos"
                        positions_dict[contig_id]["bp_pos"] = map_data[MapFile.MAP_FILE_BP]#long(map_data[MapFile.MAP_FILE_BP])
                    else:
                        positions_dict[contig_id]["bp_pos"] = -1
                        
                    contig_set.remove(contig_id)
                    
                    if len(contig_set) == 0:
                        if self._verbose: sys.stderr.write("\t\t all sequences found -->")
                        break
            
            if self._verbose: sys.stderr.write("\t\t records read: "+str(db_records_read)+"\n")
        
        #if filter_results:
        #    positions_dict = dict((contig, positions_dict[contig]) for contig in positions_dict
        #                                                if positions_dict[contig]["chr"] != -1)
        
        return positions_dict
    
    ### Obtain the contigs in the maps which are within intervals
    ### Creates a list (contig_list) of:
    ### contig_dict = {"map_file_contig":contig, "map_file_chr":chr, "map_file_pos":pos, "markers":[]}
    ###     "markers" will be filled later looking for markers associated to that contig
    def retrieve_contigs(self, map_intervals, map_sort_by):#, filter_results = True):
        #
        contig_list = []
        
        map_config = self.get_map_config()
        map_id = map_config.get_id()
        map_db_list = map_config.get_db_list()
        
        # For this genetic_map, read the info related to each database of contigs
        for db in map_db_list:
            db_records_read = 0
            
            # File with map-DB positions
            map_path = self._maps_path+map_id+"/"+map_id+"."+db
            if self._verbose: sys.stderr.write("\tMapReader: map file --> "+map_path+"\n")
            
            (sort_by, sort_sec_pos) = MapFile.get_sort_pos_contigs(map_sort_by,
                                                                   map_config.has_cm_pos(),
                                                                   map_config.has_bp_pos())
            
            interval_list = [] # To store already checked map positions; assert: map_intervals is sorted by position
            
            # Map data for this database
            for map_line in open(map_path, 'r'):
                
                map_record = map_line.strip().split()
                map_line_contig = map_record[MapFile.MAP_FILE_CONTIG]
                map_line_chr = map_record[MapFile.MAP_FILE_CHR]
                map_line_pos = map_record[sort_by]#float(map_record[sort_by])
                
                #sys.stderr.write("Map position: "+str(map_line_contig)+"\t"+str(map_line_chr)+"\t"+str(map_line_pos)+"\n")
                
                pos_in_intervals = self.__pos_in_intervals(map_intervals, interval_list, map_line_chr, map_line_pos)
                
                if pos_in_intervals:
                    contig_dict = {"map_file_contig":map_line_contig, "map_file_chr":map_line_chr, "map_file_pos":map_line_pos, "markers":[]}
                    #sys.stderr.write("Map record in intervals: "+str(map_line_contig)+"\t"+str(map_line_chr)+"\t"+str(map_line_pos)+"\n")
                    contig_list.append(contig_dict)
                
                if len(interval_list) == len(map_intervals): break
        
        return contig_list
    
    def __pos_in_intervals(self, map_intervals, interval_list, map_line_chr, map_line_pos):
        ret_value = False
        
        interval_set = set(interval_list)
        
        for i, interval in enumerate(map_intervals):
            if i in interval_set: continue
            
            int_chr = interval.get_chrom()
            int_ini_pos = interval.get_ini_pos()
            int_end_pos = interval.get_end_pos()
            
            #sys.stderr.write("Interval chr "+str(int_chr)+" ini "+str(int_ini_pos)+" end "+str(int_end_pos)+"\n")
            
            if (map_line_chr == int_chr):
                if (map_line_pos >= int_ini_pos and map_line_pos <= int_end_pos):
                    ret_value = True
                    break
                elif (map_line_pos < int_end_pos):
                    #interval_list.append(i)
                    break
                else: # map_line_pos > int_end_pos
                    interval_list.append(i)
                    #break
            elif (map_line_chr < int_chr):
                #interval_list.append(i)
                break
            else:
                interval_list.append(i)
                #break
        
        #print interval_list
        
        return ret_value
    
    
    #def __get_position_indexed_data(self, datasets_contig_index, sort_param):
    #    
    #    positions_index = {}
    #    
    #    map_config = self.get_map_config()
    #    map_id = map_config.get_id()
    #    has_cm_pos = map_config.has_cm_pos()
    #    has_bp_pos = map_config.has_bp_pos()
    #    
    #    if sort_param == MapTypes.MAP_SORT_PARAM_CM:
    #        if has_cm_pos:
    #            position_offset = MapFile.MAP_FILE_CM
    #        else:
    #            position_offset = MapFile.MAP_FILE_BP
    #    elif sort_param == MapTypes.MAP_SORT_PARAM_BP:
    #        if has_bp_pos:
    #            position_offset = MapFile.MAP_FILE_BP
    #        else:
    #            position_offset = MapFile.MAP_FILE_CM
    #    
    #    # For each database configured for this map
    #    map_db_list = map_config.get_db_list()
    #    for db in map_db_list:
    #        if self._verbose: sys.stderr.write("\tMapReader: DB: "+db+"\n")
    #        
    #        map_path = self._maps_path+map_id+"/"+map_id+"."+db
    #        if self._verbose: sys.stderr.write("\tMapReader: \t map file --> "+map_path+"\n")
    #        
    #        # For each contig position, find its chromosome index and position data
    #        for map_line in open(map_path, 'r'):
    #            map_data = map_line.strip().split("\t")
    #            contig_id = map_data[MapFile.MAP_FILE_CONTIG]
    #            
    #            # If the contig has not aligned marker associated, it wont yield markers to append
    #            if contig_id not in datasets_contig_index: continue
    #            
    #            contig_chr = int(map_data[MapFile.MAP_FILE_CHR])
    #            
    #            if has_cm_pos:
    #                contig_cm = float(map_data[MapFile.MAP_FILE_CM])
    #            else:
    #                contig_cm = -1.0
    #            
    #            if has_bp_pos:
    #                contig_bp = long(map_data[MapFile.MAP_FILE_BP])
    #            else:
    #                contig_bp = -1
    #            
    #            contig_data = [contig_id, contig_chr, contig_cm, contig_bp]
    #            contig_pos = contig_data[position_offset]
    #            
    #            if contig_chr in positions_index:
    #                chr_dict = positions_index[contig_chr]
    #                if contig_pos in chr_dict:
    #                    chr_dict[contig_pos].append(contig_data)
    #                else:
    #                    chr_dict[contig_pos] = [contig_data]
    #            else:
    #                positions_index[contig_chr] = {}
    #                positions_index[contig_chr][contig_pos] = [contig_data]
    #    
    #    return positions_index
    
    ## Obtain a double index of contigs with keys [chromosome][position] --> list of [contig_id, contig_chr, contig_pos]
    #def get_positions_index(self, datasets_contig_index, sort_param):
    #    positions_index = {}
    #    
    #    # Creates an index (by position) of the contigs of databases list (dbs_list) for this genetic map
    #    positions_index = self.__get_position_indexed_data(datasets_contig_index, sort_param)
    #    
    #    # positions_index[chromsome][poisition] --> list of [contig_id, contig_chr, contig_pos]
    #    
    #    return positions_index
    
## END