#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MappingsParser.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

from barleymapcore.maps.MappingResults import MappingResult
from barleymapcore.maps.MapInterval import MapInterval

### Class to obtain mapping results from pre-calculated datasets
### "mapping results" are those which have already map positions
### like those resulting from running bmap_align to a map
class MappingsParser(object):
    
    def parse_mapping_file(self, data_path, map_config, chrom_dict):
        mapping_results_list = []
        
        map_name = map_config.get_name()
        map_has_cm_pos = map_config.has_cm_pos()
        map_has_bp_pos = map_config.has_bp_pos()
        map_is_physical = map_config.as_physical()
        
        for hit in open(data_path, 'r'):
            if hit.startswith(">") or hit.startswith("#"): continue
            hit_data = hit.strip().split("\t")
            
            mapping_result = MappingResult.init_from_data(hit_data, map_name, chrom_dict, map_is_physical, map_has_cm_pos, map_has_bp_pos)
            mapping_results_list.append(mapping_result)
        
        return mapping_results_list
    
    def parse_mapping_file_by_id(self, query_ids_dict, data_path, map_config, chrom_dict,
                                        multiple_param, dataset_synonyms = {}, test_set = None):
        mapping_results_list = []
        
        map_name = map_config.get_name()
        map_has_cm_pos = map_config.has_cm_pos()
        map_has_bp_pos = map_config.has_bp_pos()
        map_is_physical = map_config.as_physical()
        
        for hit in open(data_path, 'r'):
            #sys.stderr.write(str(hit)+"\n")
            if hit.startswith(">") or hit.startswith("#"): continue
            hit_data = hit.strip().split("\t")
            
            mapping_result = MappingResult.init_from_data(hit_data, map_name, chrom_dict, map_is_physical, map_has_cm_pos, map_has_bp_pos)
            hit_query = mapping_result.get_marker_id()#[AlignmentResults.QUERY_ID]
            has_multiple = mapping_result.has_multiple_pos()
            
            if has_multiple and multiple_param == False: continue
            
            if test_set:
                if hit_query in dataset_synonyms:
                    hit_synonyms = dataset_synonyms[hit_query]
                    synonyms_found = test_set.intersection(hit_synonyms)
                    if len(synonyms_found) > 0:
                        mapping_result.set_marker_id("|".join(synonyms_found))
                        mapping_results_list.append(mapping_result)
                        for synonym in synonyms_found:
                            query_ids_dict[synonym] = 1
                else:
                    if hit_query in test_set:
                        #sys.stderr.write(str(hit_query)+"\n")
                        mapping_results_list.append(mapping_result)
                        query_ids_dict[hit_query] = 1 # found
                        #if hierarchical: query_ids_dict[hit_query] = 1
            else: # retrieve all mapping results
                mapping_results_list.append(mapping_result)
        
        return mapping_results_list
    
    def parse_mapping_file_by_pos(self, map_intervals, data_path, chrom_dict, map_config, map_sort_by):
        mapping_results_list = []
        
        map_name = map_config.get_name()
        map_is_physical = map_config.as_physical()
        map_has_cm_pos = map_config.has_cm_pos()
        map_has_bp_pos = map_config.has_bp_pos()
        
        current_interval_pos = 0
        current_interval = map_intervals[current_interval_pos]
        
        # Find all the hits for this map
        for hit in open(data_path, 'r'):
            if hit.startswith(">") or hit.startswith("#"): continue
            hit_data = hit.strip().split("\t")
            
            mapping_result = MappingResult.init_from_data(hit_data, map_name, chrom_dict, map_is_physical, map_has_cm_pos, map_has_bp_pos)
            
            chrom_name = mapping_result.get_chrom_name()
            
            if chrom_name != current_interval.get_chrom(): continue
            
            map_end_pos = mapping_result.get_sort_end_pos(map_sort_by)
            
            if float(map_end_pos) < float(current_interval.get_ini_pos()): continue
            
            marker_id = mapping_result.get_marker_id()
            chrom_order = mapping_result.get_chrom_order()
            map_pos = mapping_result.get_sort_pos(map_sort_by)#float(mapping_result.get_sort_pos(map_sort_by))
            
            while (float(map_pos) > float(current_interval.get_end_pos())):
                current_interval_pos += 1
                if current_interval_pos >= len(map_intervals): break
                current_interval = map_intervals[current_interval_pos]
            
            if current_interval_pos >= len(map_intervals): break
            
            if float(map_end_pos) < float(current_interval.get_ini_pos()): continue
            
            dataset_interval = MapInterval(chrom_name, map_pos, map_end_pos)
            
            #sys.stderr.write("MappingsParser: by_pos "+str(dataset_interval)+" - "+str(current_interval)+"\n")
            
            does_overlap = MapInterval.intervals_overlap(dataset_interval, current_interval)
            
            # Check if alignment overlaps with some mapping interval
            if does_overlap:
                mapping_results_list.append(mapping_result)
        
        return mapping_results_list

## END