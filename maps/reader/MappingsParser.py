#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MappingsRetriever.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

from barleymapcore.maps.MappingResults import MappingResult

### Class to obtain mapping results from pre-calculated datasets
### "mapping results" are those which have already map positions
### like those resulting from running bmap_align to a map
class MappingsParser(object):
    
    def parse_mapping_file(self, query_ids_dict, data_path, map_config, chrom_dict,
                           multiple_param, dataset_synonyms = {}, test_set = None):
        hits = []
        
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
                        hits.append(mapping_result)
                        for synonym in synonyms_found:
                            query_ids_dict[synonym] = 1
                else:
                    if hit_query in test_set:
                        #sys.stderr.write(str(hit_query)+"\n")
                        hits.append(mapping_result)
                        query_ids_dict[hit_query] = 1 # found
                        #if hierarchical: query_ids_dict[hit_query] = 1
            else: # retrieve all mapping results
                hits.append(mapping_result)
        
        return hits

## END