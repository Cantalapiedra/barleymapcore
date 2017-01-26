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
class MappingsRetriever(object):
    
    _datasets_path = ""
    _mapping_results = None
    _mapping_unmapped = None
    
    _verbose = False
    
    def __init__(self, datasets_path, verbose = False):
        self._datasets_path = datasets_path
        self._verbose = verbose
    
    def get_mapping_results(self):
        return self._mapping_results
    
    def get_mapping_unmapped(self):
        return self._mapping_unmapped
    
    #####################################################
    # Obtain the mapping results from a dataset in a given map
    #
    def retrieve_mapping_results(self, query_ids_path, dataset_list, map_config, chrom_dict, best_score_filter = False, multiple_param = True):
        results = []
        
        self._query_ids_path = query_ids_path
        
        map_id = map_config.get_id()
        
        sys.stderr.write("MappingsRetriever: searching "+query_ids_path+"...\n")
        
        # Load list of queries to search for
        initial_num_queries = 0
        query_ids_dict = {}
        for query_ids in open(query_ids_path, 'r'):
            query_ids_dict[query_ids.strip()] = 0
            initial_num_queries += 1
        
        num_results = 0
        map_results = []
        num_queries_left = initial_num_queries
        for dataset in dataset_list:
            sys.stderr.write("\t dataset: "+dataset+"\n")
            
            # This is hierarchical or not,
            # it is done because each marker has to have a uniq identifier
            # so doing this we expect to do the search faster
            temp_query_dict = dict([(query, 0) for query in query_ids_dict if query_ids_dict[query] == 0])
            num_queries_left = len(temp_query_dict)
            
            if num_queries_left == 0:
                if self._verbose: sys.stderr.write("\t\t All queries found.\n")
                break
            
            data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)
            
            if self._verbose:
                sys.stderr.write("\t\t MAP: "+map_id+"\n")
                sys.stderr.write("\t\t\t queries to search for: "+str(num_queries_left)+"\n")
                sys.stderr.write("\t\t\t path: "+data_path+"\n")
                
            test_set = set(query for query in temp_query_dict)
            
            #map_results = []
            
            map_results = self.parse_mapping_file(temp_query_dict, data_path, map_config, chrom_dict, multiple_param, test_set)
            
            num_results += len(map_results)
            if self._verbose:
                num_queries = len(set([result.get_marker_id() for result in map_results]))
                sys.stderr.write("\t\t\t hits found: "+str(len(map_results))+" for "+str(num_queries)+" queries.\n")
            
            query_ids_dict.update(temp_query_dict)
            
            results.extend(map_results)
        
        queries_found = initial_num_queries - num_queries_left
        
        # best score: there is no need to look for best score
        # since precalculated map should have taken that into account
        
        if self._verbose: sys.stderr.write("MappingsRetriever: final number of results "+str(num_results)+"\n")
        sys.stderr.write("MappingsRetriever: found "+str(queries_found)+" out of "+str(initial_num_queries)+"\n")
        
        self._mapping_results = results
        self._mapping_unmapped = [query for query in query_ids_dict.keys() if query_ids_dict[query] == 0]
        
        return results
    
    def parse_mapping_file(self, query_ids_dict, data_path, map_config, chrom_dict, multiple_param, test_set = None):
        hits = []
        
        map_name = map_config.get_name()
        map_has_cm_pos = map_config.has_cm_pos()
        map_has_bp_pos = map_config.has_bp_pos()
        
        for hit in open(data_path, 'r'):
            #sys.stderr.write(str(hit)+"\n")
            if hit.startswith(">") or hit.startswith("#"): continue
            hit_data = hit.strip().split("\t")
            
            mapping_result = MappingResult.init_from_data(hit_data, map_name, chrom_dict, map_has_cm_pos, map_has_bp_pos)
            hit_query = mapping_result.get_marker_id()#[AlignmentResults.QUERY_ID]
            has_multiple = mapping_result.has_multiple_pos()
            
            if has_multiple and multiple_param == False: continue
            
            if test_set:
                if hit_query in test_set:
                    #sys.stderr.write(str(hit_query)+"\n")
                    hits.append(mapping_result)
                    query_ids_dict[hit_query] = 1 # found
                    #if hierarchical: query_ids_dict[hit_query] = 1
            else: # retrieve all mapping results
                hits.append(mapping_result)
                query_ids_dict[hit_query]
        
        return hits

## END