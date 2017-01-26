#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AlignmentsRetriever.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

### Class to obtain alignment results from pre-calculated datasets
### "alignment results" are those which have db positions
### but no map positions yet. Therefore, these are
### like the results from bmap_align_to_db or bmap_align_to_map
class AlignmentsRetriever(object):
    
    _alignments_results = {}
    _alignments_unmapped = []
    
    _verbose = False
    
    def __init__(self, verbose = False):
        self._verbose = verbose
    
    def get_alignment_results(self):
        return self._alignments_results
    
    def get_alignment_unmapped(self):
        return self._alignments_unmapped
    
    #####################################################
    # Obtain the alignments results from a dataset in a given Map
    #
    def retrieve_alignment_results(self, query_ids_path, dataset_list, map_id,
                         best_score_filter = False):
        results = {}
        num_of_results = 0
        
        self._query_ids_path = query_ids_path
        
        sys.stderr.write("DatasetsFacade: searching "+query_ids_path+"...\n")
        
        # Load list of queries to search for
        initial_num_queries = 0
        query_ids_dict = {}
        for query_ids in open(query_ids_path, 'r'):
            query_ids_dict[query_ids.strip()] = 0
            initial_num_queries += 1
        
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
            
            data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)+".hits"
            
            if self._verbose:
                sys.stderr.write("\t\t MAP: "+map_id+"\n")
                sys.stderr.write("\t\t\t queries to search for: "+str(num_queries_left)+"\n")
                sys.stderr.write("\t\t\t path: "+data_path+"\n")
            
            test_set = set(query for query in temp_query_dict)
            
            #map_results = []
            map_results = self._get_hits(temp_query_dict, data_path, test_set)
            
            num_of_results += len(map_results)
            if self._verbose:
                num_queries = len(set([result.get_query_id() for result in map_results]))
                sys.stderr.write("\t\t\t hits found: "+str(len(map_results))+" for "+str(num_queries)+" queries.\n")
            
            for result in map_results:
                #sys.stderr.write("DatasetsFacade.py results "+str(result)+"\n")
                db_result = result.get_db_name() #[AlignmentResults.DB_NAME]
                if db_result in results:
                    results[db_result].append(result)
                else:
                    results[db_result] = [result]
            
            # end for db_list
            
            query_ids_dict.update(temp_query_dict)
        
        queries_found = initial_num_queries - num_queries_left
        
        if best_score_filter:
            results = self._best_score(results)
        
        if self._verbose: sys.stderr.write("DatasetsFacade: final number of results "+str(num_of_results)+"\n")
        sys.stderr.write("DatasetsFacade: found "+str(queries_found)+" out of "+str(initial_num_queries)+"\n")
        
        self._alignments_results = results
        self._alignments_unmapped = [query for query in query_ids_dict.keys() if query_ids_dict[query] == 0]
        
        return results
    
    def _get_hits(self, query_ids_dict, data_path, test_set):
        hits = []
        
        for hit in open(data_path, 'r'):
            #sys.stderr.write(str(hit)+"\n")
            hit_data = hit.strip().split("\t")
            alignment_result = AlignmentResult(hit_data)
            hit_query = alignment_result.get_query_id()#[AlignmentResults.QUERY_ID]
            
            if hit_query in test_set:
                #sys.stderr.write(str(hit_query)+"\n")
                hits.append(alignment_result)
                query_ids_dict[hit_query] = 1 # found
                #if hierarchical: query_ids_dict[hit_query] = 1
        
        return hits
    
    def _best_score(self, results):
        best_score_filtering = {}
        
        for db in results:
            for result in results[db]:
                query_id = result.get_query_id()#[0]
                align_score = float(result.get_align_score())#[4])
                
                if query_id in best_score_filtering:
                    query_best_score = best_score_filtering[query_id]["best_score"]
                    if align_score < query_best_score:
                        continue
                    elif align_score == query_best_score:
                        best_score_filtering[query_id]["results"].append(result)
                    else: # align_score > query_best_score
                        best_score_filtering[query_id]["results"] = [result]
                        best_score_filtering[query_id]["best_score"] = align_score
                else:
                    best_score_filtering[query_id] = {"results":[result], "best_score":align_score}
            
            results[db] = []
        
        for query_id in best_score_filtering:
            for result in best_score_filtering[query_id]["results"]:
                db = result.get_db_name()#[AlignmentResults.DB_NAME]
                results[db].append(result)
        # else: # NO FILTERING
        
        return results
    
## END