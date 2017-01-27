#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AlignmentsRetriever.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

from barleymapcore.alignment.AlignmentResult import AlignmentResult

### Class to obtain alignment results from pre-calculated datasets
### "alignment results" are those which have db positions
### but no map positions yet. Therefore, these are
### like the results from bmap_align_to_db or bmap_align_to_map
class AlignmentsParser(object):
    
    #for result in map_results:
    #    #sys.stderr.write("DatasetsFacade.py results "+str(result)+"\n")
    #    db_result = result.get_db_name() #[AlignmentResults.DB_NAME]
    #    if db_result in results:
    #        results[db_result].append(result)
    #    else:
    #        results[db_result] = [result]
    #    
    #retriever = AlignmentsRetriever(self._datasets_config, data_alignments_path, self._verbose)
    
    def parse_alignments_file(self, query_ids_dict, data_path,
                              dataset_synonyms = {}, test_set = None):
        hits = []
        
        for hit in open(data_path, 'r'):
            #sys.stderr.write(str(hit)+"\n")
            hit_data = hit.strip().split("\t")
            alignment_result = AlignmentResult(hit_data)
            hit_query = alignment_result.get_query_id()#[AlignmentResults.QUERY_ID]
            
            if test_set:
                if dataset_synonyms:
                    if hit_query in dataset_synonyms:
                        hit_synonyms = dataset_synonyms[hit_query]
                        synonyms_found = test_set.intersection(hit_synonyms)
                        if len(synonyms_found) > 0:
                            alignment_result.set_query_id("|".join(synonyms_found))
                            hits.append(alignment_result)
                            for synonym in synonyms_found:
                                query_ids_dict[synonym] = 1
                else:
                    if hit_query in test_set:
                        #sys.stderr.write(str(hit_query)+"\n")
                        hits.append(alignment_result)
                        query_ids_dict[hit_query] = 1 # found
                        #if hierarchical: query_ids_dict[hit_query] = 1
            else:
                hits.append(alignment_result)
        
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