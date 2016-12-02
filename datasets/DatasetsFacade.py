#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatasetsFacade.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os
import barleymapcore.db.DatasetsConfig
from barleymapcore.alignment.Aligners import AlignmentResults

class DatasetsFacade(object):
    
    #_datasets_conf_file = ""
    _datasets_config = None
    _datasets_path = ""
    _results = {}
    _query_ids_path = ""
    _unmapped = []
    _verbose = False
    _genes_hit_dict = {}
    _genes_hit_dict_loaded = False
    
    def __init__(self, datasets_config, datasets_path, verbose = True):
        #self._datasets_conf_file = datasets_conf_file
        self._datasets_config = datasets_config
        self._datasets_path = datasets_path
        self._verbose = verbose
    
    # Obtains the marker-->genes hits for a dataset
    def load_dataset_genes_hit(self, dataset):
        
        dataset_genes_hit_dict = {}
        
        if self._verbose: sys.stderr.write("DatasetsFacade: loading genes hits for dataset "+str(dataset)+"...\n")
        
        url = self._datasets_path+"/"+str(dataset)+"/"+str(dataset)+".genes.hits"
        dataset_genes_hit_dict["genes_hit"] = {}
        
        try:
            for line in open(url, 'r'):
                line_data = line.strip().split("\t")
                marker_id = line_data[0]
                gene_id = line_data[1]
                
                if marker_id in dataset_genes_hit_dict["genes_hit"]:
                    dataset_genes_hit_dict["genes_hit"][marker_id].append(gene_id)
                else:
                    dataset_genes_hit_dict["genes_hit"][marker_id] = [gene_id]
            
            dataset_genes_hit_dict["configured"] = True
            
        except IOError:
            dataset_genes_hit_dict["configured"] = False
        
        if self._verbose: sys.stderr.write("DatasetsFacade: loaded genes hits for dataset "+str(dataset)+".\n")
        
        return dataset_genes_hit_dict
    
    def _load_genes_hit_dict(self, dataset_list):
        self._genes_hit_dict = {}
        
        for dataset in dataset_list:
            self._genes_hit_dict[dataset] = self.load_dataset_genes_hit(dataset)
        
        return
    
    # Creates an index [contig] --> list of []
    def load_index_by_contig(self, dataset_list, db_list, hierarchical = True):
        
        contig_index = {}
        
        if self._verbose: sys.stderr.write("DatasetsFacade: loading contig's index with markers data...\n")
        
        #datasets_config = DatasetsConfig(self._datasets_conf_file, self._verbose)
        datasets_dict = self._datasets_config.get_datasets()
        
        if not self._genes_hit_dict_loaded:
            self._load_genes_hit_dict(dataset_list)
            self._genes_hit_dict_loaded = True
            
        # Look for markers for each dataset
        for dataset in dataset_list:
            
            dataset_config = self._datasets_config.get_dataset(dataset)
            dataset_name = self._datasets_config.get_dataset_name(dataset_config)#datasets_dict[dataset]["dataset_name"]
            dataset_type = self._datasets_config.get_dataset_type(dataset_config)
            
            ####### IF DATASET IS OF GENETIC MARKERS: EXCLUDE flcDNAs, HarvEST, ...
            if dataset_type != DatasetsConfig.DATASET_TYPE_GENETIC_MARKER:
                if self._verbose: sys.stderr.write("\t No markers dataset: "+dataset+"\n")
                continue
            
            if self._verbose: sys.stderr.write("\t dataset: "+dataset+"\n")
            
            dataset_genes_hit_dict = self._genes_hit_dict[dataset]
            dataset_genes_configured = dataset_genes_hit_dict["configured"]
            dataset_genes_hit_queries = dataset_genes_hit_dict["genes_hit"]
            
            # Reset the list of markers for hierarchical mode (applied over databases)
            markers_list = [] # For use in hierarchical mode
            
            for db in db_list:
                if hierarchical: markers_list_set = set(markers_list)
                if self._verbose: sys.stderr.write("\t\t db: "+db+"\n")
                
                data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(db)+".hits"
                
                tmp_markers_list = []
                
                # Find all the hits for this database
                for hit in open(data_path, 'r'):
                    hit_data = hit.strip().split("\t")
                    hit_query = hit_data[0]
                    
                    if hierarchical and hit_query in markers_list_set:
                        continue
                    
                    hit_contig = hit_data[1]
                    
                    if hit_query in dataset_genes_hit_queries:
                        hit_genes = sorted(dataset_genes_hit_queries[hit_query])
                    else:
                        hit_genes = []
                    
                    query_data = [hit_query, dataset_name, hit_genes, dataset_genes_configured]
                    
                    if hit_contig in contig_index:
                        contig_index[hit_contig].append(query_data)
                        tmp_markers_list.append(hit_query)
                    else:
                        contig_index[hit_contig] = [query_data]
                        tmp_markers_list.append(hit_query)
                
                # Update markers found for hierarchical mode
                if hierarchical: markers_list.extend(tmp_markers_list)
        
        if self._verbose: sys.stderr.write("DatasetsFacade: contig's index with markers data loaded.\n")
        
        #sys.stderr.write("DatasetsFacade: "+str(contig_index.keys()[0])+" - "+str(contig_index.values()[0])+".\n")
        
        return contig_index
    
    def retrieve_markers(self, contigs_dict):
        pass
    
    # Obtain the alignments from a dataset in a given Map
    def retrieve_ids_map(self, query_ids_path, dataset_list, map_id,
                         best_score_filter = False):
        results = {}
        num_of_results = 0
        
        self._query_ids_path = query_ids_path
        
        sys.stderr.write("DatasetsFacade: searching...\n")
        
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
            
            #for db in db_list: --> map_id
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
            if self._verbose: sys.stderr.write("\t\t\t hits found: "+str(len(map_results))+"\n")
            
            #sys.stderr.write("DatasetsFacade.py map_results "+str(map_results)+"\n")
            for result in map_results:
                #sys.stderr.write("DatasetsFacade.py results "+str(result)+"\n")
                db_result = result[AlignmentResults.DB_NAME]
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
        
        self._results = results
        self._unmapped = [query for query in query_ids_dict.keys() if query_ids_dict[query] == 0]
        
        return results
    
    def _best_score(self, results):
        best_score_filtering = {}
        
        for db in results:
            for result in results[db]:
                query_id = result[0]
                align_score = float(result[4])
                
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
                db = result[AlignmentResults.DB_NAME]
                results[db].append(result)
        # else: # NO FILTERING
        
        return results
    
    def _get_hits(self, query_ids_dict, data_path, test_set):
        hits = []
        
        for hit in open(data_path, 'r'):
            #sys.stderr.write(str(hit)+"\n")
            hit_data = hit.strip().split("\t")
            hit_query = hit_data[AlignmentResults.QUERY_ID]
            
            if hit_query in test_set:
                #sys.stderr.write(str(hit_query)+"\n")
                hits.append(hit_data)
                query_ids_dict[hit_query] = 1 # found
                #if hierarchical: query_ids_dict[hit_query] = 1
        
        return hits
    
    def get_alignment_unmapped(self):
        return self._unmapped
    
## END