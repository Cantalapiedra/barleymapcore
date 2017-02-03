#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatasetsRetriever.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

from MappingsRetriever import MappingsParser

from barleymapcore.db.DatasetsConfig import DatasetsConfig
from barleymapcore.m2p_exception import m2pException

class DatasetsRetriever(object):
    
    _datasets_config = None
    _datasets_path = None
    _verbose = False
    
    _results = None
    _unmapped = None
    
    def __init__(self, datasets_config, datasets_path, verbose = False):
        self._datasets_config = datasets_config
        self._datasets_path = datasets_path
        self._verbose = verbose
    
    def load_synonyms(self, synonyms):
        dataset_synonyms = {}
        
        if synonyms != "" and synonyms != DatasetsConfig.SYNONYMS_NO:
            for syn_line in open(synonyms, 'r'):
                syn_data = syn_line.strip().split()
                syn_key = syn_data[0]
                if syn_key in dataset_synonyms:
                    raise m2pException("Repeated synonyms entry for marker "+syn_key+".")
                else:
                    dataset_synonyms[syn_key] = syn_data
        
        return dataset_synonyms
    
    def get_results(self):
        retvalue = None
        
        if self._results != None:
            retvalue = self._results
        else:
            raise m2pException("DatasetsRetriever: error obtaining unloaded results. Call a retrieve method first.")
        
        return retvalue
    
    def get_unmapped(self):
        retvalue = None
        
        if self._unmapped != None:
            retvalue = self._unmapped
        else:
            raise m2pException("DatasetsRetriever: error obtaining unloaded list of unmapped markers. Call a retrieve method first.")
        
        return retvalue
    
    def retrieve_datasets(self, query_ids_path, dataset_list, map_config, chrom_dict,
                                 best_score_filter = False, multiple_param = True):
        self._results = []
        
        self._query_ids_path = query_ids_path
        
        map_id = map_config.get_id()
        map_db_list = map_config.get_db_list()
        
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
            
            dataset_config = self._datasets_config.get_dataset_config(dataset)
            dataset_db_list = dataset_config.get_db_list()
            
            common_dbs = set(dataset_db_list).intersection(set(map_db_list))
            
            if (DatasetsConfig.DATABASES_ANY not in dataset_db_list) and (len(common_dbs)<1):
                continue
            
            synonyms_path = dataset_config.get_synonyms()
            dataset_synonyms = self.load_synonyms(synonyms_path)
            
            # This is hierarchical or not,
            # it is done because each marker has to have a uniq identifier
            # so doing this we expect to do the search faster
            temp_query_dict = dict([(query, 0) for query in query_ids_dict if query_ids_dict[query] == 0])
            num_queries_left = len(temp_query_dict)
            
            if num_queries_left == 0:
                if self._verbose: sys.stderr.write("\t\t All queries found.\n")
                break
            
            test_set = set(query for query in temp_query_dict)
            
            data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)
            if self._verbose:
                sys.stderr.write("\t\t MAP: "+map_id+"\n")
                sys.stderr.write("\t\t\t queries to search for: "+str(num_queries_left)+"\n")
                sys.stderr.write("\t\t\t path: "+data_path+"[.hits]\n")
            
            ############ Retrieve dataset markers
            ############ either from mappings or from alignments
            data_mappings_path = data_path # mappings file
            data_alignments_path = data_path+".hits" # alignments file
            
            if os.path.exists(data_mappings_path): # mapping results are available
                
                mappings_parser = MappingsParser()
                map_results = mappings_parser.parse_mapping_file(temp_query_dict, data_mappings_path, map_config, chrom_dict,
                                                      multiple_param, dataset_synonyms, test_set)
                #retriever = MappingsRetriever(self._datasets_config, data_mappings_path, self._verbose)
                
            ### THIS COULD BE IMPLEMENTED AGAIN
            ###
            #elif os.path.exists(data_alignments_path): # alignment results are available
            #    
            #    alignments_parser = AlignmentsParser()
            #    map_results = alignments_parser.parse_alignments_file(temp_query_dict, data_alignments_path,
            #                                             dataset_synonyms, test_set)
            #    
            else:
                # TODO refactor to handled exception
                sys.stderr.write("WARNING: DatasetsRetriever: there is no available data for dataset "+dataset+"\n")
            
            #retriever.retrieve_datasets(query_ids_path, [dataset], map_config, chrom_dict, best_score_filter, multiple_param)
            
            num_results += len(map_results)
            if self._verbose:
                num_queries = len(set([result.get_marker_id() for result in map_results]))
                sys.stderr.write("\t\t\t hits found: "+str(len(map_results))+" for "+str(num_queries)+" queries.\n")
            
            query_ids_dict.update(temp_query_dict)
            
            self._results.extend(map_results)
        
        temp_query_dict = dict([(query, 0) for query in query_ids_dict if query_ids_dict[query] == 0])
        num_queries_left = len(temp_query_dict)
        
        queries_found = initial_num_queries - num_queries_left
        
        if self._verbose: sys.stderr.write("DatasetsRetriever: final number of results "+str(num_results)+"\n")
        sys.stderr.write("DatasetsRetriever: found "+str(queries_found)+" out of "+str(initial_num_queries)+"\n")
        
        self._unmapped = [query for query in query_ids_dict.keys() if query_ids_dict[query] == 0]
        
        return

## END