#!/usr/bin/env python
# -*- coding: utf-8 -*-

# DatasetsFacade.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys, os

from MappingsRetriever import MappingsRetriever
from AlignmentsRetriever import AlignmentsRetriever

from barleymapcore.db.DatasetsConfig import DatasetsConfig
from barleymapcore.alignment.AlignmentResult import *
from barleymapcore.maps.MapInterval import MapInterval
from barleymapcore.maps.MarkersBase import MarkerMapping
from barleymapcore.maps.MappingResults import MappingResult

class DatasetsFacade(object):
    
    #_datasets_conf_file = ""
    _datasets_config = None
    _datasets_path = ""
    _query_ids_path = ""
    _verbose = False
    _genes_hit_dict = {}
    _genes_hit_dict_loaded = False
    
    _mappings_retriever = None
    _alignments_retriever = None
    
    def __init__(self, datasets_config, datasets_path, verbose = True):
        #self._datasets_conf_file = datasets_conf_file
        self._datasets_config = datasets_config
        self._datasets_path = datasets_path
        self._verbose = verbose
        self._mappings_retriever = MappingsRetriever(datasets_path, verbose)
        self._alignments_retriever = AlignmentsRetriever(verbose)
    
    def get_mapping_results(self):
        return self._mappings_retriever.get_mapping_results()
    
    def get_mapping_unmapped(self):
        return self._mappings_retriever.get_mapping_unmapped()
    
    def get_alignment_results(sel):
        return self._alignments_retriever.get_results()
    
    def get_alignment_unmapped(self):
        return self._alignments_retriever.get_unmapped()
    
    #####################################################
    # Obtain the mapping results from a dataset in a given map
    #
    def retrieve_mapping_results(self, query_ids_path, dataset_list, map_config, chrom_dict,
                                 best_score_filter = False, multiple_param = True):
        
        results = self._mappings_retriever.retrieve_mapping_results(query_ids_path, dataset_list, map_config, chrom_dict,
                                                                   best_score_filter, multiple_param)
        
        return results
    
    def parse_mapping_file(self, query_ids_dict, data_path, map_config, chrom_dict, multiple_param, test_set = None):
        
        hits = self._mappings_retriever.parse_mapping_file(query_ids_dict, data_path, map_config, chrom_dict, multiple_param, test_set)
        
        return hits
    
    #####################################################
    # Obtain the alignments results from a dataset in a given Map
    #
    def retrieve_alignment_results(self, query_ids_path, dataset_list, map_id,
                         best_score_filter = False):
        
        results = self._alignments_retriever.retrieve_alignment_results(query_ids_path, dataset_list, map_id,
                                                                        best_score_filter)
        
        return results
    
    ###############################################
    
    ### Obtain markers aligned to a list of contigs
    ### and adds them to the field "markers" of each contig in contig_list
    def retrieve_markers_by_anchor(self, contig_list, map_id):
        
        if self._verbose: sys.stderr.write("DatasetsFacade: loading markers associated to contigs...\n")
        
        ## Prepare a list of contig identifiers
        ##
        contigs_dict = dict([(contig_dict["map_file_contig"],contig_dict) for contig_dict in contig_list])
        
        ## Search datasets for markers
        ## associated to those contigs
        dataset_list = self._datasets_config.get_datasets().keys()
        
        # Look for markers for each dataset
        for dataset in dataset_list:
            
            dataset_config = self._datasets_config.get_dataset_config(dataset)
            dataset_name = dataset_config.get_dataset_name()#datasets_dict[dataset]["dataset_name"]
            dataset_type = dataset_config.get_dataset_type()
            
            ####### IF DATASET IS NOT OF GENETIC MARKERS: EXCLUDE flcDNAs, ESTs, ...
            if dataset_type != DatasetsConfig.DATASET_TYPE_GENETIC_MARKER:
                if self._verbose: sys.stderr.write("\t No markers dataset: "+dataset+"\n")
                continue
            
            if self._verbose: sys.stderr.write("\t dataset: "+dataset+"\n")
            
            ########## Using the dataset generated for the whole map
            ##########
            dataset_map_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)
            if os.path.exists(dataset_map_path) and os.path.isfile(dataset_map_path):
                if self._verbose: sys.stderr.write("DatasetsFacade: loading from map data\n")
                
                self.__retrieve_markers_by_anchor(dataset_map_path, dataset_name, contigs_dict)
            
            ### The next will be deprecated. Generate always a dataset.map file during barleymap setup
            ########## In case that the positions were not precalculated for the map as a whole
            ########## but for individual DBs
            #else:
            #    if self._verbose: sys.stderr.write("DatasetsFacade: loading from DBs data\n")
            #    
            #    self.__retrieve_markers_from_db_files(markers, map_db_list, map_is_hierarchical)
        
        return
    
    # Obtain alignment results from a dataset.map file
    # and add them both to a list (markers) and to a dict of contigs (contigs_dict)
    def __retrieve_markers_by_anchor(self, data_path, dataset_name, contigs_dict):
        
        contigs_set = set(contigs_dict.keys())
        
        # Find all the hits for this map
        for hit in open(data_path, 'r'):
            
            hit_data = hit.strip().split("\t")
            alignment_result = AlignmentResult(hit_data)
            contig = alignment_result.get_subject_id()
            
            if contig in contigs_set:
                query_data = {"alignment_result":alignment_result,
                              "dataset_name":dataset_name,
                              "hit_genes":[],
                              "dataset_genes_configured":False}
                
                contigs_dict[contig]["markers"].append(query_data)
        
        return
    
    ### Obtain markers aligned to a series of alignment intervals
    ###
    def retrieve_markers_by_pos(self, map_intervals, map_config, chrom_dict, map_sort_by):
        if self._verbose: sys.stderr.write("DatasetsFacade: loading markers associated to physical positions...\n")
        
        markers = []
        
        ## Search datasets for markers
        ## associated to those contigs
        dataset_list = self._datasets_config.get_datasets().keys()
        
        # Look for markers for each dataset
        for dataset in dataset_list:
            
            dataset_config = self._datasets_config.get_dataset_config(dataset)
            dataset_name = dataset_config.get_dataset_name()#datasets_dict[dataset]["dataset_name"]
            dataset_type = dataset_config.get_dataset_type()
            
            ####### IF DATASET IS NOT OF GENETIC MARKERS: EXCLUDE flcDNAs, ESTs, ...
            if dataset_type != DatasetsConfig.DATASET_TYPE_GENETIC_MARKER:
                if self._verbose: sys.stderr.write("\t No markers dataset: "+dataset+"\n")
                continue
            
            if self._verbose: sys.stderr.write("\t dataset: "+dataset+"\n")
            
            ########## Retrieve markers within intervals
            ##########
            map_id = map_config.get_id()
            dataset_map_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)
            
            if os.path.exists(dataset_map_path) and os.path.isfile(dataset_map_path):
                if self._verbose: sys.stderr.write("DatasetsFacade: loading from map data: "+dataset_map_path+"\n")
                
                dataset_markers = self.__retrieve_markers_by_pos(dataset_map_path, dataset_name, map_intervals,
                                                                 chrom_dict, map_config, map_sort_by)
                markers.extend(dataset_markers)
        
        return markers
    
    # Obtain alignment results from a dataset.map file
    # and add them both to a list (markers) and to a dict of contigs (contigs_dict)
    def __retrieve_markers_by_pos(self, data_path, dataset_name, map_intervals, chrom_dict, map_config, map_sort_by):
        
        map_name = map_config.get_name()
        map_has_cm_pos = map_config.has_cm_pos()
        map_has_bp_pos = map_config.has_bp_pos()
        
        #contigs_set = set(contigs_dict.keys())
        markers = []
        # Find all the hits for this map
        for hit in open(data_path, 'r'):
            
            if hit.startswith(">") or hit.startswith("#"): continue
            
            hit_data = hit.strip().split("\t")
            
            mapping_result = MappingResult.init_from_data(hit_data, map_name, chrom_dict, map_has_cm_pos, map_has_bp_pos)
            marker_id = mapping_result.get_marker_id()
            chrom_name = mapping_result.get_chrom_name()
            chrom_order = mapping_result.get_chrom_order()
            map_pos = mapping_result.get_sort_pos(map_sort_by)
            
            interval = MapInterval(chrom_name, map_pos, map_pos)
            
            for map_interval in map_intervals:
                
                does_overlap = MapInterval.intervals_overlap(interval, map_interval)
                
                # Check if alignment overlaps with some mapping interval
                if does_overlap:
                    marker_mapping = MarkerMapping(marker_id, dataset_name, chrom_name, chrom_order, map_pos)
                    
                    markers.append(marker_mapping)
                    break # skip intervals, continue with next dataset record
        
        return markers
    
    ## Obtain alignment results from dataset.db files associated to a map configuration
    #def __retrieve_markers_from_db_files(self, markers, map_db_list, map_is_hierarchical):
    #    
    #    for db in map_db_list:
    #        
    #        if hierarchical: markers_list_set = set(markers)
    #        
    #        if self._verbose: sys.stderr.write("\t\t db: "+db+"\n")
    #        
    #        data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(db)+".hits"
    #        
    #        # Find all the hits for this database
    #        for hit in open(data_path, 'r'):
    #            
    #            hit_data = hit.strip().split("\t")
    #            alignment_result = AlignmentResult(hit_data)
    #            
    #            hit_query = alignment_result.get_query_id()
    #            
    #            if hierarchical and hit_query in markers_list_set:
    #                continue
    #            
    #            query_data = {"alignment_result":alignment_result,
    #                          "dataset_name":dataset_name,
    #                          "hit_genes":[],
    #                          "dataset_genes_configured":False}
    #            
    #            markers.append(query_data)
    #    
    #    return
    
    ###############################################
    # Obtains the marker-->genes hits for a dataset
    #def load_dataset_genes_hit(self, dataset):
    #    
    #    dataset_genes_hit_dict = {}
    #    
    #    if self._verbose: sys.stderr.write("DatasetsFacade: loading genes hits for dataset "+str(dataset)+"...\n")
    #    
    #    url = self._datasets_path+"/"+str(dataset)+"/"+str(dataset)+".genes.hits"
    #    dataset_genes_hit_dict["genes_hit"] = {}
    #    
    #    try:
    #        for line in open(url, 'r'):
    #            line_data = line.strip().split("\t")
    #            marker_id = line_data[0]
    #            gene_id = line_data[1]
    #            
    #            if marker_id in dataset_genes_hit_dict["genes_hit"]:
    #                dataset_genes_hit_dict["genes_hit"][marker_id].append(gene_id)
    #            else:
    #                dataset_genes_hit_dict["genes_hit"][marker_id] = [gene_id]
    #        
    #        dataset_genes_hit_dict["configured"] = True
    #        
    #    except IOError:
    #        dataset_genes_hit_dict["configured"] = False
    #    
    #    if self._verbose: sys.stderr.write("DatasetsFacade: loaded genes hits for dataset "+str(dataset)+".\n")
    #    
    #    return dataset_genes_hit_dict
    #
    #def _load_genes_hit_dict(self, dataset_list):
    #    self._genes_hit_dict = {}
    #    
    #    for dataset in dataset_list:
    #        self._genes_hit_dict[dataset] = self.load_dataset_genes_hit(dataset)
    #    
    #    return
    
    
    
    ## Creates an index [contig] --> list of []
    #def load_index_by_contig(self, map_config):
    #    
    #    contig_index = {}
    #    
    #    map_id = map_config[MapsConfig.MAP_ID]
    #    map_db_list = map_config[MapsConfig.DB_LIST]
    #    map_is_hierarchical = map_config[MapsConfig.IS_HIERARCHICAL]
    #    
    #    if self._verbose: sys.stderr.write("DatasetsFacade: loading contig's index with markers data...\n")
    #    
    #    dataset_list = self._datasets_config.get_datasets().keys()
    #    
    #    if not self._genes_hit_dict_loaded:
    #        self._load_genes_hit_dict(dataset_list)
    #        self._genes_hit_dict_loaded = True
    #        
    #    # Look for markers for each dataset
    #    for dataset in dataset_list:
    #        
    #        dataset_config = self._datasets_config.get_dataset(dataset)
    #        dataset_name = self._datasets_config.get_dataset_name(dataset_config)#datasets_dict[dataset]["dataset_name"]
    #        dataset_type = self._datasets_config.get_dataset_type(dataset_config)
    #        
    #        ####### IF DATASET IS NOT OF GENETIC MARKERS: EXCLUDE flcDNAs, ESTs, ...
    #        if dataset_type != DatasetsConfig.DATASET_TYPE_GENETIC_MARKER:
    #            if self._verbose: sys.stderr.write("\t No markers dataset: "+dataset+"\n")
    #            continue
    #        
    #        if self._verbose: sys.stderr.write("\t dataset: "+dataset+"\n")
    #        
    #        dataset_genes_hit_dict = self._genes_hit_dict[dataset]
    #        dataset_genes_configured = dataset_genes_hit_dict["configured"]
    #        dataset_genes_hit_queries = dataset_genes_hit_dict["genes_hit"]
    #        
    #        ########## Using the dataset generated for the whole map
    #        ##########
    #        data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(map_id)+".hits"
    #        if os.path.exists(data_path) and os.path.isfile(data_path):
    #            
    #            if self._verbose: sys.stderr.write("DatasetsFacade: loading from map data\n")
    #            # Find all the hits for this map
    #            for hit in open(data_path, 'r'):
    #                
    #                hit_data = hit.strip().split("\t")
    #                alignment_result = AlignmentResult(hit_data)
    #                
    #                hit_query = alignment_result.get_query_id()
    #                
    #                if hit_query in dataset_genes_hit_queries:
    #                    hit_genes = sorted(dataset_genes_hit_queries[hit_query])
    #                else:
    #                    hit_genes = []
    #                
    #                query_data = {"alignment_result":alignment_result,
    #                              "dataset_name":dataset_name,
    #                              "hit_genes":hit_genes,
    #                              "dataset_genes_configured":dataset_genes_configured}
    #                
    #                hit_contig = alignment_result.get_subject_id()#hit_data[1]
    #                
    #                if hit_contig in contig_index:
    #                    contig_index[hit_contig].append(query_data)
    #                else:
    #                    contig_index[hit_contig] = [query_data]
    #                    
    #        ########## In case that the positions were not precalculated for the map as a whole
    #        ########## but for individual DBs
    #        else:
    #            
    #            # Reset the list of markers for hierarchical mode (applied over databases)
    #            markers_list = [] # For use in hierarchical mode
    #            
    #            if self._verbose: sys.stderr.write("DatasetsFacade: loading from DBs data\n")
    #            for db in map_db_list:
    #                if hierarchical: markers_list_set = set(markers_list)
    #                if self._verbose: sys.stderr.write("\t\t db: "+db+"\n")
    #                
    #                data_path = self._datasets_path+str(dataset)+"/"+str(dataset)+"."+str(db)+".hits"
    #                
    #                tmp_markers_list = []
    #                
    #                # Find all the hits for this database
    #                for hit in open(data_path, 'r'):
    #                    
    #                    hit_data = hit.strip().split("\t")
    #                    alignment_result = AlignmentResult(hit_data)
    #                    
    #                    hit_query = alignment_result.get_query_id()
    #                    
    #                    if hierarchical and hit_query in markers_list_set:
    #                        continue
    #                    
    #                    if hit_query in dataset_genes_hit_queries:
    #                        hit_genes = sorted(dataset_genes_hit_queries[hit_query])
    #                    else:
    #                        hit_genes = []
    #                    
    #                    query_data = {"alignment_result":alignment_result,
    #                                  "dataset_name":dataset_name,
    #                                  "hit_genes":hit_genes,
    #                                  "dataset_genes_configured":dataset_genes_configured}
    #                    
    #                    hit_contig = alignment_result.get_subject_id()#hit_data[1]
    #                    
    #                    if hit_contig in contig_index:
    #                        contig_index[hit_contig].append(query_data)
    #                        tmp_markers_list.append(hit_query)
    #                    else:
    #                        contig_index[hit_contig] = [query_data]
    #                        tmp_markers_list.append(hit_query)
    #                
    #                # Update markers found for hierarchical mode
    #                if hierarchical: markers_list.extend(tmp_markers_list)
    #    
    #    if self._verbose: sys.stderr.write("DatasetsFacade: contig's index with markers data loaded.\n")
    #    
    #    #sys.stderr.write("DatasetsFacade: "+str(contig_index.keys()[0])+" - "+str(contig_index.values()[0])+".\n")
    #    
    #    return contig_index
    
## END