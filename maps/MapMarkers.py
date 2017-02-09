#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapMarkers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from SearchEngines import SearchEnginesFactory
from reader.MapReader import MapReader
from mappers.Mappers import Mappers
from enrichment.MapEnricher import MapEnricher

from barleymapcore.datasets.DatasetsFacade import DatasetsFacade
from barleymapcore.alignment.AlignmentResult import *
from barleymapcore.db.MapsConfig import MapsConfig
from barleymapcore.m2p_exception import m2pException

## Read conf file
ALIGN_ACTION = "align"
FIND_ACTION = "find"

class MapMarkers(object):
    
    #_config_path_dict = {}
    _maps_path = ""
    _map_config = None
    _facade = None
    _verbose = False
    _mapReader = None
    
    _mapping_results = None
    _maps_data = {}
    
    def __init__(self, maps_path, map_config, facade = None, verbose = False):
        self._maps_path = maps_path
        self._map_config = map_config
        self._facade = facade
        self._verbose = verbose
        # Load MapReader
        self._mapReader = MapReader(self._maps_path, map_config, self._verbose)
    
    def get_mapping_results(self):
        return self._mapping_results
    
    def get_map_config(self):
        return self._map_config
    
    # previously: setup_map
    def retrieve_mappings(self, query_ids_path, datasets_ids, sort_param, multiple_param):
        
        search_engine = SearchEnginesFactory.get_search_engine_datasets(self._maps_path, self._verbose)
        
        mapping_results = search_engine.create_map(query_ids_path, datasets_ids, self._map_config, self._facade,
                                                   sort_param, multiple_param)
        
        sys.stderr.write("MapMarkers: Map "+self._map_config.get_name()+" created.\n")
        sys.stderr.write("\n")
        
        self._mapping_results = mapping_results
        
        return mapping_results
    
    def perform_mappings(self, query_fasta_path, databases_ids, databases_config, aligner_list,
                                threshold_id, threshold_cov, n_threads,
                                best_score_param, sort_param, multiple_param, tmp_files_dir):
        
        search_type = self._map_config.get_search_type()
        
        search_engine = SearchEnginesFactory.get_search_engine(search_type, self._maps_path,
                                                               best_score_param, databases_config, aligner_list,
                                                               threshold_id, threshold_cov, n_threads, self._verbose)
        
        mapping_results = search_engine.create_map(query_fasta_path, databases_ids, self._map_config, self._facade,
                                                   sort_param, multiple_param, tmp_files_dir)
        
        sys.stderr.write("MapMarkers: Map "+self._map_config.get_name()+" created.\n")
        sys.stderr.write("\n")
        
        self._mapping_results = mapping_results
        
        ############# Perform alignments
        ## Perform alignments
        ## TODO: avoid aligning to the same DB as one of a previous map
        ## this "TODO" would need to handle correctly best_score and hierarchical
        #facade.perform_alignment(query_fasta_path, databases_ids, hierarchical, query_mode,
        #                                   threshold_id, threshold_cov, n_threads, \
        #                                   best_score)
        #
        #results = facade.get_alignment_results()
        #unaligned = facade.get_alignment_unmapped()  
        #
        ############# MAPS
        #mapMarkers = MapMarkers(maps_path, map_config, verbose_param)
        #
        #mapMarkers.create_map(results, unaligned, sort_by, multiple_param)
        #
        ## enrichment (OLD)
        #
        #mapping_results = mapMarkers.get_mapping_results()
        
        return mapping_results
    
    #def create_map(self, alignment_results, unaligned, sort_param, multiple_param):
    #    
    #    map_config = self.get_map_config()
    #    
    #    sys.stderr.write("MapMarkers: creating map: "+map_config.get_name()+"\n")
    #    
    #    map_as_physical = map_config.as_physical()
    #    
    #    mapper = Mappers.get_alignments_mapper(map_as_physical, self._mapReader, self._verbose)
    #    
    #    self._mapping_results = mapper.create_map(alignment_results, unaligned, map_config, sort_param, multiple_param)
    #    
    #    sys.stderr.write("MapMarkers: Map "+map_config.get_name()+" created.\n")
    #    sys.stderr.write("\n")
    #    
    #    return
    
    #
    def enrichment(self, annotator, show_markers, show_genes, datasets_facade, extend_window, collapsed_view, constrain_fine_mapping = False):
        if show_markers and not show_genes:
            
            # Enrich with markers
            self.enrich_with_markers(datasets_facade, extend_window,
                                            collapsed_view, constrain_fine_mapping)
            
        ########### GENES
        if show_genes:
            
            # Enrich with genes
            self.enrich_with_genes(datasets_facade, extend_window,
                                         annotator, collapsed_view, constrain_fine_mapping)
        
        return
    
    def enrich_with_markers(self, datasets_facade, extend_window, \
                            collapsed_view = False, constrain_fine_mapping = True):
        
        sys.stderr.write("MapMarkers: adding other markers...\n")
        
        mapping_results = self.get_mapping_results()
        
        if not mapping_results: raise m2pException("MapMarkers: Map has not been initiallized.")
        
        sys.stderr.write("\tMap : "+self.get_map_config().get_name()+"\n")
        
        # Markers enrichment is limited to results from a single chromosome
        # for example, in the web app
        if constrain_fine_mapping:
            fine_mapping = mapping_results.is_fine_mapping() # boolean: results from a single chromosome
            if self._verbose: sys.stderr.write("\tFine mapping: "+str(fine_mapping)+"\n")
            if not fine_mapping: return
        
        map_enricher = MapEnricher(mapping_results, self._verbose)
        
        ########## 1) Load Map Data and translate it to intervals
        ##########
        
        map_intervals = map_enricher.map_to_intervals(extend_window)
        
        ########## 2) Use those intervals to
        ##########      obtain markers within those positions (map.as_physical)
        ##########      obtain contigs within those positions and, afterwards, markers anchored to them (not map.as_physical)
        ########## and add markers to the map
        
        map_enricher.enrich_with_markers(map_intervals, datasets_facade, self._mapReader, collapsed_view)
        
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return
    
    def enrich_with_genes(self, datasets_facade, extend_window, \
                            annotator, collapsed_view = False, constrain_fine_mapping = True):
        
        sys.stderr.write("MapMarkers: adding other genes...\n")
        
        mapping_results = self.get_mapping_results()
        
        if not mapping_results: raise m2pException("MapMarkers: Map has not been initiallized.")
        
        sys.stderr.write("\tMap : "+self.get_map_config().get_name()+"\n")
        
        # Markers enrichment is limited to results from a single chromosome
        # for example, in the web app
        if constrain_fine_mapping:
            fine_mapping = mapping_results.is_fine_mapping() # boolean: results from a single chromosome
            if self._verbose: sys.stderr.write("\tFine mapping: "+str(fine_mapping)+"\n")
            if not fine_mapping: return
        
        map_enricher = MapEnricher(mapping_results, self._verbose)
        
        ########## 1) Load Map Data and translate it to intervals
        ##########
        
        map_intervals = map_enricher.map_to_intervals(extend_window)
        
        ########## 2) Use those intervals to
        ##########      obtain markers within those positions (map.as_physical)
        ##########      obtain contigs within those positions and, afterwards, markers anchored to them (not map.as_physical)
        ########## and add markers to the map
        
        map_enricher.enrich_with_genes(map_intervals, datasets_facade, self._mapReader, annotator, collapsed_view)
        
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return

## END