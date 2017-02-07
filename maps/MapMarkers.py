#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapMarkers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

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
    _verbose = False
    _mapReader = None
    _chrom_dict = None
    
    _mapping_results = None
    _maps_data = {}
    
    def __init__(self, maps_path, map_config, verbose = False):
        self._maps_path = maps_path
        self._map_config = map_config
        self._verbose = verbose
        # Load MapReader
        self._mapReader = MapReader(self._maps_path, map_config, self._verbose)
    
    def get_mapping_results(self):
        return self._mapping_results
    
    def get_map_config(self):
        return self._map_config
    
    def get_chrom_dict(self):
        return self._mapReader.get_chrom_dict()
    
    def setup_map(self, query_ids_path, datasets_ids, datasets_facade, best_score_filter, sort_param, multiple_param):
        
        map_config = self.get_map_config()
        
        sys.stderr.write("MapMarkers: setting up map: "+map_config.get_name()+"\n")
        
        map_as_physical = map_config.as_physical()
        chrom_dict = self.get_chrom_dict()
        
        ############ Retrieve pre-computed alignments
        datasets_facade.retrieve_datasets(query_ids_path, datasets_ids, map_config, chrom_dict, best_score_filter, multiple_param)
        
        mapping_results = datasets_facade.get_results()
        mapping_unmapped = datasets_facade.get_unmapped()
        
        # Obtain Mapper without MapReader
        mapper = Mappers.get_mappings_mapper(self._mapReader, self._verbose)
        
        # Translate from alignment format to map format
        #self._map = mapper.create_physical_map(markers_dict, unaligned, self._map_config, sort_param, multiple_param)
        self._mapping_results = mapper.create_map(mapping_results, mapping_unmapped, map_config, sort_param)
        
        return
    
    def create_map(self, alignment_results, unaligned, sort_param, multiple_param):
        
        map_config = self.get_map_config()
        
        sys.stderr.write("MapMarkers: creating map: "+map_config.get_name()+"\n")
        
        map_as_physical = map_config.as_physical()
        
        mapper = Mappers.get_alignments_mapper(map_as_physical, self._mapReader, self._verbose)
        
        self._mapping_results = mapper.create_map(alignment_results, unaligned, map_config, sort_param, multiple_param)
        
        sys.stderr.write("MapMarkers: Map "+map_config.get_name()+" created.\n")
        sys.stderr.write("\n")
        
        return
    
    def enrich_with_markers(self, datasets_facade, extend, extend_window, \
                            constrain_fine_mapping = True):
        
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
        
        map_intervals = map_enricher.map_to_intervals(extend, extend_window)
        
        ########## 2) Use those intervals to
        ##########      obtain markers within those positions (map.as_physical)
        ##########      obtain contigs within those positions and, afterwards, markers anchored to them (not map.as_physical)
        ########## and add markers to the map
        
        map_enricher.enrich_with_markers(map_intervals, datasets_facade, self._mapReader)
        
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return
    
    def enrich_with_genes(self, datasets_facade, extend, extend_window, \
                            annotator, constrain_fine_mapping = True):
        
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
        
        map_intervals = map_enricher.map_to_intervals(extend, extend_window)
        
        ########## 2) Use those intervals to
        ##########      obtain markers within those positions (map.as_physical)
        ##########      obtain contigs within those positions and, afterwards, markers anchored to them (not map.as_physical)
        ########## and add markers to the map
        
        map_enricher.enrich_with_genes(map_intervals, datasets_facade, self._mapReader, annotator)
        
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return

## END