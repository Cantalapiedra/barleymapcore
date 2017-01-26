#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapMarkers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from Mappers import Mapper
from MarkersBase import MarkersFields, MarkersData
from MapReader import MapReader
from MapEnricher import MapEnricher
from MarkerEnricher import MarkerEnricher

from barleymapcore.genes.GenesFacade import GenesFacade
from barleymapcore.datasets.DatasetsFacade import DatasetsFacade
from barleymapcore.alignment.AlignmentResult import *
from barleymapcore.db.MapsConfig import MapsConfig
from barleymapcore.m2p_exception import m2pException
#from barleymapcore.maps.MarkerEnricher import MarkerEnricher

## Read conf file
ALIGN_ACTION = "align"
FIND_ACTION = "find"

class MapMarkers(object):
    
    #_config_path_dict = {}
    _maps_path = ""
    _map_config = None
    _verbose = False
    _mapReader = None
    
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
    
    def create_map(self, alignment_results, unaligned, sort_param, multiple_param):
        
        map_config = self.get_map_config()
        
        sys.stderr.write("MapMarkers: creating map: "+map_config.get_name()+"\n")
        
        map_as_physical = map_config.as_physical()
        
        ### For physical maps, just re-format alignments to map format
        ###
        if map_as_physical:
            sys.stderr.write("\t"+map_config.get_name()+" is physical.\n")
            
            # Indexes the alignments by marker_id
            #markers_dict = self._get_markers_dict(alignment_results)
            
            # Obtain Mapper without MapReader
            mapper = Mapper(self._mapReader, self._verbose)
            
            # Translate from alignment format to map format
            #self._map = mapper.create_physical_map(markers_dict, unaligned, self._map_config, sort_param, multiple_param)
            self._mapping_results = mapper.create_physical_map(alignment_results, unaligned, map_config, sort_param, multiple_param)
        
        ### For anchored maps, translate positions (from alignment to anchored contigs)
        ### to map positions
        else:
            sys.stderr.write("\t"+map_config.get_name()+" is anchored.\n")
            # Indexes the alignments by marker_id
            markers_dict = self._get_markers_dict(alignment_results)
            
            # Obtain Mapper with MapReader
            mapper = Mapper(self._mapReader, self._verbose)
            
            # Create the map from the alignments
            self._mapping_results = mapper.create_anchored_map(markers_dict, unaligned, map_config, sort_param, multiple_param)
        
        sys.stderr.write("MapMarkers: Map "+map_config.get_name()+" created.\n")
        sys.stderr.write("\n")
        
        return
    
    def _get_markers_dict(self, alignment_results):
        markers_dict = {}
        # [marker_id] = set([(contig_id, local_position),...])
        # ["contig_set"] = {[contig_id,...]}
        
        # This is a temporal list of all the contigs in the alignments
        markers_dict["contig_list"] = []
        
        # Extract alignments from databases of this map
        for db in self.get_map_config().get_db_list():#dbs_list:
            if db in alignment_results:
                #alignments_list.extend(alignment_results[db])
                for alignment in alignment_results[db]:
                    # Obtain alignment data
                    marker_id = alignment.get_query_id()
                    contig_id = alignment.get_subject_id()
                    local_position = alignment.get_local_position()
                    
                    contig_tuple = (contig_id, local_position)
                    
                    # Add hit (contig, position) to list of hits of this marker
                    if marker_id in markers_dict:
                        marker_contig_list = markers_dict[marker_id]
                        marker_contig_set = set(marker_contig_list)
                        if contig_tuple not in marker_contig_set:
                            marker_contig_list.append(contig_tuple)
                    else:
                        markers_dict[marker_id] = [contig_tuple]
                    
                    # Add contig to the general list of contigs found in the alignments
                    if contig_id not in set(markers_dict["contig_list"]):
                        markers_dict["contig_list"].append(contig_id)
        
        return markers_dict
    
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
        
        ########## 1) Load Map Data and translate it to intervals
        ##########
        map_enricher = MapEnricher(mapping_results, self._verbose)
        
        map_intervals = map_enricher.map_to_intervals(extend, extend_window)
        
        ########## 2) Use those intervals to
        ##########      obtain markers within those positions (map.as_physical)
        ##########      obtain contigs within those positions and, afterwards, markers anchored to them (not map.as_physical)
        ########## and add markers to the map
        
        map_enricher.enrich_with_markers(map_intervals, datasets_facade, self._mapReader)
        
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return
    
    def enrich_with_genes(self, show_genes_option, load_annot, genes_extend,
                          genes_window_cm, genes_window_bp, sort_param, constrain_fine_mapping = True):
        
        sys.stderr.write("MapMarkers: adding genes info...\n")
        if self._verbose: sys.stderr.write("\tMode: "+str(show_genes_option)+"\n")
        
        genesManager = GenesFacade(self._config_path_dict, load_annot, show_genes_option, self._verbose)
        
        # For each map
        #for genetic_map in self._genetic_maps_list:
        
        if self._verbose: sys.stderr.write("\tMap: "+self._genetic_map+"\n")
        
        genetic_map_data = self._genetic_map_dict[self._genetic_map]
        
        if constrain_fine_mapping:
            fine_mapping = genetic_map_data.is_fine_mapping() #[MapTypes.FINE_MAPPING]
            
            if self._verbose: sys.stderr.write("\tFine mapping: "+str(fine_mapping)+"\n")
            
            if not fine_mapping: return
        
        # Load genes and annotation
        genesManager.load_data(self._genetic_map, load_annot)
        
        genetic_map_has_cm_pos = self.get_map_config().has_cm_pos() # "has_cm_pos"
        genetic_map_has_bp_pos = self.get_map_config().has_bp_pos() # "has_bp_pos"
        genes_window = self._get_genes_window(genes_extend, sort_param, genes_window_cm, genes_window_bp, \
                                              genetic_map_has_cm_pos, genetic_map_has_bp_pos)
        
        sort_by = genetic_map_data.get_sort_by() #[MapTypes.MAP_SORT_BY]
        
        genetic_map_positions = genetic_map_data.get_mapped() #[MapTypes.MAP_MAPPED]
        
        enriched_positions = genesManager.enrich_by_pos(self._genetic_map, genetic_map_positions, \
                                                        genes_extend, genes_window, \
                                                        sort_by, sort_param, load_annot, \
                                                        genetic_map_has_cm_pos, genetic_map_has_bp_pos)
        
        genetic_map_data.get_map_with_genes() #[MapTypes.MAP_WITH_GENES] = enriched_positions
        
        sys.stderr.write("MapMarkers: Genes added.\n")
        
        return

## END