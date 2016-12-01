#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapMarkers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from Mappers import Mappers
from MapsBase import MapTypes, MapFields
from MarkersBase import MarkersFields, MarkersData
from MapReader import MapReader
from barleymapcore.genes.GenesFacade import GenesFacade
from barleymapcore.datasets.DatasetsFacade import DatasetsFacade
from barleymapcore.alignment.Aligners import AlignmentResults

## Read conf file
ALIGN_ACTION = "align"
FIND_ACTION = "find"

class MapMarkers(object):
    
    #_config_path_dict = {}
    _maps_path = ""
    _maps_config = None
    _genetic_maps_list = ""
    _genetic_map_dict = {}
    _verbose = False
    _maps_data = {}
    _mapReader = None
    
    def __init__(self, maps_path, maps_config, genetic_maps_list, verbose = False):
        self._maps_path = maps_path
        self._maps_config = maps_config
        self._genetic_maps_list = genetic_maps_list
        self._verbose = verbose
    
    def get_genetic_maps(self):
        return self._genetic_map_dict
    
    def create_genetic_maps(self, markers_alignment, unmapped_list, dbs_list,
                            sort_param, multiple_param):
        
        sys.stderr.write("MapMarkers: creating maps...\n")
        
        # Indexes the alignments by marker_id
        markers_dict = self._get_markers_dict(markers_alignment, dbs_list)
        # This is a temporal list of all the contigs in the alignments
        contig_set = markers_dict["contig_set"]
        del markers_dict["contig_set"]
        
        # Load MapReader
        self._mapReader = MapReader(self._maps_path, self._maps_config, self._verbose)
        
        # Obtain Mapper
        mapper = Mappers().get_mapper(self._mapReader, enrich = False, verbose = self._verbose)
        
        # Create the genetic maps from the alignments
        for genetic_map in self._genetic_maps_list:
            self._genetic_map_dict[genetic_map] = mapper.get_genetic_map(markers_dict, contig_set,
                                                                        genetic_map, dbs_list, unmapped_list,
                                                                        sort_param, multiple_param)
        
        sys.stderr.write("MapMarkers: Maps created.\n")
        
        return
    
    def _get_markers_dict(self, markers_alignment, dbs_list):
        markers_dict = {}
        # [marker_id] = {"contigs_set":{[(contig_id, local_position)]}
        # ["contig_set"] = {[contig_id]}
        
        # Extract alignments of the several databases
        alignments_list = []
        for db in dbs_list:
            if db in markers_alignment:
                for alignment in markers_alignment[db]:
                    alignments_list.append(alignment)
        
        # Extract contigs set and marker positions dict
        contig_set = set()
        for marker_alignment in alignments_list:
            marker_id = marker_alignment[AlignmentResults.QUERY_ID]
            contig_id = marker_alignment[AlignmentResults.SUBJECT_ID]
            local_position = marker_alignment[AlignmentResults.START_POSITION]
            
            # Markers dict
            if marker_id in markers_dict:
                marker_dict = markers_dict[marker_id]
                if contig_id not in marker_dict["contigs_set"]:
                    marker_dict["contigs_set"].add((contig_id, local_position))
            else:
                markers_dict[marker_id] = {"contigs_set":set([(contig_id, local_position)])}
            
            # contigs set
            if contig_id not in contig_set:
                contig_set.add(contig_id)
        
        markers_dict["contig_set"] = contig_set
        
        return markers_dict
    
    def enrich_with_markers(self, genes_extend, genes_window_cm, genes_window_bp, sort_param, \
                            dbs_list, datasets_ids, datasets_conf_file, hierarchical, merge_maps, constrain_fine_mapping = True):
        
        sys.stderr.write("MapMarkers: adding other markers...\n")
        
        if len(self._genetic_map_dict)<=0: raise Exception("Genetic maps have not been initiallized.")
        
        # Load DatasetsFacade and datasets data indexed by contig
        datasets_path = self._config_path_dict["app_path"]+self._config_path_dict["datasets_path"]
        facade = DatasetsFacade(datasets_conf_file, datasets_path, verbose = self._verbose)
        datasets_contig_index_loaded = False
        datasets_contig_index = {}
        
        # Load Mapper
        #maps_path = self._config_path_dict["app_path"]+self._config_path_dict["maps_path"]
        maps_config_file = self._config_path_dict["app_path"]+"conf/maps.conf"
        #mapReader = MapReader(maps_path, maps_config_file, self._verbose)
        mapper = Mappers().get_mapper(self._mapReader, enrich = True, merge_maps = merge_maps, verbose = self._verbose)
        
        # For each map
        for genetic_map in self._genetic_maps_list:
            
            if self._verbose: sys.stderr.write("\tMap : "+genetic_map+"\n")
            
            genetic_map_data = self._genetic_map_dict[genetic_map]
            
            if constrain_fine_mapping:
                fine_mapping = genetic_map_data[MapTypes.FINE_MAPPING]
                
                if self._verbose: sys.stderr.write("\tFine mapping: "+str(fine_mapping)+"\n")
                
                if not fine_mapping: continue
            
            genetic_map_has_cm_pos = genetic_map_data[MapTypes.MAP_HAS_CM_POS]
            genetic_map_has_bp_pos = genetic_map_data[MapTypes.MAP_HAS_BP_POS] # "has_bp_pos"
            genetic_map_mapped = genetic_map_data[MapTypes.MAP_MAPPED]
            
            genes_window = self._get_genes_window(genes_extend, sort_param, genes_window_cm, genes_window_bp, \
                                                  genetic_map_has_cm_pos, genetic_map_has_bp_pos)
            
            map_sort_by = genetic_map_data[MapTypes.MAP_SORT_BY]
            
            #sys.stderr.write(str(genetic_map_mapped[0])+"\n")
            
            # Contigs index is only loaded once, and only if there is a genetic_map ready for fine mapping
            # Contains the relation of contigs with markers, so it loads only contigs with markers aligned to them
            if not datasets_contig_index_loaded:
                datasets_contig_index = facade.load_index_by_contig(datasets_ids, dbs_list, hierarchical)
            
            # Obtain the contigs with positions in the interval between each pair of markers
            new_positions = mapper.get_contigs(genetic_map, genetic_map_has_cm_pos, genetic_map_has_bp_pos, genetic_map_mapped, \
                                               dbs_list, datasets_contig_index, map_sort_by, sort_param, genes_extend, genes_window)
            
            # MapFields.MAP_FIELDS points out that the contig_ID, to be replaced by marker information
            # is the first field after the marker map regular information
            translated_positions = self._translate_contigs_to_markers(new_positions, MapFields.MAP_FIELDS, datasets_contig_index)
            
            genetic_map_data[MapTypes.MAP_WITH_MARKERS] = translated_positions
            
        sys.stderr.write("MapMarkers: added other markers.\n")
        
        return
    
    def _get_genes_window(self, genes_extend, sort_param, genes_window_cm, genes_window_bp, genetic_map_has_cm_pos, genetic_map_has_bp_pos):
        
        genes_window = 0
        
        if genes_extend:
            if sort_param == "cm":
                if genetic_map_has_cm_pos:
                    genes_window = genes_window_cm
                else:
                    genes_window = genes_window_bp
                
            elif sort_param == "bp":
                if genetic_map_has_bp_pos:
                    genes_window = genes_window_bp
                else:
                    genes_window = genes_window_cm
            else:
                Exception("find_other_markers: Wrong field for sorting "+str(sort_param))
            
        return genes_window
    
    # Replaces each contig for its corresponding markers
    # contigs_dict contains the contigs and their positions
    # contigs_index relates each contigs with markers
    def _translate_contigs_to_markers(self, positions_list, contig_id_pos, contigs_index):
        
        translated_positions = []
        
        if self._verbose: sys.stderr.write("Replacing "+str(len(positions_list))+" contigs...\n")
        
        num_not_found = 0
        
        for position in positions_list:
            
            marker_id = position[MapFields.MARKER_NAME_POS]
            contig_id = position[contig_id_pos]
            
            if contig_id in contigs_index:
                markers_list = contigs_index[contig_id] # Retrieve data by contig_id
                
                # For each marker associated to this contig a new line is added
                for marker_data in markers_list:
                    
                    marker_new_id = marker_data[MarkersFields.MARKER_ID_POS]
                    
                    #if marker_id == marker_new_id: continue
                    
                    marker_position = position[:contig_id_pos]+\
                    marker_data[MarkersData.MARKER_ID_POS:MarkersData.MARKER_DATASET_POS+1]+\
                    position[contig_id_pos+1:]+\
                    marker_data[MarkersData.MARKER_GENES_POS:MarkersData.MARKER_GENES_CONFIGURED_POS+1]
                    
                    translated_positions.append(marker_position)
            else:
                num_not_found += 1
            
        
        if self._verbose: sys.stderr.write("Contigs not found "+str(num_not_found)+"\n")
        
        return translated_positions
    
    def enrich_with_genes(self, show_genes_option, load_annot, genes_extend,
                          genes_window_cm, genes_window_bp, sort_param, constrain_fine_mapping = True):
        
        sys.stderr.write("MapMarkers: adding genes info...\n")
        if self._verbose: sys.stderr.write("\tMode: "+str(show_genes_option)+"\n")
        
        genesManager = GenesFacade(self._config_path_dict, load_annot, show_genes_option, self._verbose)
        
        # For each map
        for genetic_map in self._genetic_maps_list:
            
            if self._verbose: sys.stderr.write("\tMap: "+genetic_map+"\n")
            
            genetic_map_data = self._genetic_map_dict[genetic_map]
            
            if constrain_fine_mapping:
                fine_mapping = genetic_map_data[MapTypes.FINE_MAPPING]
                
                if self._verbose: sys.stderr.write("\tFine mapping: "+str(fine_mapping)+"\n")
                
                if not fine_mapping: continue
            
            # Load genes and annotation
            genesManager.load_data(genetic_map, load_annot)
            
            genetic_map_has_cm_pos = genetic_map_data[MapTypes.MAP_HAS_CM_POS] # "has_cm_pos"
            genetic_map_has_bp_pos = genetic_map_data[MapTypes.MAP_HAS_BP_POS] # "has_bp_pos"
            genes_window = self._get_genes_window(genes_extend, sort_param, genes_window_cm, genes_window_bp, \
                                                  genetic_map_has_cm_pos, genetic_map_has_bp_pos)
            
            sort_by = genetic_map_data[MapTypes.MAP_SORT_BY]
            
            genetic_map_positions = genetic_map_data[MapTypes.MAP_MAPPED]
            
            enriched_positions = genesManager.enrich_by_pos(genetic_map, genetic_map_positions, \
                                                            genes_extend, genes_window, \
                                                            sort_by, sort_param, load_annot, \
                                                            genetic_map_has_cm_pos, genetic_map_has_bp_pos)
            
            genetic_map_data[MapTypes.MAP_WITH_GENES] = enriched_positions
        
        sys.stderr.write("MapMarkers: Genes added.\n")
        
        return

## END