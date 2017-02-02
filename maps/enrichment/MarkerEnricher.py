#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MarkerEnricher.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from barleymapcore.maps.MarkerMapping import MarkerMapping
from barleymapcore.maps.MappingResults import MappingResult

ROW_TYPE_POSITION = "pos"
ROW_TYPE_MARKER = "marker"
ROW_TYPE_BOTH = "both"

## Factory
class MarkerEnricherFactory(object):
    @staticmethod
    def get_marker_enricher(mapReader, verbose = False):
        retvalue = None
        
        return PositionsMarkerEnricher(mapReader, verbose)

class MarkerEnricher(object):
    
    _mapReader = None
    _verbose = False
    
    def get_map_reader(self):
        return self._mapReader
    
    def retrieve_markers(self, map_config, map_intervals, datasets_facade, map_sort_by):
        raise m2pException("Method 'retrieve_markers' should be implemented in a class inheriting MarkerEnricher.")
    
    def sort_markers(self, markers):
        markers = sorted(markers, key=lambda marker_mapping: \
                        (int(marker_mapping.get_chrom_order()), float(marker_mapping.get_pos()),
                        marker_mapping.get_dataset_name(), marker_mapping.get_marker_id()))
        
        return markers
    
    def enrich_with_markers(self, mapping_results, markers):
        enriched_map = []
        
        mapped = mapping_results.get_mapped()
        map_sort_by = mapping_results.get_sort_by()
        
        p = 0
        num_pos = len(mapped)
        m = 0
        num_markers = len(markers)
        
        while (p<num_pos and m<num_markers):
            
            # Load position data
            #if p<num_pos:
            map_position = mapped[p]
            map_chrom_name = map_position.get_chrom_name()
            map_chrom_order = map_position.get_chrom_order()
            map_pos = float(map_position.get_sort_pos(map_sort_by))
            #print map_position
            
            #if m<num_markers:
            marker_mapping = markers[m]
            marker_chrom = marker_mapping.get_chrom_name()
            marker_chrom_order = marker_mapping.get_chrom_order()
            marker_pos = float(marker_mapping.get_pos())
            #print marker_mapping
            
            # Create rows of enriched map
            if map_chrom_order < marker_chrom_order:
                # create position
                row_type = ROW_TYPE_POSITION
                p+=1
                
            elif marker_chrom_order < map_chrom_order:
                # create marker
                row_type = ROW_TYPE_MARKER
                m+=1
                
            else: # marker_chrom_order == map_chrom_order
                #print "SAME CHROM"
                #print str(map_pos)+"\t"+str(marker_pos)
                if map_pos < marker_pos:
                    # create position
                    row_type = ROW_TYPE_POSITION
                    p+=1
                    
                elif marker_pos < map_pos:
                    # create marker
                    row_type = ROW_TYPE_MARKER
                    m+=1
                    
                else: # marker_pos == map_pos
                    #print "SAME POS"
                    # create position-marker
                    row_type = ROW_TYPE_BOTH
                    p+=1
                    m+=1
                #print str(row_type)+"\n"
            
            row = self._create_row(map_position, marker_mapping, row_type=row_type)
            
            enriched_map.append(row)
        
        while (p<num_pos):
            # create position
            map_position = mapped[p]
            row = self._create_row(map_position, None, row_type=ROW_TYPE_POSITION)
            enriched_map.append(row)
            p+=1
        
        while (m<num_markers):
            # create marker
            marker_mapping = markers[m]
            row = self._create_row(None, marker_mapping, row_type=ROW_TYPE_MARKER)
            enriched_map.append(row)
            m+=1
        
        mapping_results.set_map_with_markers(enriched_map)
        
        return
    
    def _create_row(self, map_position, marker_mapping, row_type):
        row = None
        
        if row_type == ROW_TYPE_POSITION:
            row = self._create_row_position(map_position)
            
        elif row_type == ROW_TYPE_MARKER:
            row = self._create_row_marker(marker_mapping)
            
        elif row_type == ROW_TYPE_BOTH:
            row = self._create_row_position_marker(map_position, marker_mapping)
            
        else:
            raise m2pException("MapEnricher: unrecognized row type "+str(row_type)+".")
        
        if self._verbose: sys.stderr.write("MapEnricher: new enriched row created: "+str(row)+"\n")
        
        return row
    
    def _create_row_position(self, map_position):
        marker_mapping = MarkerMapping.get_empty()
        
        new_map_position = map_position.clone()
        new_map_position.set_feature(marker_mapping)
        
        return new_map_position
    
    def _create_row_marker(self, marker_mapping):
        map_position = MappingResult.get_empty()
        
        new_map_position = map_position.clone()
        new_map_position.set_feature(marker_mapping)
        
        return new_map_position
    
    def _create_row_position_marker(self, map_position, marker_mapping):
        
        new_map_position = map_position.clone()
        new_map_position.set_feature(marker_mapping)
        
        return new_map_position

class PositionsMarkerEnricher(MarkerEnricher):
    
    def __init__(self, mapReader, verbose = False):
        self._mapReader = mapReader
        self._verbose = verbose
        return
    
    def retrieve_markers(self, map_config, map_intervals, datasets_facade, map_sort_by):
        markers = []
        
        sys.stderr.write("PhysicalEnricher: retrieve markers...\n")
        
        # 1) Obtain the translation to numeric chromosome (for sorting purposes)
        # of chromosome names
        chrom_dict = self._mapReader.get_chrom_dict()
        
        # 2) Obtain the markers in the intervals
        #
        markers = datasets_facade.retrieve_markers_by_pos(map_intervals, map_config, chrom_dict, map_sort_by)
        
        # 3) Sort the list by chrom and position
        markers = self.sort_markers(markers)
        
        return markers

######## ContigsMarkerEnricher will be useful when we want to show
######## markers which hit specific Contigs instead of by map positions.
######## For example, to show the data associated to contigs or to show
######## the markers from alignment results (which report contigs)
######## instead of from mapping results (which report map positions)

class ContigsMarkerEnricher(MarkerEnricher):
    
    def __init__(self, mapReader, verbose):
        self._mapReader = mapReader
        self._verbose = verbose
        return
    
    def retrieve_markers(self, map_config, map_intervals, datasets_facade, map_sort_by):
        
        sys.stderr.write("AnchoredEnricher: retrieve markers...\n")
        
        # 1) Obtain the contigs found in those map intervals
        # and also the translation to numeric chromosome (for sorting purposes)
        # of chromosome names
        contig_list = self._mapReader.retrieve_contigs(map_intervals, map_sort_by)
        
        chrom_dict = self._mapReader.get_chrom_dict()
        
        # 2) Obtain the markers which hit to those contigs and add them to each contig
        #   in contig_list (field "markers" of each contig)
        
        datasets_facade.retrieve_markers_by_anchor(contig_list, map_config)
        
        # 3) Reformat to have the markers but with the map positions of the contigs
        # and sort the list by chrom and position
        markers = self.__get_list_of_markers(contig_list, chrom_dict)
        
        markers = self.sort_markers(markers)
        
        return markers
    
    def __get_list_of_markers(self, contig_list, chrom_dict):
        markers = []
        
        for contig_dict in contig_list:
            
            if len(contig_dict["markers"]) == 0: continue
            
            pos = contig_dict["map_file_pos"]
            chrom = contig_dict["map_file_chr"]
            chrom_order = chrom_dict[chrom]
            for marker in contig_dict["markers"]:
                alignment_result = marker["alignment_result"]
                dataset_name = marker["dataset_name"]
                query_id = alignment_result.get_query_id()
                
                marker_mapping = MarkerMapping(query_id, dataset_name, chrom, chrom_order, pos)
                markers.append(marker_mapping)
        
        return markers

## END