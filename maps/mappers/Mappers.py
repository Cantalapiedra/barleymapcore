#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Mappers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from MapsBase import MapPosition, MappingResults
from barleymapcore.db.MapsConfig import MapsConfig

NUM_FIELDS = 7

class Mapper(object):
    
    _mapReader = None
    _verbose = False
    
    def __init__(self, mapReader, verbose = False):
        self._mapReader = mapReader
        self._verbose = verbose
    
    ## Obtain a finished map of markers from alignments to physical map sequences (genome)
    def create_physical_map(self, alignment_results, unaligned_markers, map_config, sort_param, multiple_param):
        
        # Re-format positions of aligned markers to map positions
        #
        markers_positions = self._reformatPositions(alignment_results, map_config)
        
        finished_map = self._finish_map(markers_positions, unaligned_markers, \
                                        sort_param, multiple_param, map_config)
        
        return finished_map
    
    # Translates positions of alignments to map format
    def _reformatPositions(self, alignment_results, map_config):
        
        # Get full configuration for genetic map
        map_name = map_config.get_name()
        
        markers_positions = {}
        # marker_id --> {"positions":[], "hits_no_position":[], "genetic_map":map_name}
        
        # Extract alignments from databases of this map
        for db in map_config.get_db_list():
            if db in alignment_results:
                
                ## Read the alignment
                for alignment in alignment_results[db]:
                    # Obtain alignment data
                    marker_id = alignment.get_query_id()
                    
                    if marker_id in markers_positions:
                        marker_pos = markers_positions[marker_id]["positions"]
                    else:
                        marker_pos = []
                        markers_positions[marker_id] = {"positions":marker_pos, "hits_no_position":[],
                                                        "genetic_map":map_name}
                    
                    contig_id = alignment.get_subject_id()
                    local_position = alignment.get_local_position()
                    
                    new_pos = {"chr":contig_id, "cm_pos":-1, "bp_pos":local_position}
                    
                    if not self._existPosition(marker_pos, new_pos):
                        marker_pos.append(new_pos)
        
        return markers_positions
    
    ## Obtain a finished map of markers from alignments to anchored sequences
    def create_anchored_map(self, markers_dict, unaligned_markers, map_config, sort_param, multiple_param):
        
        # get the full list of contigs found in the alignments
        contig_list = markers_dict["contig_list"]
        del markers_dict["contig_list"] # This is essential to avoid this key to be used as a marker_id
        
        # Get full configuration for genetic map
        map_name = map_config.get_name()
        
        if self._verbose:
            sys.stderr.write("Mapper: Map: "+str(map_name)+"\n")
            sys.stderr.write("\t total contigs --> "+str(len(contig_list))+"\n")
        
        # Obtain the positions of contigs in the map: contig_id --> chr, cM, bp
        map_contigs_positions = self._mapReader.obtain_map_positions(contig_list)
        
        if self._verbose: sys.stderr.write("\t positions: "+str(len(map_contigs_positions))+"\n")
        
        # Translate positions of aligned markers to map positions,
        # and detect markers which hit to contigs without map position
        markers_positions = self._resolvePositions(map_contigs_positions, markers_dict, map_name)
        
        #sys.stderr.write(str(markers_positions)+"\n")
        
        finished_map = self._finish_map(markers_positions, unaligned_markers, \
                                        sort_param, multiple_param, map_config)
        
        return finished_map
    
    # Obtain a double index of contigs with keys [chromosome][position] --> list of [contig_id, contig_chr, contig_pos]
    def get_positions_index(self, datasets_contig_index, sort_param):
        
        positions_index = self._mapReader.get_positions_index(datasets_contig_index, sort_param)
        
        return positions_index
    
    # Translates positions of aligned markers to map positions,
    # and detects markers which hit to contigs without map position
    def _resolvePositions(self, map_contigs_positions, markers_dict, map_name):
        
        markers_positions = {}
        # marker_id --> {"positions":[], "hits_no_position":[], "genetic_map":map_name}
        
        # For each marker in the alignment results
        for marker_id in markers_dict:
            #sys.stderr.write(marker_id+"\n")
            if marker_id in markers_positions:
                marker_pos = markers_positions[marker_id]["positions"]
            else:
                marker_pos = []
                markers_positions[marker_id] = {"positions":marker_pos, "hits_no_position":[],
                                                "genetic_map":map_name}
            
            #sys.stderr.write("Current "+str(marker_pos)+"\n")
            
            # Associate contigs positions on the map, to the markers
            for contig_tuple in markers_dict[marker_id]:
                
                contig_id = contig_tuple[0]
                local_position = contig_tuple[1]
                #sys.stderr.write("Local position: "+str(local_position)+"\n")
                
                if contig_id in map_contigs_positions:
                    contig_pos = map_contigs_positions[contig_id]
                    
                    #sys.stderr.write("Contig pos: "+str(contig_pos["bp_pos"])+"\n")
                    
                    final_contig_pos = self._clone_contig_pos(contig_pos)
                    
                    #sys.stderr.write("Final position: "+str(final_contig_pos["bp_pos"])+"\n")
                    
                    # Avoid adding twice a position
                    if not self._existPosition(marker_pos, final_contig_pos):
                        marker_pos.append(final_contig_pos)
                else:
                    # Alignments without map position
                    markers_positions[marker_id]["hits_no_position"].append(contig_id)
            
        return markers_positions
    
    def _clone_contig_pos(self, contig_pos):
        
        contig_chr = contig_pos["chr"]
        contig_cm_pos = contig_pos["cm_pos"]
        contig_bp_pos = contig_pos["bp_pos"]
        
        new_contig_pos = {"chr":contig_chr, "cm_pos":contig_cm_pos, "bp_pos":contig_bp_pos}
        
        return new_contig_pos
    
    def _existPosition(self, pos_list, test_pos):
        retValue = False
        
        test_chr = test_pos["chr"]
        test_cm = test_pos["cm_pos"]
        test_bp = test_pos["bp_pos"]
        
        for pos in pos_list:
            pos_chr = pos["chr"]
            pos_cm = pos["cm_pos"]
            pos_bp = pos["bp_pos"]
            
            if pos_chr == test_chr and pos_cm == test_cm and pos_bp == test_bp:
                retValue = True
                break
            
        return retValue
    
    # Creates the final dictionary of markers with map positions
    # including sorting the map, and creating lists of unaligned and unmapped markers
    def _finish_map(self, markers_positions, unaligned_markers, \
                    sort_param, multiple_param, map_config):
        
        ### Create a list of positions from the dictionary of positions of markers (markers_positions)
        # (creates the structure of each position (marker_id, chr, cM, ...))
        chrom_dict = self._mapReader.get_chromosomes_dict()
        positions_list = self._createPositions(markers_positions, multiple_param, chrom_dict)
        
        # Sort the list
        sorted_positions = self._sort_positions_list(positions_list, sort_param)
        
        ### Create a list of markers with alignments but NO map position (unmapped markers),
        # with the final format to display
        unmapped_markers = self._get_unmapped_markers(markers_positions)
        
        ### This is the list of markers without alignment at all
        #
        unaligned_markers = sorted(unaligned_markers, key=lambda sorting: sorting)
        
        # Creates the data structure for the genetic map results
        mapping_results = MappingResults()
        mapping_results.set_mapped(sorted_positions)
        mapping_results.set_fine_mapping(self._single_chromosome(sorted_positions))
        mapping_results.set_unmapped(unmapped_markers)
        mapping_results.set_unaligned(unaligned_markers)
        mapping_results.set_sort_by(sort_param)
        mapping_results.set_map_config(map_config)
        
        # MAP_WITH_GENES and MAP_WITH_MARKERS are assigned in other procedures
        
        return mapping_results
    
    def _createPositions(self, markers_positions, multiple_param, chrom_dict):
        positions_list = []
        
        for marker_id in markers_positions:
            #sys.stderr.write(marker_id+"\n")
            positions = markers_positions[marker_id]["positions"]
            #sys.stderr.write(str(positions)+"\n")
            num_marker_pos = len(positions)
            
            if num_marker_pos == 0: continue # contigs without position
            if (not multiple_param) and (num_marker_pos > 1):
                if self._verbose: sys.stderr.write("Mappers: discarded multiple pos marker: "+str(marker_id)+"\n")
                continue # Multiple positions
            
            num_contig_no_pos = len(markers_positions[marker_id]["hits_no_position"])
            genetic_map_name = markers_positions[marker_id]["genetic_map"]
            
            for pos in positions:
                # marker - chr - cm_pos - bp_pos - multiple - has_contigs_with_no_pos - map_name
                chrom_order = chrom_dict[pos["chr"]] # Numeric value of chromsome (for sorting purposes)
                
                map_position = MapPosition(marker_id, pos["chr"], chrom_order, pos["cm_pos"], pos["bp_pos"],
                                       num_marker_pos > 1, num_contig_no_pos > 0, genetic_map_name)
                positions_list.append(map_position)
        
        return positions_list
    
    def _sort_positions_list(self, positions_list, sort_param):
        sorted_list = []
        
        sorted_list = sorted(positions_list, key=lambda map_position: \
                             (map_position.get_chrom_order(), map_position.get_sort_pos(sort_param),
                              map_position.get_sort_sec_pos(sort_param), map_position.get_marker_id()))
        
        return sorted_list
    
    def _get_unmapped_markers(self, markers_positions):
        positions_list = []
        
        for marker_id in markers_positions:
            hits_no_pos = markers_positions[marker_id]["hits_no_position"]
            num_contig_no_pos = len(hits_no_pos)
            
            if num_contig_no_pos == 0: continue # if has all the contigs with map position, continue
            
            # This is to point out if the marker has other alignments which YES have map position
            num_marker_pos = len(markers_positions[marker_id]["positions"])
            
            for contig in hits_no_pos:
                # marker contig has_other_contigs_with_map_position
                positions_list.append([marker_id, contig, num_marker_pos > 0])
        
        positions_list = sorted(positions_list, key=lambda sorting: (sorting[0], sorting[1], sorting[2]))
        
        return positions_list
    
    # Returns True if only one chromosome found in results
    def _single_chromosome(self, map_positions):
        retValue = True
        
        prev_chr = ""
        for position in map_positions:
            pos_chr = position.get_chrom_name() 
            
            if prev_chr != "":
                if pos_chr != prev_chr:
                    retValue = False
                    break
            
            prev_chr = pos_chr
        
        return retValue
    
## END