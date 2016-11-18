#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Mappers.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from barleymapcore.maps.MarkerEnricher import MarkerEnricher
from MapsBase import MapTypes, MapFields

NUM_FIELDS = 7

class Mappers(object):
    
    def __init__(self):
        return
    
    # Factory
    def get_mapper(self, mapReader, enrich = False, merge_maps = False, verbose = False):
        mapper = None
        
        if merge_maps:
            sys.stderr.write("Mappers: Warning: merge mode is deprecated, so regular Mapper will be used instead.\n")
            #mapper = MergedMapper(mapReader, enrich, verbose)
            mapper = Mapper(mapReader, enrich, verbose)
        else:
            mapper = Mapper(mapReader, enrich, verbose)
        
        return mapper

class Mapper(object):
    
    _mapReader = None
    _markerEnricher = None
    _verbose = False
    
    def __init__(self, mapReader, enrich, verbose = False):
        self._mapReader = mapReader
        if enrich:
            self._markerEnricher = MarkerEnricher()
        self._verbose = verbose
    
    def get_genetic_map(self, markers_dict, contig_set, genetic_maps, dbs_list, unmapped_list, sort_param, multiple_param):
        return self._get_genetic_maps(markers_dict, contig_set, genetic_maps, dbs_list, unmapped_list, sort_param, multiple_param)
    
    def _get_genetic_maps(self, markers_dict, contig_set, genetic_map, dbs_list, unmapped_list, sort_param, multiple_param):
        genetic_map_dict = {}
        
        maps_data = self._mapReader.get_maps_data()
        
        #sys.stderr.write(str(maps_data)+"\n")
        
        # Get full configuration for genetic map
        genetic_map_name = maps_data[genetic_map][MapTypes.MAP_NAME]
        genetic_map_as_physical = maps_data[genetic_map][MapTypes.MAP_AS_PHYSICAL]
        genetic_map_has_cm_pos = maps_data[genetic_map][MapTypes.MAP_HAS_CM_POS]
        genetic_map_has_bp_pos = maps_data[genetic_map][MapTypes.MAP_HAS_BP_POS]
        (sort_by, sort_sec_pos) = self._get_sort_pos_map(sort_param, genetic_map_has_cm_pos, genetic_map_has_bp_pos)
        
        if self._verbose: sys.stderr.write("Genetic map: "+str(genetic_map)+"\n")
        if self._verbose: sys.stderr.write("\t parameters: total contigs --> "+str(len(contig_set))+"; genetic_map --> "+genetic_map+"; databases --> "+str(dbs_list)+"\n")
        
        # Obtain the positions of contigs: contig_id --> chr, cM, bp
        positions_dict = self._mapReader.obtain_positions(contig_set, genetic_map, dbs_list)
        
        if self._verbose: sys.stderr.write("\t positions: "+str(len(positions_dict))+"\n")
        
        # Associates positions to markers, detects contigs without map position and resolves redundancy
        markers_positions = self._resolvePositions(positions_dict, markers_dict, genetic_map_name, genetic_map_as_physical)
        
        #sys.stderr.write(str(markers_positions)+"\n")
        
        finished_map = self._finish_map(markers_positions, unmapped_list, sort_by, sort_sec_pos, multiple_param, genetic_map_name, \
                                        genetic_map_as_physical, genetic_map_has_cm_pos, genetic_map_has_bp_pos)
        
        return finished_map
    
    # This function exists in GeneEnricher
    def _get_sort_pos_map(self, sort_param, map_has_cm, map_has_bp):
        sort_by = -1
        sec_pos = -1
        
        # sort type
        if sort_param == "cm":
            if map_has_cm:
                sort_by = MapFields.MARKER_CM_POS
                sec_pos = MapFields.MARKER_BP_POS
            else:
                sort_by = MapFields.MARKER_BP_POS
                sec_pos = MapFields.MARKER_CM_POS
                
        elif sort_param == "bp":
            if map_has_bp:
                sort_by = MapFields.MARKER_BP_POS
                sec_pos = MapFields.MARKER_CM_POS
            else:
                sort_by = MapFields.MARKER_CM_POS
                sec_pos = MapFields.MARKER_BP_POS
        else:
            raise Exception("Mapper: Wrong field for sorting "+str(sort_param))
        
        return (sort_by, sec_pos)
    
    ## Returns a list of [contig_id, contig_chr, contig_pos] between each two consecutive positions in map
    def get_contigs(self, genetic_map, genetic_map_has_cm_pos, genetic_map_has_bp_pos, sorted_positions, \
                    dbs_list, datasets_contig_index, map_sort_by, sort_param, genes_extend, genes_window):
        
        new_positions = []
        
        # Obtain a double index of contigs with keys [chromosome][position] --> list of [contig_id, contig_chr, contig_pos]
        positions_index = self._mapReader.get_positions_index(genetic_map, dbs_list, datasets_contig_index, sort_param)
        
        new_positions = self._markerEnricher.enrich(positions_index, genetic_map_has_cm_pos, genetic_map_has_bp_pos, \
                                                    sorted_positions, map_sort_by, sort_param, genes_extend, genes_window)
        
        return new_positions
    
    def _get_new_contig_pos(self, chrom, cm_pos, bp_pos):
        new_contig_pos = {"chr":chrom, "cm_pos":cm_pos, "bp_pos":bp_pos}
        return new_contig_pos
    
    def _resolvePositions(self, positions_dict, markers_dict, genetic_map_name, genetic_map_as_physical):
        markers_positions = {}
        
        for marker_id in markers_dict:
            #sys.stderr.write(marker_id+"\n")
            if marker_id in markers_positions:
                marker_pos = markers_positions[marker_id]["positions"]
            else:
                marker_pos = []
                markers_positions[marker_id] = {"positions":marker_pos, "hits_no_position":[],
                                                "genetic_map":genetic_map_name}
            
            #sys.stderr.write("Current "+str(marker_pos)+"\n")
            
            # Associate contigs positions to markers
            for contig_data in markers_dict[marker_id]["contigs_set"]:
                contig_id = contig_data[0]
                local_position = contig_data[1]
                #sys.stderr.write("Local position: "+str(local_position)+"\n")
                if contig_id in positions_dict:
                    contig_pos = positions_dict[contig_id]
                    #sys.stderr.write("Contig pos: "+str(contig_pos["bp_pos"])+"\n")
                    
                    final_contig_pos = self._get_new_contig_pos(contig_pos["chr"], contig_pos["cm_pos"],
                                           contig_pos["bp_pos"])
                    
                    if genetic_map_as_physical:
                        #sys.stderr.write("LOCAL POS\n")
                        #sys.stderr.write(str(local_position)+"\n")
                        #sys.stderr.write(str(contig_pos["bp_pos"])+"\n")
                        final_contig_pos["bp_pos"] = contig_pos["bp_pos"] + long(local_position) # physical map
                    
                    #sys.stderr.write("Final position: "+str(final_contig_pos["bp_pos"])+"\n")
                    
                    # Avoid adding twice a position
                    if not self._existPosition(marker_pos, final_contig_pos):
                        marker_pos.append(final_contig_pos)
                else:
                    # Alignments without map position
                    markers_positions[marker_id]["hits_no_position"].append(contig_id)
            
        return markers_positions
    
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
    
    # Returns True if only one chromosome found in results
    def _single_chromosome(self, positions):
        retValue = True
        
        prev_chr = -1
        for position in positions:
            pos_chr = position[MapFields.MARKER_CHR_POS]
            
            if prev_chr != -1:
                if pos_chr != prev_chr:
                    retValue = False
                    break
            
            prev_chr = pos_chr
        
        return retValue
    
    def _finish_map(self, markers_positions, unmapped_list, sort_by, sort_sec_pos, multiple_param, genetic_map_name, \
                    genetic_map_as_physical, genetic_map_has_cm_pos, genetic_map_has_bp_pos):
        genetic_map = {}
        
        # Create a list of positions from the dictionary of positions of markers (markers_positions)
        # (creates the structure of each position (marker_id, chr, cM, ...))
        positions_list = self._splitPositions(markers_positions, multiple_param)
        
        # Sort the list of positions
        sorted_positions = self._sort_positions_list(positions_list, sort_by, sort_sec_pos)
        
        # Create a list of markers with alignments but NO map position, with the final format to display
        marker_no_pos_list = self._get_no_pos_list(markers_positions)
        
        # This is the list of markers without alignment at all
        unmapped_list = sorted(unmapped_list, key=lambda sorting: sorting)
        
        # Creates the data structure for the genetic map results
        genetic_map[MapTypes.MAP_MAPPED] = sorted_positions
        genetic_map[MapTypes.FINE_MAPPING] = self._single_chromosome(sorted_positions)
        genetic_map[MapTypes.MAP_UNMAPPED] = marker_no_pos_list
        genetic_map[MapTypes.MAP_UNALIGNED] = unmapped_list
        genetic_map[MapTypes.MAP_NAME] = genetic_map_name
        genetic_map[MapTypes.MAP_AS_PHYSICAL] = genetic_map_as_physical
        genetic_map[MapTypes.MAP_HAS_CM_POS] = genetic_map_has_cm_pos
        genetic_map[MapTypes.MAP_HAS_BP_POS] = genetic_map_has_bp_pos
        genetic_map[MapTypes.MAP_SORT_BY] = sort_by
        genetic_map[MapTypes.MAP_SORT_SEC_POS] = sort_sec_pos
        
        # MAP_WITH_GENES and MAP_WITH_MARKERS are assigned in other procedures
        
        return genetic_map
    
    def _splitPositions(self, markers_positions, multiple_param):
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
                positions_list.append([marker_id, pos["chr"], pos["cm_pos"], pos["bp_pos"], num_marker_pos > 1, num_contig_no_pos > 0, genetic_map_name])
        
        return positions_list
    
    def _sort_positions_list(self, positions_list, sort_by, sort_sec_pos):
        sorted_list = []
        
        sorted_list = sorted(positions_list, key=lambda sorting: \
                             (sorting[MapFields.MARKER_CHR_POS], sorting[sort_by], sorting[sort_sec_pos], sorting[MapFields.MARKER_NAME_POS]))
        
        return sorted_list
    
    def _get_no_pos_list(self, markers_positions):
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

## Currently discontinued
class MergedMapper(Mapper):
    
    def _get_genetic_maps(self, markers_dict, contig_set, genetic_maps_list, dbs_list, unmapped_list, sort_param, multiple_param):
        genetic_map_dict = {}
        
        markers_to_search_for = markers_dict
        
        final_markers_positions = {}
        markers_positions = {}
        for genetic_map in genetic_maps_list:
            
            if self._verbose: sys.stderr.write("Genetic map: "+str(genetic_map)+"\n")
            if self._verbose: sys.stderr.write("\t parameters: "+str(len(contig_set))+" "+genetic_map+" "+str(dbs_list)+"\n")
            positions_dict = self._mapReader.obtain_positions(contig_set, genetic_map, dbs_list)
            
            if self._verbose: sys.stderr.write("\t positions: "+str(len(positions_dict))+"\n")
            
            markers_positions = self._resolvePositions(positions_dict, markers_to_search_for, genetic_map)
            
            for marker_id in markers_positions:
                marker_dict = markers_positions[marker_id]
                
                # Keep the unmapped ones to next genetic map
                if len(marker_dict["positions"]) == 0 and len(marker_dict["hits_no_position"]) > 0:
                    pass
                # Register mapped ones positions
                else:
                    final_markers_positions[marker_id] = markers_positions[marker_id]
                    del markers_to_search_for[marker_id]
        
        # Recover the unmapped ones
        for marker_id in markers_positions:
            marker_dict = markers_positions[marker_id]
            if len(marker_dict["positions"]) == 0 and len(marker_dict["hits_no_position"]) > 0:
                final_markers_positions[marker_id] = markers_positions[marker_id]
        
        markers_positions = final_markers_positions
        
        finished_map = self._finish_map(markers_positions, unmapped_list, sort_param, multiple_param, "merged", genetic_map_has_bp_pos)
        genetic_map_dict["merged"] = finished_map
        
        return genetic_map_dict
    
## END