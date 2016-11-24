#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapsBase.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

# MapFile: fields of original maps of contigs (IBSC2012 Morex contigs, IBSC2016 chromosomes, ...)
class MapFile(object):
    
    # These are the fields in the map files
    MAP_FILE_MARKER = 0
    MAP_FILE_CHR = 1
    MAP_FILE_CM = 2
    MAP_FILE_BP = 3

# MapFields in table of results from barleymap
class MapFields(object):
    
    MAP_FIELDS = 7
    # These are the fields of master rows (markers on map)
    MARKER_NAME_POS = 0
    MARKER_CHR_POS = 1
    MARKER_CM_POS = 2
    MARKER_BP_POS = 3
    MULTIPLE_POS = 4
    OTHER_ALIGNMENTS_POS = 5
    MAP_NAME = 6
    
class MapTypes(object):
    
    # The different result types
    MAP_MAPPED = "sorted_map"
    MAP_UNMAPPED = "marker_no_pos_list"
    MAP_UNALIGNED = "marker_unmapped_list"
    MAP_WITH_GENES = "map_with_genes"
    MAP_WITH_MARKERS = "map_with_markers"
    
    RESULTS_LIST = [MAP_MAPPED, MAP_WITH_GENES, MAP_WITH_MARKERS, MAP_UNMAPPED, MAP_UNALIGNED]
    RESULTS_DICT = {MAP_MAPPED:"Map", MAP_WITH_GENES:"Map with genes", MAP_WITH_MARKERS:"Map with markers", \
                    MAP_UNMAPPED:"Alignments without map position", MAP_UNALIGNED:"Unaligned markers"}
    
    FINE_MAPPING = "fine_mapping" # boolean: false if more than one chromosome found in MAP_MAPPED
    MAP_NAME = "map_name"
    MAP_HAS_CM_POS = "has_cm_pos"
    MAP_HAS_BP_POS = "has_bp_pos"
    MAP_AS_PHYSICAL = "as_physical" # whether add local_position from alignment to position of contig in map
    MAP_SORT_BY = "sort_by"
    MAP_SORT_SEC_POS = "sort_sec_pos"
    MAP_SORT_PARAM_CM = "cm"
    MAP_SORT_PARAM_BP = "bp"
    

class MapHeaders(object):
    
    OUTPUT_HEADERS = ["Marker", "chr", "cM", "base pairs", "multiple positions", "other alignments", "Map"]
    GENES_HEADERS = ["Gene", "HC/LC", "Map", "chr", "cM", "bp"]
    ANNOT_HEADERS = ["Description", "InterPro", "Signatures", "PFAM server", "GO terms"]
    MARKERS_HEADERS = ["Marker", "Dataset", "chr", "cM", "bp", "genes"]

## END
