#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MapsBase.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class MapTypes(object):
    MAP_SORT_PARAM_CM = "cm"
    MAP_SORT_PARAM_BP = "bp"

class MapHeaders(object):
    
    MARKER_NAME_POS = 0
    MARKER_CHR_POS = 1
    MARKER_CM_POS = 2
    MARKER_BP_POS = 3
    MULTIPLE_POS = 4
    OTHER_ALIGNMENTS_POS = 5
    MAP_NAME = 6
    
    OUTPUT_HEADERS = ["Marker", "chr", "cM", "base_pairs", "multiple_positions", "other_alignments", "Map"]
    GENES_HEADERS = ["Gene", "HC/LC", "Map", "chr", "cM", "bp"]
    ANNOT_HEADERS = ["Description", "InterPro", "Signatures", "PFAM server", "GO terms"]
    MARKERS_HEADERS = ["Marker", "Dataset", "chr", "cM", "bp"]

## END
