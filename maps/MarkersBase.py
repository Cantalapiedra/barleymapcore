#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MarkersBase.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class MarkerMapping(object):
    _marker_id = ""
    _dataset_name = ""
    _chrom_name = ""
    _chrom_order = -1
    _pos = -1
    
    def __init__(self, marker_id, dataset_name, chrom_name, chrom_order, pos):
        self._marker_id = marker_id
        self._dataset_name = dataset_name
        self._chrom_name = chrom_name
        self._chrom_order = chrom_order
        self._pos = pos
    
    @staticmethod
    def get_empty():
        return MarkerMapping("-", "-", "-", -1, -1)
    
    def __str__(self):
        return " - ".join([self._marker_id, self._dataset_name, str(self._chrom_name)+"/"+str(self._chrom_order), str(self._pos)])
    
    def get_marker_id(self):
        return self._marker_id
    
    def get_dataset_name(self):
        return self._dataset_name
    
    def get_chrom_name(self):
        return self._chrom_name
    
    def get_chrom_order(self):
        return self._chrom_order
    
    def get_pos(self):
        return self._pos
    
class MarkersFields(object):
    
    MARKER_ID_POS = 0
    MARKER_DATASET_POS = 1
    MARKER_CHR_POS = 2
    MARKER_CM_POS = 3
    MARKER_BP_POS = 4
    MARKER_GENES_POS = 5
    MARKER_GENES_CONFIGURED_POS = 6
    
    MARKERS_FIELDS = 7
    
class MarkersData(object):
    
    MARKER_ID_POS = 0
    MARKER_DATASET_POS = 1
    MARKER_GENES_POS = 2
    MARKER_GENES_CONFIGURED_POS = 3
    
## END