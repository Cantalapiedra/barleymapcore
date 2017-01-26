#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MappingResults.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

## This class represents the map position of a marker which has been aligned first to a DB
##
class MappingResult(object):
    _marker_id = ""
    _chrom_name = ""
    _chrom_order = -1
    _cm_pos = -1.0
    _bp_pos = -1
    _multiple_pos = False
    _other_alignments = False
    _map_name = ""
    _feature = None
    
    MAP_FIELDS = 7
    
    def __init__(self, marker_id, chrom_name, chrom_order, cm_pos, bp_pos, has_multiple_pos, has_other_alignments, _map_name):
        self._marker_id = marker_id
        self._chrom_name = chrom_name
        self._chrom_order = chrom_order
        self._cm_pos = cm_pos
        self._bp_pos = bp_pos
        self._multiple_pos = has_multiple_pos
        self._other_alignments = has_other_alignments
        self._map_name = _map_name
    
    @staticmethod
    def get_empty():
        return MapPosition("-", "-", -1, -1.0, -1, False, False, "")
    
    def get_marker_id(self):
        return self._marker_id
    
    def get_chrom_name(self):
        return self._chrom_name
    
    def get_chrom_order(self):
        return self._chrom_order
    
    def get_cm_pos(self):
        return self._cm_pos
    
    def get_bp_pos(self):
        return self._bp_pos
    
    def has_multiple_pos(self):
        return self._multiple_pos
    
    def has_other_alignments(self):
        return self._other_alignments
    
    def get_map_name(self):
        return self._map_name
    
    def get_sort_pos(self, sort_by):
        ret_value = -1
        
        if sort_by == MapTypes.MAP_SORT_PARAM_CM:
            ret_value = float(self._cm_pos)
        elif sort_by == MapTypes.MAP_SORT_PARAM_BP:
            ret_value = long(self._bp_pos)
        else:
            raise m2pException("Unrecognized sort field "+str(sort_by)+".")
        
        return ret_value
    
    def get_sort_sec_pos(self, sort_by):
        ret_value = -1
        
        if sort_by == MapTypes.MAP_SORT_PARAM_CM:
            ret_value = self._bp_pos
        elif sort_by == MapTypes.MAP_SORT_PARAM_BP:
            ret_value = self._cm_pos
        else:
            raise m2pException("Unrecognized sort field "+str(sort_by)+".")
        
        return ret_value
    
    def __str__(self):
        return " - ".join([self._marker_id, str(self._chrom_name)+"/"+str(self._chrom_order), str(self._cm_pos), str(self._bp_pos), str(self._feature)])
    
    def set_feature(self, feature):
        self._feature = feature
    
    def get_feature(self):
        return self._feature
    
##############################
## A class with the results of barleymap
## including the mapped, unmapped and unaligned markers
class MappingResults(object):
    _mapped = None
    _unmapped = None
    _unaligned = None
    _fine_mapping = False
    _sort_by = ""
    _map_config = None
    
    _map_with_genes = None
    _map_with_markers = None
    
    def __init__(self):
        return
    
    def get_mapped(self):
        return self._mapped
    
    def set_mapped(self, mapped):
        self._mapped = mapped
    
    def get_unmapped(self):
        return self._unmapped
    
    def set_unmapped(self, unmapped):
        self._unmapped = unmapped
    
    def get_unaligned(self):
        return self._unaligned
    
    def set_unaligned(self, unaligned):
        self._unaligned = unaligned
    
    def is_fine_mapping(self):
        return self._fine_mapping
    
    def set_fine_mapping(self, fine_mapping):
        self._fine_mapping = fine_mapping
    
    def get_sort_by(self):
        return self._sort_by
    
    def set_sort_by(self, sort_by):
        self._sort_by = sort_by
    
    def get_map_config(self):
        return self._map_config
    
    def set_map_config(self, map_config):
        self._map_config = map_config
    
    def set_map_with_genes(self, map_with_genes):
        self._map_with_genes = map_with_genes
    
    def get_map_with_genes(self):
        return self._map_with_genes
    
    def set_map_with_markers(self, map_with_markers):
        self._map_with_markers = map_with_markers
    
    def get_map_with_markers(self):
        return self._map_with_markers
    
## END