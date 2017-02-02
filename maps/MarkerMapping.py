#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MarkerMapping.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class MarkerMapping(object):
    _marker_id = ""
    _dataset_name = ""
    _chrom_name = ""
    _chrom_order = "-1"
    _pos = "-1"
    _empty = False
    
    def __init__(self, marker_id, dataset_name, chrom_name, chrom_order, pos, empty = False):
        self._marker_id = marker_id
        self._dataset_name = dataset_name
        self._chrom_name = chrom_name
        self._chrom_order = chrom_order
        self._pos = pos
        self._empty = empty
    
    # An empty MarkerMapping is used in enriched maps
    # for those mapping positions without markers associated
    @staticmethod
    def get_empty():
        return MarkerMapping("-", "-", "-", "", "", empty = True)
    
    def is_empty(self):
        return self._empty
    
    def set_empty(self, empty):
        self._empty = empty
    
    def __str__(self):
        return " - ".join([self._marker_id, self._dataset_name, str(self._chrom_name)+"/"+str(self._chrom_order), str(self._pos), str(self._empty)])
    
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
    
## END