#!/usr/bin/env python
# -*- coding: utf-8 -*-

# FeatureMapping.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).
    
# This is a class with features which can be attached to
# a MappingResult, including genes or markers in the same position, etc.
class FeatureMapping(object):
    _feature_id = ""
    _dataset_name = ""
    _chrom_name = ""
    _chrom_order = "-1"
    _pos = "-1"
    _feature_type = "" # genetic_marker, gene, ... see DatasetsConfig
    _empty = False
    
    def __init__(self, feature_id, dataset_name, chrom_name, chrom_order, pos, feature_type, empty = False):
        self._feature_id = feature_id
        self._dataset_name = dataset_name
        self._chrom_name = chrom_name
        self._chrom_order = chrom_order
        self._pos = pos
        self._feature_type = feature_type
        self._empty = empty
    
    # An empty FeatureMapping is used in enriched maps
    # for those mapping positions without features associated
    @staticmethod
    def get_empty():
        return FeatureMapping("-", "-", "-", "", "", "", empty = True)
    
    def is_empty(self):
        return self._empty
    
    def set_empty(self, empty):
        self._empty = empty
    
    def __str__(self):
        return " - ".join([self._feature_id, self._dataset_name,
                           str(self._chrom_name)+"/"+str(self._chrom_order), str(self._pos),
                           str(self._feature_type), str(self._empty)])
    
    def get_feature_id(self):
        return self._feature_id
    
    def get_dataset_name(self):
        return self._dataset_name
    
    def get_chrom_name(self):
        return self._chrom_name
    
    def get_chrom_order(self):
        return self._chrom_order
    
    def get_pos(self):
        return self._pos
    
    def get_feature_type(self):
        return self._feature_type
    
class GeneMapping(FeatureMapping):
    _annots = {}
    
    def get_annots(self):
        return self._annots
    
    def set_annot(self, annots):
        self._annots = annots
    
    def add_annot(self, annot):
        self._annots[annot.get_annot_id()] = annot
    
    def get_annot(self, annot_id):
        return self._annots[annot_id]

class GeneAnnotation(object):
    _id = ""
    _value = ""
    _type = ""
    
    def __init__(self, annot_id, annot_value, annot_type):
        self._annot_id = annot_id
        self._annot_value = annot_value
        self._annot_type = annot_type
    
    def get_annot_id(self):
        return self._id
    
## END