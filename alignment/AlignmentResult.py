#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AlignmentResult.py is part of Barleymap.
# Copyright (C)  2017  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class AlignmentResult(object):
    #QUERY_ID = 0
    #SUBJECT_ID = 1
    #ALIGN_IDENTITY = 2
    #QUERY_COVERAGE = 3
    #ALIGNMENT_SCORE = 4
    #STRAND = 5
    #QSTART_POS = 6
    #QEND_POS = 7
    #START_POSITION = 8
    #END_POSITION = 9
    #DB_NAME = 10
    #ALGORITHM = 11
    
    _query_id = ""
    _subject_id = ""
    _align_ident = 0.0
    _query_cov = 0.0
    _align_score = 0
    _strand = "+"
    _local_position = 0
    _end_position = 0
    _qstart_pos = 0
    _qend_pos = 0
    _db_name = ""
    _algorithm = ""
    
    def __init__(self, alignment_data):
        
        self._query_id = alignment_data[0]
        self._subject_id = alignment_data[1]
        self._align_ident = alignment_data[2]
        self._query_cov = alignment_data[3]
        self._align_score = alignment_data[4]
        self._strand = alignment_data[5]
        self._local_position = alignment_data[6]
        self._end_position = alignment_data[7]
        
        self._qstart_pos = alignment_data[8]
        self._qend_pos = alignment_data[9]
        
        self._db_name = alignment_data[10]
        self._algorithm = alignment_data[11]
    
    def get_query_id(self):
        return self._query_id
    
    def get_subject_id(self):
        return self._subject_id
    
    def get_align_ident(self):
        return self._align_ident
    
    def get_query_cov(self):
        return self._query_cov
    
    def get_align_score(self):
        return self._align_score
    
    def get_strand(self):
        return self._strand
    
    def get_local_position(self):
        return self._local_position
    
    def get_end_position(self):
        return self._end_position
    
    def get_qstart_pos(self):
        return self._qstart_pos
    
    def get_qend_pos(self):
        return self._qend_pos
    
    def get_db_name(self):
        return self._db_name
    
    def get_algorithm(self):
        return self._algorithm
    
## END