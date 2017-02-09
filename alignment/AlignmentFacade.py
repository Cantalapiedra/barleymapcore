#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AlignmentFacade.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# Copyright (C)  2016-2017 Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import os, sys

from AlignmentEngines import AlignmentEnginesFactory
from barleymapcore.db.DatabasesConfig import REF_TYPE_STD

class AlignmentFacade():
    
    _paths_config = None
    
    _alignment_results = None
    
    _verbose = False
    
    def __init__(self, paths_config, verbose = False):
        self._paths_config = paths_config
        self._verbose = verbose
    
    # Performs the alignment of fasta sequences different DBs
    def perform_alignment(self, query_fasta_path, dbs_list, databases_config, search_type, aligner_list, \
                          threshold_id = 98, threshold_cov = 95, n_threads = 1, ref_type_param = REF_TYPE_STD):
        
        ## Create the SearchEngine (greedy, hierarchical, exhaustive searches on top of splitblast, gmap,...)
        alignment_engine = AlignmentEnginesFactory.get_alignment_engine(search_type, aligner_list, self._paths_config, 
                                                               ref_type_param, n_threads, self._verbose)
        
        ## Perform the search and alignments
        alignment_results = alignment_engine.perform_alignment(query_fasta_path, dbs_list, databases_config, threshold_id, threshold_cov)
        
        self._alignment_results = alignment_results
        
        return alignment_results
    
    def get_alignment_results(self):
        return self._alignment_results
    
## END