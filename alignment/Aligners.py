#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aligners.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import os, sys

import m2p_split_blast
import m2p_gmap
import m2p_hsblastn
import AlignmentResult
import barleymapcore.utils.alignment_utils as alignment_utils
from barleymapcore.m2p_exception import m2pException

class BaseAligner(object):
    
    _app_path = ""
    _n_threads = 1
    _dbs_path = ""
    _results_hits = []
    _results_unmapped = []
    _verbose = False
    
    def __init__(self, app_path, n_threads, dbs_path, verbose = False):
        self._app_path = app_path
        self._n_threads = n_threads
        self._dbs_path = dbs_path
        self._verbose = verbose
    
    def align(self, fasta_path, dbs_path, db, threshold_id, threshold_cov):
        pass
    
    def get_hits(self):
        return self._results_hits
    
    def get_unmapped(self):
        return self._results_unmapped
    
class SplitBlastnAligner(BaseAligner):
    _split_blast_path = ""
    
    def __init__(self, app_path, n_threads, dbs_path, split_blast_path, verbose = False):
        BaseAligner.__init__(self, app_path, n_threads, dbs_path, verbose)
        self._split_blast_path = split_blast_path
        
    def align(self, fasta_path, db, threshold_id, threshold_cov):
        
        sys.stderr.write("\n")
        
        fasta_headers = alignment_utils.get_fasta_headers(fasta_path)
        
        sys.stderr.write("SplitBlastnAligner: DB --> "+str(db)+"\n")
        sys.stderr.write("SplitBlastnAligner: to align "+str(len(fasta_headers))+"\n")
        
        # get_best_score_hits from m2p_split_blast.py
        self._results_hits = m2p_split_blast.get_best_score_hits(self._split_blast_path, self._app_path, self._n_threads, \
                                                 fasta_path, self._dbs_path, db, threshold_id, threshold_cov, \
                                                 self._verbose)
        
        query_list = [a.get_query_id() for a in self._results_hits]
        
        sys.stderr.write("SplitBlastnAligner: aligned "+str(len(set([a.split(" ")[0] for a in query_list])))+"\n")
        
        self._results_unmapped = alignment_utils.filter_list(fasta_headers, query_list)
        
        sys.stderr.write("SplitBlastnAligner: no hits "+str(len(self._results_unmapped))+"\n")
    
class GMAPAligner(BaseAligner):
    def align(self, fasta_path, db, threshold_id, threshold_cov):
        
        sys.stderr.write("\n")
        
        fasta_headers = alignment_utils.get_fasta_headers(fasta_path)
        
        sys.stderr.write("GMAPAligner: DB --> "+str(db)+"\n")
        sys.stderr.write("GMAPAligner: to align "+str(len(fasta_headers))+"\n")
        
        # get_hits from m2p_gmap.py
        self._results_hits = m2p_gmap.get_hits(self._app_path, self._n_threads, fasta_path, self._dbs_path, db,
                                      threshold_id, threshold_cov, \
                                      self._verbose)
        
        query_list = [a.get_query_id() for a in self._results_hits]
        
        sys.stderr.write("GMAPAligner: aligned "+str(len(set([a.split(" ")[0] for a in query_list])))+"\n")
        
        self._results_unmapped = alignment_utils.filter_list(fasta_headers, query_list)
        
        sys.stderr.write("GMAPAligner: no hits "+str(len(self._results_unmapped))+"\n")

class HSBlastnAligner(BaseAligner):
    
    def __init__(self, app_path, n_threads, dbs_path, verbose = False):
        BaseAligner.__init__(self, app_path, n_threads, dbs_path, verbose)
        
    def align(self, fasta_path, db, threshold_id, threshold_cov):
        
        sys.stderr.write("\n")
        
        fasta_headers = alignment_utils.get_fasta_headers(fasta_path)
        
        sys.stderr.write("HSBlastnAligner: DB --> "+str(db)+"\n")
        sys.stderr.write("HSBlastnAligner: to align "+str(len(fasta_headers))+"\n")
        
        # get_best_score_hits from m2p_hs_blast.py
        self._results_hits = m2p_hsblastn.get_best_score_hits(self._app_path, self._n_threads, fasta_path, self._dbs_path, db, \
                                                 threshold_id, threshold_cov, \
                                                 self._verbose)
        
        query_list = [a.get_query_id() for a in self._results_hits]
        
        sys.stderr.write("HSBlastnAligner: aligned "+str(len(set([a.split(" ")[0] for a in query_list])))+"\n")
        
        self._results_unmapped = alignment_utils.filter_list(fasta_headers, query_list)
        
        sys.stderr.write("HSBlastnAligner: no hits "+str(len(self._results_unmapped))+"\n")

class ListAligner(BaseAligner):
    _aligner_list = []
    _blastn_hits = []
    _gmap_hits = []
    _tmp_files_dir = ""
    
    def __init__(self, aligner_list, tmp_files_dir):
        self._aligner_list = aligner_list
        self._tmp_files_dir = tmp_files_dir
        
    def align(self, fasta_path, db, threshold_id, threshold_cov):
        fasta_to_align = fasta_path
        
        prev_aligner_to_align = fasta_to_align
        fasta_created = False
        
        try:
            for aligner in self._aligner_list:
                if self._verbose: sys.stderr.write("ListAligner: "+str(aligner)+"\n")
                
                try:
                    aligner.align(prev_aligner_to_align, db, threshold_id, threshold_cov)
                except m2pException as m2pe:
                    sys.stderr.write("\t"+m2pe.msg+"\n")
                    sys.stderr.write("\tContinuing with next aligner...\n")
                    continue
                
                prev_aligner_to_align = alignment_utils.extract_fasta_headers(fasta_path, \
                                                                              aligner.get_unmapped(), self._tmp_files_dir)
                fasta_created = True
                
                self._results_hits = self._results_hits + aligner.get_hits()
                self._results_unmapped = aligner.get_unmapped()
                if len(self._results_unmapped) == 0: break # CPCantalapiedra 201701
        
        except Exception:
            raise
        finally:
            if fasta_created: os.remove(prev_aligner_to_align)
 
##