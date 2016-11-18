#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aligners.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import os, sys

SELECTION_BEST_SCORE = "best_score"
SELECTION_NONE = "none"

from m2p_split_blast import get_best_score_hits
from m2p_gmap import get_hits
import barleymapcore.utils.alignment_utils as alignment_utils

class AlignmentResults(object):
    QUERY_ID = 0
    SUBJECT_ID = 1
    ALIGN_IDENTITY = 2
    QUERY_COVERAGE = 3
    ALIGNMENT_SCORE = 4
    STRAND = 5
    LOCAL_POSITION = 6
    DB_NAME = 7
    ALGORITHM = 8

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
    
    def align(self, fasta_path, dbs_path, db, threshold_id, threshold_cov, selection):
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
        
    def align(self, fasta_path, db, threshold_id, threshold_cov, selection):
        fasta_headers = alignment_utils.get_fasta_headers(fasta_path)
        sys.stderr.write("SplitBlastnAligner: num of seqs to align to "+str(db)+" "+str(len(fasta_headers))+"\n")
        
        # get_best_score_hits from m2p_split_blast.py
        self._results_hits = get_best_score_hits(self._split_blast_path, self._app_path, self._n_threads, \
                                                 fasta_path, self._dbs_path, db, threshold_id, threshold_cov, \
                                                 selection, self._verbose)
        
        self._results_unmapped = alignment_utils.filter_list(fasta_headers, [a[0] for a in self._results_hits])
        
        sys.stderr.write("SplitBlastnAligner: num of seqs unmapped "+str(len(self._results_unmapped))+"\n")
    
class GMAPAligner(BaseAligner):
    def align(self, fasta_path, db, threshold_id, threshold_cov, selection):
        fasta_headers = alignment_utils.get_fasta_headers(fasta_path)
        
        sys.stderr.write("GMAPAligner: num of seqs to align to "+str(db)+" "+str(len(fasta_headers))+"\n")
        
        # get_hits from m2p_gmap.py
        self._results_hits = get_hits(self._app_path, self._n_threads, fasta_path, self._dbs_path, db, threshold_id, threshold_cov, \
                                      selection, self._verbose)
        
        #sys.stderr.write("Results: "+str(len([a[0] for a in self._results_hits]))+"\n")
        #sys.stderr.write("ResultsB: "+str(len(set([a.split(" ")[0] for a in [a[0] for a in self._results_hits]])))+"\n")
        
        self._results_unmapped = alignment_utils.filter_list(fasta_headers, [a[0] for a in self._results_hits])
        
        sys.stderr.write("GMAPAligner: num of seqs aligned "+str(len(set([a.split(" ")[0] for a in [a[0] for a in self._results_hits]])))+"\n")
        sys.stderr.write("GMAPAligner: num of seqs without hits "+str(len(self._results_unmapped))+"\n")

class DualAligner(BaseAligner):
    _blastn_aligner = None
    _gmap_aligner = None
    _blastn_hits = []
    _gmap_hits = []
    _tmp_files_dir = ""
    
    def __init__(self, blastn_aligner, gmap_aligner, tmp_files_dir):
        self._blastn_aligner = blastn_aligner
        self._gmap_aligner = gmap_aligner
        self._tmp_files_dir = tmp_files_dir
        
    def align(self, fasta_path, db, threshold_id, threshold_cov, selection):
        fasta_to_align = fasta_path
        self._blastn_aligner.align(fasta_to_align, db, threshold_id, threshold_cov, selection)
        
        gmap_fasta_to_align = alignment_utils.extract_fasta_headers(fasta_path, self._blastn_aligner.get_unmapped(), self._tmp_files_dir)
        try:
            self._gmap_aligner.align(gmap_fasta_to_align, db, threshold_id, threshold_cov, selection)
        except Exception:
            raise
        finally:
            os.remove(gmap_fasta_to_align)
            
        self._gmap_hits = self._gmap_aligner.get_hits()
        
        self._results_hits = self._blastn_aligner.get_hits()+self._gmap_hits
        self._results_unmapped = self.get_gmap_unmapped()
        
    def get_blastn_hits(self):
        return self._blastn_aligner.get_hits()
    
    def get_blastn_unmapped(self):
        return self._blastn_aligner.get_unmapped()
    
    def get_gmap_hits(self):
        return self._gmap_aligner.get_hits()
    
    def get_gmap_unmapped(self):
        return self._gmap_aligner.get_unmapped()

class ListAligner(BaseAligner):
    _aligner_list = []
    _blastn_hits = []
    _gmap_hits = []
    _tmp_files_dir = ""
    
    def __init__(self, aligner_list, tmp_files_dir):
        self._aligner_list = aligner_list
        self._tmp_files_dir = tmp_files_dir
        
    def align(self, fasta_path, db, threshold_id, threshold_cov, selection):
        fasta_to_align = fasta_path
        
        prev_aligner_to_align = fasta_to_align
        fasta_created = False
        
        try:
            for aligner in self._aligner_list:
                #sys.stderr.write("ListAligner: "+str(aligner)+"\n")
                aligner.align(prev_aligner_to_align, db, threshold_id, threshold_cov, selection)
                
                prev_aligner_to_align = alignment_utils.extract_fasta_headers(fasta_path, \
                                                                              aligner.get_unmapped(), self._tmp_files_dir)
                fasta_created = True
                
                self._results_hits = self._results_hits + aligner.get_hits()
                self._results_unmapped = aligner.get_unmapped()
                
        except Exception:
            raise
        finally:
            
            if fasta_created: os.remove(prev_aligner_to_align)
    
##