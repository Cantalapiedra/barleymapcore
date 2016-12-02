#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AlignmentFacade.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import os, sys

from Aligners import AlignmentResults, SplitBlastnAligner, GMAPAligner, DualAligner, ListAligner
#, SELECTION_BEST_SCORE
import barleymapcore.utils.alignment_utils as alignment_utils
from barleymapcore.db.DatabasesConfig import REF_TYPE_BIG, REF_TYPE_STD

QUERY_MODE_GENOMIC = "genomic"
QUERY_MODE_CDNA = "cdna"
QUERY_MODE_DUAL = "auto"

class AlignmentFacade():
    
    _split_blast_path = ""
    _blastn_app_path = ""
    _gmap_app_path = ""
    _blastn_dbs_path = ""
    _gmap_dbs_path = ""
    _gmapl_app_path = ""
    _tmp_files_dir = ""
    _databases_config = None
    
    _alignment_results = {}
    _alignment_unmapped = []
    
    _verbose = False
    
    def __init__(self, split_blast_path, blastn_app_path, gmap_app_path,
                 blastn_dbs_path, gmap_dbs_path, gmapl_app_path, 
                 tmp_files_dir, databases_config, verbose = False):
        self._split_blast_path = split_blast_path
        self._blastn_app_path = blastn_app_path
        self._gmap_app_path = gmap_app_path
        self._blastn_dbs_path = blastn_dbs_path
        self._gmap_dbs_path = gmap_dbs_path
        self._gmapl_app_path = gmapl_app_path
        self._tmp_files_dir = tmp_files_dir
        self._verbose = verbose
        
        self._databases_config = databases_config
    
    # Performs the alignment of fasta sequences different DBs
    def perform_alignment(self, query_fasta_path, dbs_list, hierarchical, query_type = "genomic", \
                          threshold_id = 98, threshold_cov = 95, n_threads = 1, \
                          best_score_filter = False, ref_type_param = REF_TYPE_STD):
        results = {} # A list of hits for each db
        
        fasta_to_align = query_fasta_path
        
        # Create a record for each DB
        for db in dbs_list:
            results[db] = []
        
        tmp_files_list = []
        try:
            for db in dbs_list:
                # CPCantalapiedra 2016-11
                # Obtain ref_type of current database
                
                if self._databases_config.database_exists(db):
                    database_config = self._databases_config.get_database(db)
                    ref_type = self._databases_config.get_ref_type(database_config)
                else:
                    ref_type = ref_type_param
                
                # CPCantalapiedra 2016-11
                # Obtain suitable aligner for each database
                aligner = self._get_aligner(query_type, n_threads, self._tmp_files_dir, ref_type)
                
                aligner.align(fasta_to_align, db, threshold_id, threshold_cov)
                
                results[db] = aligner.get_hits()
                
                ## Recover unmapped queries if needed
                if len(dbs_list) > 1 and hierarchical:
                    unmapped = aligner.get_unmapped()
                    
                    if len(unmapped) > 0:
                        fasta_to_align = alignment_utils.extract_fasta_headers(fasta_to_align, unmapped, self._tmp_files_dir)
                        tmp_files_list.append(fasta_to_align)
                    else:
                        break
                # else: fasta_to_align = fasta_path
            
            #unmapped_number = len(alignment_utils.filter_list(fasta_headers, [a[0] for a in self._results_hits]))
            #sys.stderr.write("AlignmentFacade: total number of seqs unmapped "+str(unmapped_number)+"\n")
            
        except Exception:
            raise
        finally:
            for tmp_file in tmp_files_list:
                os.remove(tmp_file)
        
        results = self._best_score(results, best_score_filter)
        
        self._alignment_results = results
        
        ## Recover unmapped queries
        if hierarchical:
            self._alignment_unmapped = aligner.get_unmapped()
        else:
            fasta_headers = alignment_utils.get_fasta_headers(query_fasta_path)
            no_redundant_results = set()
            for db in self._alignment_results:
                for result in self._alignment_results[db]:
                    if result[AlignmentResults.QUERY_ID] not in no_redundant_results:
                        no_redundant_results.add(result[AlignmentResults.QUERY_ID])
                
            self._alignment_unmapped = alignment_utils.filter_list(fasta_headers, no_redundant_results)
        
        return results
    
    def _best_score(self, results, best_score_filter):
        ###### best score filtering
        if best_score_filter:
            best_score_filtering = {}
            for db in results:
                for result in results[db]:
                    query_id = result[AlignmentResults.QUERY_ID]
                    align_score = float(result[AlignmentResults.ALIGNMENT_SCORE])
                    
                    if query_id in best_score_filtering:
                        query_best_score = best_score_filtering[query_id]["best_score"]
                        if align_score < query_best_score:
                            continue
                        elif align_score == query_best_score:
                            best_score_filtering[query_id]["results"].append(result)
                        else: # align_score > query_best_score
                            best_score_filtering[query_id]["results"] = [result]
                            best_score_filtering[query_id]["best_score"] = align_score
                    else:
                        best_score_filtering[query_id] = {"results":[result], "best_score":align_score}
                
                results[db] = []
            
            #results = {}
            #
            ## Create a record for each DB
            #for db in dbs_list:
            #    results[db] = []
            
            for query_id in best_score_filtering:
                for result in best_score_filtering[query_id]["results"]:
                    db = result[AlignmentResults.DB_NAME]
                    results[db].append(result)
        
        return results
    
    # Returns a new aligner based on the query_type supplied
    def _get_aligner(self, query_type, n_threads, tmp_files_dir = "./", ref_type = REF_TYPE_STD): # This is an AlignerFactory
        
        aligner = None
        
        if query_type == QUERY_MODE_GENOMIC:
            aligner = SplitBlastnAligner(self._blastn_app_path, n_threads, self._blastn_dbs_path, self._split_blast_path, self._verbose)
            
        elif query_type == QUERY_MODE_CDNA:
            # CPCantalapiedra 2016-11
            if ref_type == REF_TYPE_BIG: # When sequence DB is too big to use gmap, instead gmapl has to be used
                aligner = GMAPAligner(self._gmapl_app_path, n_threads, self._gmap_dbs_path, self._verbose)
            else:
                aligner = GMAPAligner(self._gmap_app_path, n_threads, self._gmap_dbs_path, self._verbose)
                
        elif query_type == QUERY_MODE_DUAL:
            blastn_aligner = SplitBlastnAligner(self._blastn_app_path, n_threads, self._blastn_dbs_path, self._split_blast_path, self._verbose)
            # CPCantalapiedra 2016-11
            if ref_type == REF_TYPE_BIG:
                gmap_aligner = GMAPAligner(self._gmapl_app_path, n_threads, self._gmap_dbs_path, self._verbose)
            else:
                gmap_aligner = GMAPAligner(self._gmap_app_path, n_threads, self._gmap_dbs_path, self._verbose)
                
            aligner = DualAligner(blastn_aligner, gmap_aligner, tmp_files_dir)
            
        elif "," in query_type:
            aligner_list = []
            for aligner_type in query_type.split(","):
                if self._verbose: sys.stderr.write("AlignmentFacade: aligner "+str(aligner_type)+"\n")
                if aligner_type == QUERY_MODE_GENOMIC:
                    current_aligner = SplitBlastnAligner(self._blastn_app_path, n_threads, self._blastn_dbs_path, self._split_blast_path, self._verbose)
                    aligner_list.append(current_aligner)
                elif aligner_type == QUERY_MODE_CDNA:
                    # CPCantalapiedra 2016-11
                    if ref_type == REF_TYPE_BIG:
                        current_aligner = GMAPAligner(self._gmapl_app_path, n_threads, self._gmap_dbs_path, self._verbose)
                    else:
                        current_aligner = GMAPAligner(self._gmap_app_path, n_threads, self._gmap_dbs_path, self._verbose)
                    aligner_list.append(current_aligner)
                else:
                    sys.stderr.write("WARNING: AlignmentFacade: Unknown aligner type "+aligner_type+".\nSkipping to next aligner.\n")
                
            aligner = ListAligner(aligner_list, tmp_files_dir)
            
        else:
            raise Exception("Unknown query type "+query_type+" when requesting aligner.")
        
        return aligner
    
    def get_alignment_results(self):
        return self._alignment_results
    
    def get_alignment_unmapped(self):
        return self._alignment_unmapped
    
    
    
