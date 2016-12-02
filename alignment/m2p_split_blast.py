#!/usr/bin/env python
# -*- coding: utf-8 -*-

# m2p_split_blast.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from subprocess import Popen, PIPE

#from Aligners import SELECTION_BEST_SCORE, SELECTION_NONE

def __split_blast(split_blast_path, blast_app_path, n_threads, query_fasta_path, blast_dbs_path, db_name, verbose = False):
    results = []
    
    ###### Retrieve num of fasta seqs to calculate necessary bins
    retValue = 0
    p = Popen(" ".join(["cat", query_fasta_path, " | grep -c \"^>\""]), \
              shell=True, stdout=PIPE, stderr=sys.stderr)
    output = p.communicate()[0]
    retValue = p.returncode
    if retValue == 0:
        num_of_seqs = int(output)
        split_blast_bins = (num_of_seqs / 50) + 1
    else:
        split_blast_bins = 100
    
    ###### Split blast
    blast_command = " ".join([split_blast_path+"split_blast.pl", str(n_threads), str(split_blast_bins), blast_app_path, \
                "-dust no -soft_masking false -task megablast", \
                '-outfmt \\"6 qseqid qlen sseqid slen length qstart qend sstart send bitscore evalue pident mismatch gapopen\\"'])
    
    blast_db = "".join(["-db ", blast_dbs_path, db_name]) # blast_db = "".join(["-db ", blast_dbs_path, db_name , ".fa"]) # 
    blast_query = " ".join(["-query ", query_fasta_path])
    #blast_cmd = " ".join([ResourcesMng.get_deploy_dir()+blast_command, blast_db, blast_query])
    blast_cmd = " ".join([blast_command, blast_db, blast_query])
    
    if verbose: sys.stderr.write("m2p_split_blast: Executing '"+blast_cmd+"'\n")
    
    retValue = 0
    p = Popen(blast_cmd, shell=True, stdout=PIPE, stderr=sys.stderr)
    com_list = p.communicate()
    output = com_list[0]
    output_err = com_list[1]
    retValue = p.returncode
    
    if retValue != 0: raise Exception("m2p_split_blast: Blast return != 0. "+blast_cmd+"\n"+str(output)+"\n"+str(output_err)+"\n")
    if "error" in output or "Error" in output or "ERROR" in output:
        sys.stderr.write("m2p_split_blast: error in blast output. We will report 0 results for this alignment.\n")
        sys.stderr.write(output+"\n")
        sys.stderr.write(str(output_err)+"\n")
        results = []
    else:
        if verbose: sys.stderr.write("m2p_split_blast: Blast return value "+str(retValue)+"\n")
        
        [results.append(line) for line in output.strip().split("\n") if line != "" and not line.startswith("#")]
        # startswith("#") due to split_blast.pl printing in stdout comments, warnings and so on
    
    return results

def __filter_blast_results(results, threshold_id, threshold_cov, db_name, verbose = False):
    
    filtered_results = []
    
    # assert: if an int is passed to the method, the comparison fails
    thres_id = float(threshold_id)
    thres_cov = float(threshold_cov)
    
    filter_dict = {}
    
    for line in results:
        line_data = line.split("\t")
        
        # filter: based on identity
        align_ident = float(line_data[11])
        if align_ident < thres_id:
            continue
        
        # filter: based on query coverage of alignment
        align_len = int(line_data[4])
        query_len = int(line_data[1])
        
        query_cov = (align_len/float(query_len))*100
        if query_cov < thres_cov:
            continue
        
        query_id = line_data[0]
        subject_id = line_data[2]
        align_score = float(line_data[9])
        
        if query_id == "2249-1308": print "FOUND "+subject_id
        
        # strand and local position
        if line_data[7]>line_data[8]:
            strand = "-"
            local_position = long(line_data[8])
            end_position = long(line_data[7])
        else:
            strand = "+"
            local_position = long(line_data[7])
            end_position = long(line_data[8])
        
        qstart_pos = long(line_data[5])
        qend_pos = long(line_data[6])
        
        result_tuple = (subject_id, align_ident, query_cov, align_score,
                        strand, local_position, end_position, qstart_pos, qend_pos)
        
        # For a given DB, keep always the best score
        #if selection == SELECTION_BEST_SCORE:
        if query_id in filter_dict:
            prev_max_score = filter_dict[query_id]["max_score"]
            
            if align_score > prev_max_score:
                filter_dict[query_id]["query_list"] = [result_tuple]
                filter_dict[query_id]["max_score"] = align_score
                
            elif align_score == prev_max_score:
                filter_dict[query_id]["query_list"].append(result_tuple)
            #else:
            #    print "FILTERED"
        else:
            filter_dict[query_id] = {"query_list":[result_tuple], "max_score":align_score}
        
        #elif selection == SELECTION_NONE:
        #    if query_id in filter_dict:
        #        filter_dict[query_id]["query_list"].append(result_tuple)
        #    else:
        #        filter_dict[query_id] = {"query_list":[result_tuple], "max_score":-1}
        #else:
        #    raise Exception("m2p_split_blast: unknown value "+str(selection)+" for selection parameter.")
    
    algorithm = "blastn"
    # Recover filtered results
    for query_id in filter_dict:
        for alignment_data in filter_dict[query_id]["query_list"]:
            subject_id = alignment_data[0]
            align_ident = alignment_data[1]
            query_cov = alignment_data[2]
            align_score = alignment_data[3]
            strand = alignment_data[4]
            local_position = alignment_data[5]
            end_position = alignment_data[6]
            qstart_pos = alignment_data[7]
            qend_pos = alignment_data[8]
            # This MUST coincide with Aligners.AlignmentResults fields
            filtered_results.append([query_id, subject_id, align_ident, query_cov, align_score, strand,
                                     qstart_pos, qend_pos, local_position, end_position,  
                                     db_name, algorithm])
    
    return filtered_results

def get_best_score_hits(split_blast_path, blast_app_path, n_threads, query_fasta_path, blast_dbs_path, db_name, \
                        threshold_id, threshold_cov, verbose = False):
    
    if verbose: sys.stderr.write("m2p_split_blast: "+query_fasta_path+" against "+db_name+"\n")
    
    results = __split_blast(split_blast_path, blast_app_path, n_threads, query_fasta_path, blast_dbs_path, db_name, verbose)
    
    if verbose: sys.stderr.write("m2p_split_blast: raw results --> "+str(len(results))+"\n")
    
    if len(results)>0:
        results = __filter_blast_results(results, threshold_id, threshold_cov, db_name, verbose)
        
    if verbose: sys.stderr.write("m2p_split_blast: pass-filter results --> "+str(len(results))+"\n")
    #sys.stderr.write(str(len(results))+"\n")
    #sys.stderr.write(str(results)+"\n")
    
    return results