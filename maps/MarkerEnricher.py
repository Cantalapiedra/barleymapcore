#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MarkerEnricher.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from MapsBase import MapFile, MapFields

def get_sort_pos_contigs(sort_param, map_has_cm, map_has_bp):
    sort_by = -1
    sec_pos = -1
    
    # sort type
    if sort_param == "cm":
        if map_has_cm:
            sort_by = MapFile.MAP_FILE_CM
            sec_pos = MapFile.MAP_FILE_BP
        else:
            sort_by = MapFile.MAP_FILE_BP
            sec_pos = MapFile.MAP_FILE_CM
    elif sort_param == "bp":
        if map_has_bp:
            sort_by = MapFile.MAP_FILE_BP
            sec_pos = MapFile.MAP_FILE_CM
        else:
            sort_by = MapFile.MAP_FILE_CM
            sec_pos = MapFile.MAP_FILE_BP
    else:
        raise Exception("get_sort_pos_contigs: Wrong field for sorting "+str(sort_param))
    
    return (sort_by, sec_pos)

## All these functions are very similar to those in GeneEnricher.GenesBetween
class MarkerEnricher(object):
    
    def __init__(self):
        return
    
    def _get_contigs_before(self, positions_index, pos_chr, pos_pos, genes_window):
        contigs = []
        
        if pos_chr in positions_index:
            chr_dict = positions_index[pos_chr]
            
            for contigs_list in [chr_dict[position] for position in chr_dict \
                                   if position<pos_pos and \
                                   position>=pos_pos-genes_window]:
                    contigs.extend(contigs_list)  
        
        return contigs
    
    def _get_contigs_around(self, index_of_genes, chrom, start_pos, end_pos):
        contigs = []
        
        if chrom in index_of_genes:
            chr_dict = index_of_genes[chrom]
            
            for contigs_list in [chr_dict[position] for position in chr_dict \
                               if position>=start_pos \
                               and position<=end_pos]:
                contigs.extend(contigs_list)
        
        return contigs
    
    def _get_contigs_after(self, positions_index, pos_chr, pos_pos, genes_window):
        contigs = []
        
        if pos_chr in positions_index:
            chr_dict = positions_index[pos_chr]
            
            for contigs_list in [chr_dict[position] for position in chr_dict \
                               if position>pos_pos and \
                               position<=pos_pos+genes_window]:
                contigs.extend(contigs_list)
        
        return contigs
    
    def _get_contigs_between(self, index_of_genes, chrom, start_pos, end_pos):
        
        contigs = []
        
        if chrom in index_of_genes:
            
            chr_dict = index_of_genes[chrom]
            
            for contigs_list in [chr_dict[position] for position in chr_dict \
                               if position>=start_pos and \
                               position<end_pos]:
                
                contigs.extend(contigs_list)
        
        #sys.stderr.write("contigs "+str(len(contigs))+"\n")
        
        return contigs
    
    def _empty_position(self, chrom):
        return ["extended"]+[str(chrom)]+(MapFields.MAP_FIELDS-2)*["-"]
    
    def _empty_contig(self, map_has_cm, map_has_bp):
        if map_has_bp and map_has_cm:
            ret_gene = ["empty_marker", "-", "-", "-"]
        else:
            ret_gene = ["empty_marker", "-", "-"]
        
        return ret_gene
    
    def _append_empty(self, empty_contig, position, ready_positions):
        
        new_position = [a for a in position]
        new_position.extend(empty_contig)
        ready_positions.append(new_position)
        
        return
    
    def _append_contigs(self, list_of_contigs, position, sorted_by, sec_pos, ready_positions, map_has_cm, map_has_bp):
        
        if map_has_bp and map_has_cm:
            sorted_list = sorted(list_of_contigs, key=lambda sorting: \
                             (float(sorting[sorted_by]), float(sorting[sec_pos]), sorting[MapFile.MAP_FILE_MARKER]))
        else:
            sorted_list = sorted(list_of_contigs, key=lambda sorting: \
                             (float(sorting[sorted_by]), sorting[MapFile.MAP_FILE_MARKER]))
        
        for contig in sorted_list:
            new_position = [a for a in position]
            new_position.extend(contig)
            ready_positions.append(new_position)
        
        return
    
    def enrich(self, positions_index, genetic_map_has_cm_pos, genetic_map_has_bp_pos, sorted_positions, map_sort_by, sort_param, genes_extend, genes_window):
        
        new_positions = []
        
        (contigs_sort_by, contigs_sec_pos) = get_sort_pos_contigs(sort_param, genetic_map_has_cm_pos, genetic_map_has_bp_pos)
        
        contigs_found = 0
        
        # For each two consecutive positions in the same chromosome in the map, find the contigs between positions
        prev_position = None
        prev_chr = -1
        prev_pos = -1
        for position in sorted_positions:
            
            #sys.stderr.write("Positions to compare "+str(prev_position)+" - "+str(position)+"\n")
            
            pos_chr = int(position[MapFields.MARKER_CHR_POS])
            pos_pos = float(position[map_sort_by])
            
            if prev_chr != pos_chr:
                
                # Append extension for previous position and insert it for next first position
                if prev_position:
                    
                    prev_chr = int(prev_position[MapFields.MARKER_CHR_POS])
                    prev_pos = long(prev_position[map_sort_by])
                    
                    if genes_extend:
                        contigs_after = self._get_contigs_after(positions_index, prev_chr, prev_pos, genes_window)
                        new_position = self._empty_position(prev_chr)
                    else:
                        contigs_after = []
                        new_position = prev_position
                    
                    contigs_around = self._get_contigs_around(positions_index, prev_chr, prev_pos, prev_pos)
                    
                    if len(contigs_around) == 0:
                        self._append_empty(self._empty_contig(genetic_map_has_cm_pos, genetic_map_has_bp_pos), prev_position, new_positions)
                    else:
                        self._append_contigs(contigs_around, prev_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                             genetic_map_has_cm_pos, genetic_map_has_bp_pos)
                    
                    if len(contigs_after) != 0:
                        self._append_contigs(contigs_after, new_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                             genetic_map_has_cm_pos, genetic_map_has_bp_pos)
                    
                    contigs_found += len(contigs_around)+len(contigs_after)
                    
                # Insert extension for first position
                if genes_extend:
                    contigs_before = self._get_contigs_before(positions_index, pos_chr, pos_pos, genes_window)
                    
                    new_position = self._empty_position(pos_chr)
                    if len(contigs_before) == 0:
                        self._append_empty(self._empty_contig(genetic_map_has_cm_pos, genetic_map_has_bp_pos), new_position, new_positions)
                    else:
                        self._append_contigs(contigs_before, new_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                             genetic_map_has_cm_pos, genetic_map_has_bp_pos)
                    
                    contigs_found += len(contigs_before)
                    
            elif pos_pos != prev_pos:
                
                if prev_pos > pos_pos: raise Exception("Map positions should have been already sorted, but they are not.")
                
                # Append genes between the previous position and this one
                contigs_between = self._get_contigs_between(positions_index, pos_chr, prev_pos, pos_pos)
                
                if len(contigs_between) == 0:
                    self._append_empty(self._empty_contig(genetic_map_has_cm_pos, genetic_map_has_bp_pos), prev_position, new_positions)
                else:
                    self._append_contigs(contigs_between, prev_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                         genetic_map_has_cm_pos, genetic_map_has_bp_pos)
                
                contigs_found += len(contigs_between)
                
            else:
                self._append_empty(self._empty_contig(genetic_map_has_cm_pos, genetic_map_has_bp_pos), prev_position, new_positions)
            
            prev_position = position
            prev_chr = pos_chr
            prev_pos = pos_pos
        
        # Append extension for last position
        if prev_position: # ie.: if there was at least one position in sorted_positions
            if genes_extend:
                contigs_after = self._get_contigs_after(positions_index, prev_chr, prev_pos, genes_window)
                new_position = self._empty_position(prev_chr)
            else:
                contigs_after = []
                new_position = prev_position
            
            contigs_around = self._get_contigs_around(positions_index, prev_chr, prev_pos, prev_pos)
            
            if len(contigs_around) == 0:
                self._append_empty(self._empty_contig(genetic_map_has_cm_pos, genetic_map_has_bp_pos), prev_position, new_positions)
            else:
                self._append_contigs(contigs_around, prev_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                     genetic_map_has_cm_pos, genetic_map_has_bp_pos)
            
            if len(contigs_after) != 0:
                self._append_contigs(contigs_after, new_position, contigs_sort_by, contigs_sec_pos, new_positions, \
                                     genetic_map_has_cm_pos, genetic_map_has_bp_pos)
            contigs_found += len(contigs_around)+len(contigs_after)
        
        #sys.stderr.write("Contigs found "+str(contigs_found)+"\n")
        
        return new_positions
    
## END