#!/usr/bin/env python
# -*- coding: utf-8 -*-

# GeneEnricher.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from barleymapcore.genes.GenesBase import GenesFields

from barleymapcore.maps.MappingResults import MappingResult

def get_sort_pos_genes(sort_param, map_has_cm, map_has_bp):
    sort_by = -1
    sec_pos = -1
    
    # sort type
    if sort_param == "cm":
        if map_has_cm:
            sort_by = GenesFields.GENES_CM_POS
            sec_pos = GenesFields.GENES_BP_POS
        else:
            sort_by = GenesFields.GENES_BP_POS
            sec_pos = GenesFields.GENES_CM_POS
    elif sort_param == "bp":
        if map_has_bp:
            sort_by = GenesFields.GENES_BP_POS
            sec_pos = GenesFields.GENES_CM_POS
        else:
            sort_by = GenesFields.GENES_CM_POS
            sec_pos = GenesFields.GENES_BP_POS
    else:
        raise Exception("get_sort_pos_genes: Wrong field for sorting "+str(sort_param))
    
    return (sort_by, sec_pos)
    
class GeneEnricher(object):
    
    _enricher = None
    _annotator = None
    _verbose = False
    
    def __init__(self, show_genes_option, annotator, verbose = True):
        if show_genes_option == "marker":
            self._enricher = GenesOnMarker(annotator, verbose)
        elif show_genes_option == "between":
            self._enricher = GenesBetween(annotator, verbose)
        else:
            raise Exception("Unknown option for show genes: "+str(show_genes_option))
        
        self._verbose = verbose
        
    def enrich(self, list_of_genes, sorted_positions, map_sort_by, sort_param, genes_extend, genes_window, load_annot, map_has_cm, map_has_bp):
        return self._enricher.enrich(list_of_genes, sorted_positions, map_sort_by, sort_param, genes_extend, genes_window, load_annot, map_has_cm, map_has_bp)
        
    def _empty_gene(self, load_annot):
        ret_gene = ["empty_gene"]+(GenesFields.GENES_FIELDS - 1)*["-"]
        
        if load_annot: self._annotator.annotate_void(ret_gene)
        
        return ret_gene
    
    def _empty_position(self, chrom):
        return None#["extended"]+[str(chrom)]+(MapPosition.MAP_FIELDS-2)*["-"]
    
    ## This function creates and indexed list of genes
    ## First index is chromosome
    ## Second index is based on position (cM or bp, depending on sorted_by)
    ## Each index position contains a list of genes on that chromosome-position
    def _get_index_of_genes(self, list_of_genes, sorted_by, sec_pos):
        
        if self._verbose: sys.stderr.write("GeneEnricher: creating index...\n")
        
        genes_index = {}
        for gene in list_of_genes:
            gene_id = gene[GenesFields.GENES_ID_POS]
            gene_chr = int(gene[GenesFields.GENES_CHR_POS])
            gene_pos = float(gene[sorted_by])
            
            if gene_chr in genes_index:
                chr_dict = genes_index[gene_chr]
                if gene_pos in chr_dict:
                    list_on_position = chr_dict[gene_pos]
                else:
                    list_on_position = []
                    chr_dict[gene_pos] = list_on_position
                    
            else:
                chr_dict = {}
                list_on_position = []
                chr_dict[gene_pos] = list_on_position
                genes_index[gene_chr] = chr_dict
            
            list_on_position.append(gene)
        
        if self._verbose: sys.stderr.write("GeneEnricher: index ready.\n")
        
        return genes_index
    
    ## This function is like _append_genes
    ## but for empty genes, because "-" position can not be used for sorting :)
    def _append_empty(self, empty_gene, position, ready_positions):
        
        new_position = [a for a in position]
        new_position.extend(empty_gene)
        ready_positions.append(new_position)
        
        return
    
    def _append_genes(self, list_of_genes, position, sorted_by, sec_pos, ready_positions, map_has_cm, map_has_bp):
        
        #if self._verbose: sys.stderr.write("GeneEnricher: appending genes: sorting...\n")
        
        if map_has_bp and map_has_cm:
            sorted_list = sorted(list_of_genes, key=lambda sorting: \
                             (float(sorting[sorted_by]), float(sorting[sec_pos]), sorting[GenesFields.GENES_ID_POS]))
        else:
            sorted_list = sorted(list_of_genes, key=lambda sorting: \
                             (float(sorting[sorted_by]), sorting[GenesFields.GENES_ID_POS]))
        
        #if self._verbose: sys.stderr.write("GeneEnricher: appending genes: parsing...\n")
        
        for gene in sorted_list:
            new_position = [a for a in position]
            new_position.extend(gene)
            ready_positions.append(new_position)
        
        #if self._verbose: sys.stderr.write("GeneEnricher: appending genes: ready.\n")
        
        return
    
    def _get_genes_around(self, index_of_genes, chrom, start_pos, end_pos):
        genes = []
        
        if chrom in index_of_genes:
            chr_dict = index_of_genes[chrom]
            
            for genes_list in [chr_dict[position] for position in chr_dict \
                               if position>=start_pos \
                               and position<=end_pos]:
                genes.extend(genes_list)
        
        return genes
    
    def _get_genes_after(self, index_of_genes, chrom, genes_pos, genes_window):
        genes = []
        
        if chrom in index_of_genes:
            chr_dict = index_of_genes[chrom]
            
            for genes_list in [chr_dict[position] for position in chr_dict \
                               if position>genes_pos and \
                               position<=genes_pos+genes_window]:
                genes.extend(genes_list)
        
        return genes
    
    def _get_genes_before(self, index_of_genes, chrom, genes_pos, genes_window):
        genes = []
       
        if chrom in index_of_genes:
            chr_dict = index_of_genes[chrom]
            
            for genes_list in [chr_dict[position] for position in chr_dict \
                               if position<genes_pos and \
                               position>=genes_pos-genes_window]:
                genes.extend(genes_list)
        
        return genes
    
    def _get_genes_between(self, index_of_genes, chrom, start_pos, end_pos):
        genes = []
        
        if chrom in index_of_genes:
            chr_dict = index_of_genes[chrom]
            
            for genes_list in [chr_dict[position] for position in chr_dict \
                               if position>=start_pos and \
                               position<end_pos]:
                genes.extend(genes_list)
        
        return genes

class GenesBetween(GeneEnricher):
    
    def __init__(self, annotator, verbose):
        self._annotator = annotator
        self._verbose = verbose
    
    # Obtain genes around a Marker
    def enrich(self, list_of_genes, sorted_positions, map_sort_by, sort_param, genes_extend, genes_window, load_annot, map_has_cm, map_has_bp):
        ready_positions = []
        
        if self._verbose: sys.stderr.write("GenesBetween: sort by: "+str(sort_param)+"\n")
        
        (genes_sort_by, genes_sec_pos) = get_sort_pos_genes(sort_param, map_has_cm, map_has_bp)
        
        index_of_genes = self._get_index_of_genes(list_of_genes, genes_sort_by, genes_sec_pos)
        
        prev_position = None
        prev_chrom = -1
        prev_genes_pos = -1
        for position in sorted_positions:
            chrom = position.get_chrom()
            genes_pos = float(position[map_sort_by])
            
            if chrom != prev_chrom: # New chromosome
                
                # Append genes of the last position on the previous chromosome
                if prev_position:
                    
                    if genes_extend:
                        genes_after = self._get_genes_after(index_of_genes, prev_chrom, prev_genes_pos, genes_window)
                        new_position = self._empty_position(prev_chrom)
                    else:
                        genes_after = []
                        new_position = prev_position
                    
                    genes_around = self._get_genes_around(index_of_genes, prev_chrom, prev_genes_pos, prev_genes_pos)
                    
                    if len(genes_around) == 0:
                        self._append_empty(self._empty_gene(load_annot), prev_position, ready_positions)
                    else:
                        self._append_genes(genes_around, prev_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
                    
                    if len(genes_after) != 0:
                        self._append_genes(genes_after, new_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
                
                # Insert genes before current position of the new chromosome
                if genes_extend:
                    genes_before = self._get_genes_before(index_of_genes, chrom, genes_pos, genes_window)
                    
                    new_position = self._empty_position(chrom)
                    if len(genes_before) == 0:
                        self._append_empty(self._empty_gene(load_annot), new_position, ready_positions)
                    else:
                        self._append_genes(genes_before, new_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
            
            elif genes_pos != prev_genes_pos:
                # Append genes between the previous position and this one
                genes_between = self._get_genes_between(index_of_genes, chrom, prev_genes_pos, genes_pos)
                
                if len(genes_between) == 0:
                    self._append_empty(self._empty_gene(load_annot), prev_position, ready_positions)
                else:
                    self._append_genes(genes_between, prev_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
                
            else:
                self._append_empty(self._empty_gene(load_annot), prev_position, ready_positions)
            
            prev_position = position
            prev_chrom = chrom
            prev_genes_pos = genes_pos
        
        # Last marker of all
        if prev_position: # ie.: if there was at least one position in sorted_positions
            if genes_extend:
                genes_after = self._get_genes_after(index_of_genes, prev_chrom, prev_genes_pos, genes_window)
                new_position = self._empty_position(chrom)
            else:
                genes_after = []
                new_position = prev_position
            
            genes_around = self._get_genes_around(index_of_genes, prev_chrom, prev_genes_pos, prev_genes_pos)
            
            if len(genes_around) == 0:
                self._append_empty(self._empty_gene(load_annot), prev_position, ready_positions)
            else:
                self._append_genes(genes_around, prev_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
            
            if len(genes_after) != 0:
                self._append_genes(genes_after, new_position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
        
        return ready_positions

class GenesOnMarker(GeneEnricher):
    
    def __init__(self, annotator, verbose):
        self._annotator = annotator
        self._verbose = verbose
    
    # Obtain genes around a Marker
    def enrich(self, list_of_genes, sorted_positions, map_sort_by, sort_param, genes_extend, genes_window, load_annot, map_has_cm, map_has_bp):
        ready_positions = []
        
        if self._verbose: sys.stderr.write("GenesOnMarker: sort by: "+str(sort_param)+"\n")
        
        (genes_sort_by, genes_sec_pos) = get_sort_pos_genes(sort_param, map_has_cm, map_has_bp)
        
        index_of_genes = self._get_index_of_genes(list_of_genes, genes_sort_by, genes_sec_pos)
        
        for position in sorted_positions:
            #if self._verbose: sys.stderr.write("GenesOnMarker: pos: "+str(position)+"\n")
            chrom = position.get_chrom()
            
            genes_pos = float(position[map_sort_by])
            
            if genes_extend:
                start_pos = genes_pos - genes_window
                end_pos = genes_pos + genes_window
            else:
                start_pos = genes_pos
                end_pos = genes_pos
            
            genes = self._get_genes_around(index_of_genes, chrom, start_pos, end_pos)
            
            if len(genes) == 0:
                self._append_empty(self._empty_gene(load_annot), position, ready_positions)
            else:
                self._append_genes(genes, position, genes_sort_by, genes_sec_pos, ready_positions, map_has_cm, map_has_bp)
        
        return ready_positions

## END