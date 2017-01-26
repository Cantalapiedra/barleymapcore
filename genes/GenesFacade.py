#!/usr/bin/env python
# -*- coding: utf-8 -*-

# GenesFacade.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys
from GeneAnnotator import GeneAnnotator
from barleymapcore.maps.enrichment.GeneEnricher import GeneEnricher, get_sort_pos_genes

from barleymapcore.maps.MapsBase import MapTypes
from GenesBase import AnnotFile

class GenesFacade(object):

    _config_path_dict = {}
    _genes_path = ""
    _annot_path = ""
    _annot_filename = ""
    _genes_dict = {}
    _annot_dict = {}
    _verbose = False
    
    _sort_by = -1
    _sec_pos = -1
    
    #_genetic_map_dict = {}
    _annotator = GeneAnnotator()
    _geneEnricher = None
    
    def __init__(self, config_path_dict, load_annot, show_genes_option, verbose = True):
        
        self._config_path_dict = config_path_dict
        
        genes_path = config_path_dict["app_path"]+config_path_dict["genes_path"]
        annot_path = config_path_dict["app_path"]+config_path_dict["annot_path"]
        annot_filename = config_path_dict["annot_filename"]
        
        self._genes_path = genes_path
        self._annot_path = annot_path
        self._annot_filename = annot_filename
        
        #self._genetic_map_dict = genetic_map_dict
        
        if load_annot:
            self.load_annot(self._annot_path, self._annot_filename)
        
        self._geneEnricher = GeneEnricher(show_genes_option, self._annotator, self._verbose)
        
        self._verbose = verbose
        
        return
    
    def enrich_by_pos(self, genetic_map, genetic_map_positions, genes_extend, genes_window, sort_by, sort_param, load_annot, \
                      map_has_cm_pos, map_has_bp_pos):
        
        if self._verbose: sys.stderr.write("GenesFacade: enrich by pos...\n")
        
        enriched_positions = self._geneEnricher.enrich(self._genes_dict[genetic_map], genetic_map_positions, sort_by, sort_param, \
                                                genes_extend, genes_window, load_annot, \
                                                map_has_cm_pos, map_has_bp_pos)
        
        return enriched_positions
    
    def load_data(self, genetic_map, load_annot):
        
        if self._verbose: sys.stderr.write("GenesFacade: loading data...\n")
        
        ### LOAD DATA
        self.load_genes(genetic_map, self._genes_path)
        
        if load_annot:
            self.append_annot(genetic_map, self._annot_dict)
        
        if self._verbose: sys.stderr.write("GenesFacade: data loaded.\n")
        
        return
    
    def load_genes(self, genetic_map, genes_path):
        
        genes_list = []
        
        # Load genes positions associated to this genetic map
        genes_data_path = genes_path+genetic_map+"_genes.tab"
        for line in open(genes_data_path):
            genes_list.append([a.strip() for a in line.strip().split("\t")])
        
        #map_has_bp = self._genetic_map_dict[genetic_map][MapTypes.MAP_HAS_BP_POS]
        
        ###### IS THIS REALLY NECESSARY??
        # Sorts the list of genes for this map
        #self._genes_dict[genetic_map] = self._sort_genes_list(genes_list, sort_param, map_has_bp)
        ## Instead I try just this:
        self._genes_dict[genetic_map] = genes_list
        
        if self._verbose: sys.stderr.write("\t# of genes "+str(len(genes_list))+"\n")
        
        return
    
    #def _sort_genes_list(self, positions_list, sort_param, map_has_bp):
    #    sorted_list = []
    #    
    #    (sort_by, sec_pos) = get_sort_pos_genes(sort_param, map_has_bp)
    #    
    #    if map_has_bp:
    #        sorted_list = sorted(positions_list, key=lambda sorting: (int(sorting[GenesFields.GENES_CHR_POS]), \
    #                                                                  float(sorting[sort_by]), float(sorting[sec_pos]), sorting[GenesFields.GENES_ID_POS]))
    #    else:
    #        sorted_list = sorted(positions_list, key=lambda sorting: (int(sorting[GenesFields.GENES_CHR_POS]), \
    #                                                                  float(sorting[sort_by]), sorting[GenesFields.GENE_ID_POS]))
    #    
    #    return sorted_list
    
    def load_annot(self, annot_path, annot_filename):
        self._annot_dict = {}
        
        ## pfam_src points out that not only PFAM can be the source of the data (Panther and other DBs as well)
        ## In fact, "annot_db_" would be a more general and better prefix for both "pfam_" fields
        
        annot_file = annot_path+"/"+annot_filename
        for annot_line in open(annot_file, 'r'):
            annot_line_data = annot_line.strip().split("\t")
            
            # PFAM
            if len(annot_line_data) > AnnotFile.PFAM_ID_POS:
                line_pfam_src = annot_line_data[AnnotFile.PFAM_SRC_POS]
                line_pfam_id = annot_line_data[AnnotFile.PFAM_ID_POS]
            else:
                continue # NO ANNOTATION, SO SKIP THIS ID
            
            # Gene ID
            gene_id = annot_line_data[AnnotFile.GENE_ID_POS]
            
            # Human readable description
            if len(annot_line_data) > AnnotFile.READABLE_POS: line_readable = annot_line_data[AnnotFile.READABLE_POS]
            else: line_readable = ''
            
            # InterPro
            if len(annot_line_data) > AnnotFile.IPR_POS: line_ipr = annot_line_data[AnnotFile.IPR_POS]
            else: line_ipr = ''
            
            # GOs
            if len(annot_line_data) > AnnotFile.GO_POS: line_go = annot_line_data[AnnotFile.GO_POS]
            else: line_go = ''
            
            if gene_id in self._annot_dict:
                gene_annot_data = self._annot_dict[gene_id]
                if line_ipr != '' and line_ipr not in gene_annot_data['ipr']:
                    gene_annot_data['ipr'].append(line_ipr)
                
                if line_readable != '' and line_readable not in gene_annot_data['readable']:
                    gene_annot_data['readable'].append(line_readable)
                
                if line_pfam_id != '' and line_pfam_id not in gene_annot_data['pfam_id']:
                    gene_annot_data['pfam_src'].append(line_pfam_src)
                    gene_annot_data['pfam_id'].append(line_pfam_id)
                
                if line_go != '' and line_go not in gene_annot_data['go']:
                    gene_annot_data['go'].append(line_go)
                
            else:
                gene_annot_data = {'readable':[], 'ipr':[], 'go':[], \
                                        'pfam_src':[], 'pfam_id':[]}
                self._annot_dict[gene_id] = gene_annot_data
                if line_ipr != '':
                    gene_annot_data['ipr'].append(line_ipr)
                    gene_annot_data['readable'].append(line_readable)
                
                if line_pfam_id != '':
                    gene_annot_data['pfam_src'].append(line_pfam_src)
                    gene_annot_data['pfam_id'].append(line_pfam_id)
                
                if line_go != '':
                    gene_annot_data['go'].append(line_go)
        
        if self._verbose: sys.stderr.write("GenesFacade: # annotated genes --> "+str(len(self._annot_dict))+"\n")
        
        return

    def append_annot(self, genetic_map, annot_dict):
        
        num_annotated = 0
        num_void = 0
        gene_list = self._genes_dict[genetic_map]
        
        for gene_data in gene_list:
            gene_id = gene_data[0]
            
            if gene_id in annot_dict:
                annot_data = annot_dict[gene_id]
                self._annotator.annotate_gene(gene_data, annot_data)
                num_annotated += 1
            else:
                self._annotator.annotate_void(gene_data)
                num_void += 1
        
        if self._verbose: sys.stderr.write("GenesFacade: Map: "+genetic_map+" --> # annotated genes "+str(num_annotated)+"\n")
        if self._verbose: sys.stderr.write("GenesFacade: Map: "+genetic_map+" --> # no annot genes "+str(num_void)+"\n")
        
        return

## END