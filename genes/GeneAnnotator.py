#!/usr/bin/env python
# -*- coding: utf-8 -*-

# GeneAnnotator.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

from GenesBase import AnnotFields

class GeneAnnotator(object):
    
    def annotate_gene(self, gene_data, annot_data):
        gene_data.append(";".join(annot_data['readable']))
        gene_data.append(",".join(annot_data['ipr']))
        gene_data.append(",".join(annot_data['pfam_id']))
        gene_data.append(",".join(annot_data['pfam_src']))
        gene_data.append(",".join(annot_data['go']))
        
        return
    
    def annotate_void(self, gene_data):
        gene_data.extend(AnnotFields.GENES_ANNOT_FIELDS*["-"])
        return

## END
