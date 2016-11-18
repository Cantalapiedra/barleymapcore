#!/usr/bin/env python
# -*- coding: utf-8 -*-

# GenesBase.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

class GenesFields(object):
    
    GENES_ID_POS = 0
    GENES_TYPE_POS = 1
    GENES_MAP_POS = 2
    GENES_CHR_POS = 3
    GENES_CM_POS = 4
    GENES_BP_POS = 5
    
    GENES_FIELDS = 6
    
class AnnotFields(object):
    
    GENES_ANNOT_DESC = 0
    GENES_ANNOT_INTERPRO = 1
    GENES_ANNOT_PFAM = 2
    GENES_ANNOT_SERVER = 3
    GENES_ANNOT_GO = 4
    
    GENES_ANNOT_FIELDS = 5
    
class AnnotFile(object):

    GENE_ID_POS = 0
    PFAM_SRC_POS = 3
    PFAM_ID_POS = 4
    READABLE_POS = 12
    GO_POS = 13
    IPR_POS = 11

## END