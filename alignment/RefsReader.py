#!/usr/bin/env python
# -*- coding: utf-8 -*-

# RefsReader.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

class RefsReader(object):
    
    # Fields in references.conf file
    REF_CONFIG_REF_NAME = 0
    REF_CONFIG_REF_ID = 1
    REF_CONFIG_REF_TYPE = 2
    
    _config_file = ""
    _verbose = False
    
    _config_dict = {}
    # [ref_id] = {"ref_name", "ref_type"}
    
    def __init__(self, config_file, verbose = True):
        self._config_file = config_file
        self._load_config()
        self._verbose = verbose
        
    
    def _load_config(self):
        for config_line in open(self._config_file, 'r'):
            config_data = config_line.strip().split(" ")
            ref_id = config_data[self.REF_CONFIG_REF_ID]
            ref_name = config_data[self.REF_CONFIG_REF_NAME]
            ref_type = config_data[self.REF_CONFIG_REF_TYPE]
            # ref_type:
            # ----big (AlignmentFacade._REF_TYPE_BIG)
            # ----std (AlignmentFacade._REF_TYPE_NORMAL)
            
            self._config_dict[ref_id] = {"ref_name":ref_name, "ref_type":ref_type}
        
    
    # Obtain type (big or std size) of fasta DB
    def get_ref_type(self, db):
        return self._config_dict[db]["ref_type"]

## END