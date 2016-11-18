#!/usr/bin/env python
# -*- coding: utf-8 -*-

# data_utils.py is part of Barleymap.
# Copyright (C)  2013-2014  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

import sys

from barleymapcore.m2p_exception import m2pException

def read_paths(config_file_path): # TODO pass this to utils package
    config_path_dict = {}
    
    sys.stderr.write("Reading paths from config file...\n")
    
    for config_line in open(config_file_path, 'r'):
        if config_line.startswith("#"): continue
        config_data = config_line.strip().split(" ")
        config_path_dict[config_data[0]] = config_data[1]
    
    sys.stderr.write("Config file read.\n")
    
    return config_path_dict

## NEXT CODE IS FIDDLY. WOULD BE BETTER WITH A SORTED DICT
# Can be used to load databases, datasets, everything that is
# configured in a conf_file with key-->value
# and the user choice is a list of "," or default all
def load_data(conf_file, users_list = None, verbose = False):
    databases_list = []
    databases_param = []
    
    if verbose: sys.stderr.write("Configured data:\n")
    
    databases_param_default = []
    databases_configured = {}
    for database in open(conf_file, 'r'):
        database_data = database.strip().split(" ")
        databases_param_default.append(database_data[0])
        databases_configured[database_data[0]] = database_data[1]
        
        if verbose: sys.stderr.write("\t"+database_data[0]+" --> ID: "+database_data[1]+"\n")
    
    if users_list:
        databases = users_list.split(",")
        for database in databases:
            if database in databases_configured:
                databases_list.append(databases_configured[database])
                databases_param.append(database)
            else:
                raise m2pException("Warning: data "+database+" is not present in config file.")
        
        databases_param = ",".join(databases_param)
    else:
        databases_list = [database for database in [databases_configured[db] for db in databases_param_default]]
        databases_param = ",".join(databases_param_default)
        
    return (databases_param, databases_list)

def load_datasets(conf_file, verbose = False):
    datasets_dict = {}
    
    if verbose: sys.stderr.write("Loading datasets configuration...\n")
    
    for line in open(conf_file, 'r'):
        dataset_data = line.strip().split(" ")
        if verbose: sys.stderr.write("\t dataset: "+str(dataset_data)+"\n")
        dataset_name = dataset_data[0]
        dataset_id = dataset_data[1]
        dataset_type = dataset_data[2]
        
        if dataset_id in datasets_dict:
            raise Exception("Duplicated dataset "+str(dataset_id)+" in config file.")
        else:
            datasets_dict[dataset_id] = {"dataset_id":dataset_id, "dataset_name":dataset_name, "dataset_type":dataset_type}
    
    if verbose: sys.stderr.write("Datasets configuration loaded.\n")
    
    return datasets_dict

## END