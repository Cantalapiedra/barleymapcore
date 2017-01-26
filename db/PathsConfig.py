#!/usr/bin/env python
# -*- coding: utf-8 -*-

# PathsConfig.py is part of Barleymap.
# Copyright (C)  2016  Carlos P Cantalapiedra.
# (terms of use can be found within the distributed LICENSE file).

from barleymapcore.utils.data_utils import read_paths
from barleymapcore.db.ConfigBase import ConfigBase

class PathsConfig(object):
    
    _config_path_dict = None
    
    # Keys in config file
    _APP_PATH = "app_path"
    _DATASETS_PATH = "datasets_path"
    
    # Values read from config file
    _app_path = ""
    _datasets_path = ""
    
    def __init__(self, app_abs_path):
        self.load_config(app_abs_path)
    
    def load_config(self, app_abs_path):
        paths_conf_file = app_abs_path+ConfigBase.PATH_FILE
        self._config_path_dict = read_paths(paths_conf_file)
        
        self._app_path = self._config_path_dict[self._APP_PATH]
        self._datasets_path = self._config_path_dict[self._DATASETS_PATH]

    def get_app_path(self):
        return self._app_path
    
    def get_datasets_path(self, ):
        return self._datasets_path
    