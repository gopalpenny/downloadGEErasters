#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 17:21:39 2023

@author: gopal
"""

from gee_raster_download_functions import download_gee_rasters

# OVERVIEW
# This script downloads raster data within a polygon specified by shapefile_name.
# The function identifies all rasters within the ee Image Collection (ic_name) that
# overlap the shapefile. Because there may be an unnecessarily large number
# of rasters, the set of rasters is re-sampled to select only those that
# are closest to a regular interval of days each year. For instance, if 
# doy_increment = 5, the functions selects the raster closest to doy = 1 (Jan 1), 
# the raster closest to doy = 6 (Jan 6), doy = 11 (Jan 11), and so on. 

# GUIDE
# Create a directory in the Google drive on your computer to manage data
# and download the raster. In the example below this directory is named 
# 'download_gee_rasters'. Set the other configuration parameters accordingly.
# Then run this script by loading a conda environment and typing:
    
# python gee_raster_config_and_download.py

# The dates can be set to a narrow range at first, for testing. Then the range
# can be expanded. Images will not be re-downloaded. 

# The first time a configuration file is run with a particular shapefile,
# a new directory is created as a subdirectory of [google_drive_path].
# This directory contains a csv file to keep track of the status of all images,
# [ic-name]_table.csv. This file keeps track of which rasters should be downloaded
# (the "sample" column), which rasters have been sent to Google Earth Engine
# (the "gee_task_sent" column) and which rasters have downloaded ("downloaded").
# The directory will also have a subdirectory where GEE will save the downloaded 
# rasters.

# The rasters are saved in float32 format (opposed to float64) to limit
# harddrive space.

# CONFIGURATION PARAMETERS

path_params = {
    'shapefile_name' : "SquareBuffer1km", # name of shapefile (without .shp) (contains a single polygon)
    'shapefile_directory' : "./UbonRatchathani", # local path to shapefile directory
    # path to google drive folder
    'google_drive_path' : '/Users/gopal/Google Drive/_Research/Research projects/ML/download_gee_rasters'
    }

date_params = {
    'start_date' : '2019-05-15', # start date for imagery download
    'end_date' : '2019-06-01', # end date for imagery download
    'doy_increment' : 6 # days which to sample image collection (days)
    }

ic_s1_params = {
    'ic_name' : 's1' # currently can only be 's1' (can be "s2" in future)
    }

ic_s2_params = {
    'ic_name' : 's2', # currently can only be 's1' (can be "s2" in future)
    'tile_option' : 0, # index of tile, ordered alphabetically
    # 'tile_option' : '55PWB', # Or, can be name of tile or "any"
    'max_cloud_pct' : 50 # exclude tiles with > 50% clouds
    }

ic_oli8_params = {
    'ic_name' : 'oli8', # currently can only be 's1' (can be "s2" in future)
    'rowpath_option' : 0, # can be index of row/path (0, 1, etc, must not exceed number of tiles)
    # 'rowpath_option' : '55PWB' # Or, can be name of row/path or "any"
    'max_cloud_pct' : 60 # exclude tiles with > 60% clouds
    }


# download_gee_rasters(ic_s1_params, path_params, date_params)
# download_gee_rasters(ic_s2_params, path_params, date_params)
download_gee_rasters(ic_oli8_params, path_params, date_params)