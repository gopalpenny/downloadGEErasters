#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 17:21:39 2023

@author: gopal
"""


from gee_raster_download_functions import download_gee_rasters

shp_name = "SquareBuffer1km" # name of shapefile (without .shp)
shp_dir = "./UbonRatchathani" # local path to shapefile 

# path to google drive folder
google_drive_path = '/Users/gopal/Google Drive/_Research/Research projects/ML/downloadsentinel'
start_date = '2016-01-01' # start date for imagery download
end_date = '2017-01-01' # end date for imagery download
ic_name = 's1' # currently can only be 's1' (can be "s2" in future)
doy_increment = 6 # days which to sample image collection (days)

download_gee_rasters(ic_name, shp_name, shp_dir, google_drive_path, 
                         start_date, end_date, doy_increment)