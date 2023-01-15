#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 17:21:18 2023

@author: gopal
"""


import numpy as np
import pandas as pd
import geopandas as gpd
import geemap
import ee
import re
# import plotnine as p9
# import sys
# from datetime import timedelta
import os
from datetime import datetime

def download_gee_rasters(ic_name, shp_name, shp_dir, google_drive_path, 
                         start_date, end_date, doy_increment):

    ee.Initialize()
    
    boundary_gpd = gpd.read_file(os.path.join(shp_dir, shp_name + ".shp"))
    boundary_ee = geemap.geopandas_to_ee(boundary_gpd)
    
    assert boundary_ee.size().getInfo() == 1 # must only have 1 feature in the shapefile
    
    # convert from feature collection to feature to geometry
    ee_geometry = boundary_ee.first().geometry()
    
    shp_gdrive_path = os.path.join(google_drive_path,shp_name)
    shp_gee_folder_name = ic_name + '_' + shp_name + '_GEE'
    shp_gee_folder_path = os.path.join(shp_gdrive_path, shp_gee_folder_name)
    
    
    if not os.path.exists(google_drive_path):
        os.mkdir(google_drive_path)
    if not os.path.exists(shp_gdrive_path):
        os.mkdir(shp_gdrive_path)
    if not os.path.exists(shp_gee_folder_path):
        os.mkdir(shp_gee_folder_path)
    
    # Map = geemap.Map()
    # Map.addLayer(boundary_ee,{},'boundary_ee')
    # Map.centerObject(boundary_ee, 13)
    
    if ic_name == "s1":
        ic_ee_identifier = "COPERNICUS/S1_GRD"
        ic_table_path = os.path.join(shp_gdrive_path, 's1_table.csv')
    
    else:
        raise Exception("Sorry, ic_name can only be s1 at the moment")
    
    ic = ee.ImageCollection(ic_ee_identifier) \
        .filterBounds(ee_geometry) # \
        # .filterDate(start_date, end_date) # not necessary here -- use only for downloading
    
    if not os.path.exists(ic_table_path):
        ic_table = create_ic_table(ic, ic_table_path, doy_increment)
    else:
        ic_table = pd.read_csv(ic_table_path)
    
    idx_download = ((ic_table['date'] >= start_date) & 
                    (ic_table['date'] <= end_date) & 
                    (ic_table['sample']))
    ic_names_to_download = ic_table.names[idx_download].values
    
    # %%
    for im_name in ic_names_to_download:
        update_ic_table_all(ic_table_path, shp_gee_folder_path)
        download_shp_ic(ic, shp_name, im_name, ee_geometry, shp_gee_folder_name, ic_table_path)
        
    return None


def create_ic_table(ic, ic_table_path, doy_increment):
    
    print("Creating ic_table...")
    
    doy_select = np.arange(1, 370, doy_increment)

    
    print("Aggregating all image names...")
    ic_image_names = ic.aggregate_array("system:index").getInfo()
    
    ic_table = pd.DataFrame({'names' : ic_image_names})
    ic_table['datestr'] = [re.sub(".*_.*_.*_.*_([0-9]+)T[0-9]+_.*","\\1",x) for x in ic_table["names"]]
    ic_table['date'] = pd.to_datetime(ic_table['datestr'], format = "%Y%m%d")
    ic_table = ic_table.sort_values('date')
    
    # ic_table['day_diff'] = 1
    # print(np.array(ic_table.loc[1:2,'date']))
    # print(np.array(ic_table.loc[0:1,'date']))
    # print(ic_table['date'].shift(0))
    ic_table['days_until_next'] = (ic_table['date'].shift(-1) - ic_table['date']).dt.days
    
    # np.array(ic_table.loc[1:2,'date']) - np.array(ic_table.loc[0:1,'date'])
    ic_table['doy'] = [int(datetime.strftime(x, '%j')) for x in ic_table['date']]
    ic_table['year'] = [int(datetime.strftime(x, '%Y')) for x in ic_table['date']]
    
    # sample ic_table names to those closest to doy_select (i.e., doy increment)
    # get names in each year
    
    print(f"Resampling images dates based on {doy_increment}-day increment...")
    resample_names = []
    for year in pd.unique(ic_table['year']):
        # year = 2016
        days_orig = np.array(ic_table[ic_table['year'] == year]['doy'])
        days_resampled = resample_nearest_days(doy_select, days_orig)
        # days_resampled.shape
        ic_im_year = ic_table[ic_table['year'] == year]
        # ic_im_year = ic_im_year[~ic_im_year['doy'].duplicated()]
        ic_im_year_resampled_names = ic_im_year.names[
            (ic_table['doy'].isin(days_resampled)) &
            (~ic_im_year['doy'].duplicated())]
        resample_names = resample_names + ic_im_year_resampled_names.values.tolist()
    
    
    ic_table['sample'] = ic_table['names'].isin(resample_names)
    
    # Print number of sampled images
    print("Total number of images:",ic_table.shape[0])
    print("Number of selected images:",np.sum(ic_table['sample'].values))
    
    ic_table['downloaded'] = 0
    ic_table['gee_task_sent'] = ""
    
    # sae the image table
    print("Saving ic_table...")
    ic_table.to_csv(ic_table_path, index = False)

    return ic_table


def resample_nearest_days(doy_select, days_orig):
    """
    This function returns an np array, days_nearest, which reflects the unique
    values of days_orig that are closest to each value of doy_select
    doy_select: np.array
    days_orig: np.array
    """
    doy_select_mat = np.expand_dims(days_orig, 0).repeat(len(doy_select), 0)
    days_orig_mat = np.expand_dims(doy_select, 1).repeat(len(days_orig), 1)
    doy_nearest_idx = np.argmin(np.abs(days_orig_mat - doy_select_mat),axis = 1)
    days_nearest = days_orig[np.unique(doy_nearest_idx)]
    return days_nearest


def check_ic_table_image(ic_table_path, im_name):
    # check if we can download file

    # 1 read ic_table
    ic_table = pd.read_csv(ic_table_path, keep_default_na = False)

    # 2 subset to image
    ic_table_im = ic_table[ic_table['names'] == im_name]
    downloaded = ic_table_im.downloaded.values.item()
    gee_task_sent = ic_table_im.gee_task_sent.values.item()
    # 3 return status as "gee_task_sent", "downloaded", or "not_downloaded"
    if downloaded:
        im_status = 'downloaded'
    elif not gee_task_sent == '':
        im_status = 'gee_task_sent'
    else:
        im_status = 'not_downloaded'

    return im_status

def update_ic_table_gee_task_sent(ic_table_path, im_name, val):
    # update specific image in ic_table with val
    # 1 read ic_table
    ic_table = pd.read_csv(ic_table_path)

    # 2 insert val to talbe
    ic_table.loc[ic_table['names'] == im_name, 'gee_task_sent'] = val
    
    ic_table.to_csv(ic_table_path, index = False)

    return ic_table


def update_ic_table_all(ic_table_path, shp_gee_folder_path):
    # check all downloaded files

    # 1 read ic_table
    ic_table = pd.read_csv(ic_table_path, keep_default_na = False)
    
    # 2 check/updated all downloaded files
    images_downloaded = [re.sub("\\.tif$","",x) for x in os.listdir(shp_gee_folder_path)]

    ic_table['downloaded'] = ic_table['names'].isin(images_downloaded)
    ic_table.to_csv(ic_table_path, index = False)
    
    return ic_table


def download_shp_ic(ic, shp_name, im_name, ee_geometry, shp_gee_folder_name, ic_table_path):
    
    im_status = check_ic_table_image(ic_table_path, im_name)
    
    if im_status == "not_downloaded":
        ic_im_ic = ic.filterMetadata('system:index','equals',im_name)
        assert ic_im_ic.size().getInfo() == 1
        ic_im = ic_im_ic.first()
        ic_im_shp_name = im_name + '_' + shp_name

        print(f'Sending {shp_name}/{im_name} to GEE')
        task = ee.batch.Export.image.toDrive(
            image = ic_im.float(),
            description = ic_im_shp_name,
            folder = shp_gee_folder_name,
            fileNamePrefix = im_name,
            region = ee_geometry,
            scale = 10,
            maxPixels = 1000000,
            fileFormat = 'GeoTIFF')
        task.start()
        
        val = datetime.strftime(datetime.today(), "%Y-%m-%d")
        update_ic_table_gee_task_sent(ic_table_path, im_name, val)
        
    else:
        print(f'Image {shp_name}/{im_name} skipped, due to status: {im_status}')
    
    return 1