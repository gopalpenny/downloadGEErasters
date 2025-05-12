import ee
import numpy as np
import pandas as pd
import os
# import rs

ee.Initialize()

# Load subdistricts of Punjab
watershed = ee.FeatureCollection('projects/ee-gopalpenny/assets/india/cauvery_from_raster_250m_buffer_simplified')
# watershed.getInfo()

# punjab_subdist_all = subdist_punjab.aggregate_array('NAME').getInfo()
# len(punjab_subdist_all)
  # .filter(ee.Filter.eq('NAME','Abohar'))

# Set path to google drive
gee_dirname = 'GEE_cauvery_download'

# region_dir_name = gee_dirname + '_evi'
gee_path_local = os.path.join('/Users/gopal/Google Drive/_Research/Research projects/ML/download_gee_rasters/cauvery', gee_dirname)
# os.mkdir(gee_path_local)
# region_path_local = os.path.join(gee_path_local, region_dir_name)
# if not os.path.exists(region_path_local):
#   os.mkdir(region_path_local)

# Load MODIS EVI data
def get_date_from_system_index(img):
  date = img.get('system:index')
  return img.set('date',date)


modis_8day_et = ee.ImageCollection("MODIS/061/MOD16A2GF") \
  .filter(ee.Filter.calendarRange(2016, 2022, 'year'))
# years.get(0).add(1).getInfo()


years = np.arange(2016,2022)

  
# monsoon_year = years[0]
def get_annual_et(monsoon_year):
    start_date = ee.Date.fromYMD(int(monsoon_year),6,1)
    end_date = ee.Date.fromYMD(int(monsoon_year)+1,5,31)
    modis_annual_et = modis_8day_et \
        .filterDate(start_date, end_date) \
        .sum() \
        .set('system:time_start', start_date.millis()) \
        .set('system:time_end', end_date.millis()) \
        .set('system:index', str(monsoon_year))

    task = ee.batch.Export.image.toDrive(
        image = modis_annual_et.select(['ET','PET']).float(),
        description = 'MODIS_PET_' + str(monsoon_year),
        folder = gee_dirname,
        region = watershed.geometry().bounds(),
        scale = 250,
        crs = 'EPSG:4326',
        maxPixels = 1e13
    )
    task.start()
    return(task)


for i, monsoon_year in enumerate(years):
    if not os.path.exists(os.path.join(gee_path_local, 'MODIS_PET_' + str(monsoon_year) + '.tif')):
        task = get_annual_et(monsoon_year)
        print('Exporting MODIS_PET_ for monsoon year: ', monsoon_year)
    else:
        print('MODIS_PET_ for monsoon year: ', monsoon_year, ' already exists')

# task.status()
