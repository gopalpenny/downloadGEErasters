import ee
import numpy as np
import pandas as pd
import os
import rs

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

chirps_precip = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD") \
  .filter(ee.Filter.calendarRange(2016, 2022, 'year'))
# years.get(0).add(1).getInfo()
# chirps_precip.first().getInfo()

years = np.arange(2016,2022)

  
# monsoon_year = years[0]
def get_annual_et(monsoon_year):
    start_date = ee.Date.fromYMD(int(monsoon_year),6,1)
    end_date = ee.Date.fromYMD(int(monsoon_year)+1,5,31)
    chirps_annual_precip = chirps_precip \
        .filterDate(start_date, end_date) \
        .sum() \
        .set('system:time_start', start_date.millis()) \
        .set('system:time_end', end_date.millis()) \
        .set('system:index', str(monsoon_year))

    task = ee.batch.Export.image.toDrive(
        image = chirps_annual_precip.select('precipitation'),
        description = 'MODIS_precip_' + str(monsoon_year),
        folder = gee_dirname,
        region = watershed.geometry().bounds(),
        scale = 5566,
        maxPixels = 1e13
    )
    task.start()
    return(task)

# get_annual_et(2016)

for i, monsoon_year in enumerate(years):
    if not os.path.exists(os.path.join(gee_path_local, 'MODIS_precip_' + str(monsoon_year) + '.tif')):
        task = get_annual_et(monsoon_year)
        print('Exporting precipitation for monsoon year: ', monsoon_year)
    else:
        print('Precip for monsoon year: ', monsoon_year, ' already exists')

# task.status()
