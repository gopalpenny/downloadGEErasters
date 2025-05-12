import ee
import numpy as np
import pandas as pd
import os

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

srtm = ee.Image(ee.Image("CGIAR/SRTM90_V4"))
# years.get(0).add(1).getInfo()
# chirps_precip.first().getInfo()

task = ee.batch.Export.image.toDrive(
        image = srtm.select('elevation'),
        description = 'SRTM_elevation',
        folder = gee_dirname,
        region = watershed.geometry().bounds(),
        scale = 90,
        maxPixels = 1e13
    )
task.start()