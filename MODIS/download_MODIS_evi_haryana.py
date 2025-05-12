import ee
import ee.batch
import numpy as np
import pandas as pd
import os

print(os.getcwd())
print(os.listdir())
import rs

ee.Initialize(project = 'myresearch-421903')

# Load subdistricts of Punjab
# subdist = ee.FeatureCollection('projects/ee-gopalpenny/assets/INDIA_SUBDISTRICT_11')
# subdist_punjab = subdist \
#   .filter(ee.Filter.eq('STATE_UT','Punjab'))

haryana_fc = ee.FeatureCollection('projects/myresearch-421903/assets/haryana')
haryana_ft = ee.Feature(haryana_fc.union().first())
haryana_ft_buffer = haryana_ft.buffer(15000)
# haryana_ft_buffer.getInfo()
# .buffer(10) # buffer to ensure all pixels are included

# task3 = ee.batch.Export.table.toDrive(
#     collection = ee.FeatureCollection(haryana_ft_buffer),
#     description = 'haryana_ft_buffer4',
#     folder = 'GEE_DIV_MODIS_EVI_haryana',
#     fileFormat = 'SHP',
#     fileNamePrefix= 'haryana_ft_buffer4'
# )
# task3.start()

# haryana_ft_buffer.getInfo()
# punjab_subdist_all = subdist_punjab.aggregate_array('NAME').getInfo()
# len(punjab_subdist_all)
  # .filter(ee.Filter.eq('NAME','Abohar'))

# Set path to google drive
gee_dirname = 'GEE_DIV_MODIS_EVI_haryana'

region_dir_name = gee_dirname + '_evi'
gee_div_paths = ['/Users/gopal/Google Drive/_Research/Research projects/G_D/GEE_Div',
                 '/Users/gpenny/Google Drive/My Drive/_Research/Research projects/G_D/GEE_Div']
gee_div_path = [path for path in gee_div_paths if os.path.exists(path)][0]
gee_path_local = os.path.join(gee_div_path, gee_dirname)
region_path_local = os.path.join(gee_path_local, region_dir_name)
if not os.path.exists(region_path_local):
  os.mkdir(region_path_local)

# Load MODIS EVI data
def get_date_from_system_index(img):
  date = img.get('system:index')
  return img.set('date',date)
modis_terra_evi = ee.ImageCollection("MODIS/006/MOD13Q1").map(get_date_from_system_index)
modis_aqua_evi = ee.ImageCollection("MODIS/061/MYD13Q1").map(get_date_from_system_index)

# modis_terra_evi.first().getInfo()


modis_evi_merged = modis_terra_evi.merge(modis_aqua_evi)
modis_evi_merged = modis_evi_merged \
  .filter(ee.Filter.calendarRange(2000, 2023, 'year')) \
  .filter(ee.Filter.calendarRange(4, 12, 'month')) #.filterDate('2020-04-01','2020-04-10')

# convert SummaryQA band from bitmask 0-1 uint to 0-3 int
def prep_modis_bands(img):
    qa_summary_band = rs.extractQABits(img.select('SummaryQA'), 0, 1).int16()
    # qa_detailed_VIquality = rs.extractQABits(img.select('DetailedQA'), 0, 1).int16()
    # qa_detailed_VIusefulness = rs.extractQABits(img.select('DetailedQA'), 2, 5).int16()
    # qa_detailed_aerosol_quantity = rs.extractQABits(img.select('DetailedQA'), 6, 7).int16()
    # qa_detailed_adjacent_cloud_detected = rs.extractQABits(img.select('DetailedQA'), 8, 8).int16()
    # qa_detailed_atmospheric_correction = rs.extractQABits(img.select('DetailedQA'), 9, 9).int16()
    # qa_detailed_mixed_clouds = rs.extractQABits(img.select('DetailedQA'), 10, 10).int16()
    # qa_detailed_land_water_mask = rs.extractQABits(img.select('DetailedQA'), 11, 13).int16()
    # qa_detailed_possible_snow_ice = rs.extractQABits(img.select('DetailedQA'), 14, 14).int16()
    # qa_detailed_possible_shadow = rs.extractQABits(img.select('DetailedQA'), 15, 15).int16()
    # img_out = img.select(['EVI']).int16() \
    #     .addBands(qa_summary_band.rename('SummaryQA')) \
    #     .addBands(qa_detailed_VIquality.rename('DetailedQA_VIquality')) \
    #     .addBands(qa_detailed_VIusefulness.rename('DetailedQA_VIusefulness')) \
    #     .addBands(qa_detailed_aerosol_quantity.rename('DetailedQA_aerosol_quantity')) \
    #     .addBands(qa_detailed_adjacent_cloud_detected.rename('DetailedQA_adjacent_cloud_detected')) \
    #     .addBands(qa_detailed_atmospheric_correction.rename('DetailedQA_atmospheric_correction')) \
    #     .addBands(qa_detailed_mixed_clouds.rename('DetailedQA_mixed_clouds')) \
    #     .addBands(qa_detailed_land_water_mask.rename('DetailedQA_land_water_mask')) \
    #     .addBands(qa_detailed_possible_snow_ice.rename('DetailedQA_possible_snow_ice')) \
    #     .addBands(qa_detailed_possible_shadow.rename('DetailedQA_possible_shadow'))
    img_out = img.select('EVI').updateMask(qa_summary_band.lte(1)) # return only good quality pixels
    return img_out

modis_evi = ee.ImageCollection(modis_evi_merged.map(prep_modis_bands))

# modis_evi.first().getInfo()
# modis_evi.size().getInfo()

# Get list of images of modis_evi
modis_images = modis_evi.aggregate_array('date').getInfo()
modis_images.sort()
# modis_evi.filter(ee.Filter.eq('date',modis_images[32])).first().projection().getInfo()
# modis_evi.filter(ee.Filter.eq('date',modis_images[33])).first().projection().getInfo()

# modis_images = modis_images[0:1]
# modis_evi.first().getInfo()

# Export each image in modis_images to google drive
proj_region = modis_terra_evi.first().projection() # get single projection to set as primary. 
for image in modis_images:
  img_name = 'modis_evi_haryana_' + image
  if not os.path.exists(os.path.join(region_path_local, img_name + '.tif')):
    task = ee.batch.Export.image.toDrive(
        image = modis_evi.filter(ee.Filter.eq('date',image)).first().reproject(proj_region), # ensure single projection
        description = image,
        folder = img_name,
        fileNamePrefix = img_name,
        region = haryana_ft_buffer.geometry(),
        scale = 250,
        maxPixels = 1e13
    )
    task.start()
    print(task.status())
  else:
    print('Image' + img_name + 'already exists')



# task = ee.batch.Export.image.toDrive(


# haryana_ft_buffer.getInfo()
