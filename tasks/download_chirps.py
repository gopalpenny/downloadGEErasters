#### NOT WORKING FOR SOME REASON ####
# %%
import numpy as np
import ee
ee.Initialize()

# %%

start_year = 2013
end_year = 2022

years = np.arange(start_year, end_year + 1, 1)

# %%
year = years[0]
start_date = str(years[0]) + '-01-01'
end_date = str(years[-1]) + '-12-31'
num_years = years[-1] - years[0] + 1

# %%

countries = ee.FeatureCollection("FAO/GAUL/2015/level0") \
    .filter(ee.Filter.inList('ADM0_NAME',['India','Bangladesh','Nepal']))

rect = ee.Geometry.Rectangle([68, 7, 97.5, 34]) #.getInfo()

# %%

chirps_annual_avg = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD') \
    .filterDate(start_date, end_date) \
    .select('precipitation') \
    .reduce(ee.Reducer.sum()) 

# %%
chirps_annual_avg.getInfo()

# %%
filename = f'chirps_precip_sum_{start_year}_{end_year}'

# %%
# NOT WORKING FOR SOME REASON
task_chirps = ee.batch.Export.image.toDrive(
    image= chirps_annual_avg,
    description= filename,
    fileNamePrefix= filename,
    folder= 'GEE_Div',
    dimensions= 1000,
    region= rect
)

# task_chirps.start()

# %% GFSAD
gfsad1k = ee.Image("USGS/GFSAD1000_V1")

# %%
gfsad1k.getInfo()

# %%
# NOT WORKING FOR SOME REASON
task_gfsad1k = ee.batch.Export.image.toDrive(
    image = gfsad1k,
    description = 'gfsad1k_Sasia',
    fileNamePrefix = 'gfsad1k_Sasia',
    folder = 'GEE_Div',
    region = rect
)

task_gfsad1k.start()

# %%
ee.batch.Task.list()[0]
# %%
