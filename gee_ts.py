
# %%
import numpy as np
import geemap
import ee
ee.Initialize()


# %%

# %%
regions = (ee.FeatureCollection('FAO/GAUL/2015/level1')
             .filter(ee.Filter.eq('ADM0_NAME', 'United States of America'))
             .filter(ee.Filter.eq('ADM1_NAME', 'California')))
# us_states.getInfo()

ee_reducer_fn = ee.Reducer.sum()


ic_name = 'UCSB-CHG/CHIRPS/PENTAD'
years = np.arange(2005, 2006 + 1, 1)
years = [int(x) for x in years]
yr = years[0]
mo_start = 10

# %%
# Map = geemap.Map()

# Map.addLayer(california,{},'california')
# Map.centerObject(california, 6)
# Map

# %%

ee_years = ee.List.sequence(years[0], years[-1])
ic = ee.ImageCollection(ic_name)

def get_im_yr(yr):
    ee_yr_start = ee.Date.fromYMD(yr, mo_start, 1)
    ee_yr_end = ee_yr_start.advance(1, 'year')

    im_yr = (ic
            .filterDate(ee_yr_start, ee_yr_end)
            .reduce(ee_reducer_fn))
    yr_region = im_yr.reduceRegions(collection = regions, reducer = ee.Reducer.mean())

    return yr_region

yrs_out = ee_years.map(get_im_yr)

# %%

ee.FeatureCollection(ee_years).getInfo()


# %%

Map.addLayer(im_yr.clip(california),{'min':0, 'max':2000},'im_yr')
Map.addLayerControl()

# %%


yr_region.getInfo()
# %%

yr_region.getInfo()