"""
RS functions
1. cloudScore:
    --Cloud scoring function
2. fill_l7_slc:    
    --Fill LANDSAT 7 SLC-off pixels with successive dilation
3. get_raster_pts: 
    --Get feature collection of gridded points of raster pixels within a region
4. get_pixel_timeseries:
    --Get timeseries at all points in a feature collection of all images in an image collection
5. extractQABits:
    --Extract QA bits using the position of the bit
6. Get monthly average raster values for a region
7. Get approximate monthly water balance for a region using CHIRPS and MODIS ET
8. Get a timeseries plot of monthly water balance for a region using CHIRPS and MODIS ET


# CALL FUNCTIONS:
# rs = require('users/gopalpenny/default:functions/RS_functions.js')
"""

import ee


def cloudScore(img):
  """1. Cloud scoring function
  Compute a cloud score.  This expects the input image to have the common
  band names: ["red", "blue", etc], so it can work across sensors.
  A helper to apply an expression and linearly rescale the output.
  """
  def rescale(img, exp, thresholds):
    return img.expression(exp, {img: img}) \
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

  # Compute several indicators of cloudyness and take the minimum of them.
  score = ee.Image(1.0)
  # Clouds are reasonably bright in the blue band.
  score = score.min(rescale(img, 'img.blue', [0.1, 0.3]))

  # Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]))

  # Clouds are reasonably bright in all infrared bands.
  score = score.min(
      rescale(img, 'img.nir + img.swir1 + img.swir2', [0.3, 0.8]))

  # Clouds are reasonably cool in temperature.
  score = score.min(rescale(img, 'img.temp', [300, 290]))

  # However, clouds are not snow.
  ndsi = img.normalizedDifference(['green', 'swir1'])
  return score.min(rescale(ndsi, 'img', [0.8, 0.6]))


def dilate_l7_slc(im):
  """2. FILL LANDSAT 7 SLC-off pixels with successive dilation
  """ 
  im_dilate = im.focal_median({'kernel' : ee.Kernel.plus(30,"meters"),'iterations' : 1})
  im_dilate_only = im_dilate.updateMask(im.mask().add(-1).multiply(-1).multiply(im_dilate.mask())) # .add(-1).multiply(-1) is a "not", .multiply(x) is "and"
  im_update = im.unmask().add(im_dilate_only.unmask()).updateMask(im_dilate.mask()) \
      .set('system:time_start',im.get('system:time_start'))
  return im_update

def fill_l7_slc(im):
  im_update = dilate_l7_slc(im)
  im_update2 = dilate_l7_slc(im_update)
  im_update3 = dilate_l7_slc(im_update2)
  im_update4 = dilate_l7_slc(im_update3)
  im_update5 = dilate_l7_slc(im_update4)
  im_update6 = dilate_l7_slc(im_update5)
  im_update7 = dilate_l7_slc(im_update6)
  im_update8 = dilate_l7_slc(im_update7)
  return im_update8


def get_raster_pts(raster, region, scale = 0):
  """
  # 3. Get feature collection of gridded points of raster pixels within a region
  # call as: 
  # pts = rs.get_raster_pts(raster, region, scale)
  # get a geometry containing coordinates for for all gridded points within region of a raster
  # masked pixels in the raster are excluded from the output
  # scale parameter is optional -- by default, the scale is the nominalScale of the image.
  # The nominalScale can be overridden by providing a scale that is greater.
  """

  raster_clip = raster.select([0]).clip(region)
  # get image projection
  # scale = scale | 0 # not needed python
  proj = raster_clip.projection()
  proj_scale = ee.Number(proj.nominalScale()).max(scale)
  
  # get coordinates image
  latlon = ee.Image.pixelLonLat().reproject(proj)
  
  # put each lon lat in a list
  coords = latlon.select(['longitude', 'latitude']) \
    .updateMask(raster_clip.mask()) \
                   .reduceRegion(
    reducer = ee.Reducer.toList(),
    geometry = region,
    scale = proj_scale
  )
  
  # get lat & lon
  lat = ee.List(coords.get('latitude'))
  lon = ee.List(coords.get('longitude'))
  
  # zip them. Example: zip([1, 3],[2, 4]) --> [[1, 2], [3,4]]
  points_list = lon.zip(lat)
  
  def genFeature(latlon):
    return ee.Feature(ee.Geometry.Point(latlon))
  def setPt(pt):
    return pt.set('pt_id',ee.Number.parse(pt.id()).int()) 

  points_features = points_list.map(genFeature)
  points_fc = ee.FeatureCollection(points_features).map(setPt)
  
  return points_fc


def get_pixel_timeseries(pts_fc, image_collection, bands, ic_property_id, scale = 0):
  """
  Get timeseries at all points in a feature collection of all images in an image collection
  call as: 
  
  time_series = rs.get_pixel_timeseries(
    pts_fc = pts,
    image_collection = image_collection,
    bands = ['precipitation'],
    ic_property_id = 'system:index',
    scale = 30) # for Landsat resolution
  
  
  This function returns a feature collection of points, sampling an image collection 
  at each point in each image. 
  
  Parameters
  ----------
  pts_fc : ee.FeatureCollection
  - feature collection of points to sample
  image_collection : ee.ImageCollection
  - image collection to sample frome
  bands : list of strings
  - index (numbers) or names of bands, as a list -- [0] OR [NDVI','EVI'] # NOTE: if only 1 band, 'first' and band name are included in output (as duplicates)
  ic_property_id : str
  - name of the field containing image ID (probably a date) -- e.g., 'system:index'
  scale : int
  - scale at which to reduce pixels. defaults to nominalScale
  """
  ic = image_collection \
    .filterBounds(pts_fc) \
    .select(bands)
  proj = image_collection.first().select([0]).projection()
  proj_scale = ee.Number(proj.nominalScale()).max(scale)
  
  # # prep variables to change 'first' to band name if only one band
  bands = ee.List(bands)
  only_one_band = bands.size().eq(1) # 1 if one band, 0 otherwise
  name_first_result = ee.List([bands.get(0),'first']).get(only_one_band.int())
  
  def get_ts(img):
    pts_fc_im = pts_fc.filterBounds(img.geometry())
    feature_ts =  img.reduceRegions(pts_fc_im,ee.Reducer.first(),proj_scale) ##.rename('bands')

    def get_ts_image(feat):
      feat2 = feat \
        .set(bands.get(0),feat.get(name_first_result)) \
        .set('image_id',img.get(ic_property_id))
      return feat2

    ts_image = feature_ts.map(get_ts_image)
    return ts_image
  ts = ic.map(get_ts).flatten()
  
  # old_names = ts.first().propertyNames()
  # new_names = old_names.replace('first',ee.List(bands).get(0))

  return ts #.select(orig_names,new_names)

def get_pixel_ts_allbands(pts_fc, image_collection, ic_property_id, scale = 0):
  """
  Get timeseries at all points in a feature collection of all images in an image collection
  call as: 
  
  time_series = rs.get_pixel_ts_allbands(
    pts_fc = pts,
    image_collection = image_collection,
    ic_property_id = 'system:index',
    scale = 30) # for Landsat resolution
  
  
  This function returns a feature collection of points, sampling an image collection 
  at each point in each image. 
  
  Parameters
  ----------
  pts_fc : ee.FeatureCollection
  - feature collection of points to sample
  image_collection : ee.ImageCollection
  - image collection to sample frome
  ic_property_id : str
  - name of the field containing image ID (probably a date) -- e.g., 'system:index'
  scale : int
  - scale at which to reduce pixels. defaults to nominalScale
  """
  ic = image_collection \
    .filterBounds(pts_fc)
  proj = image_collection.first().select([0]).projection()
  proj_scale = ee.Number(proj.nominalScale()).max(scale)
  
  def get_ts(img):
    pts_fc_im = pts_fc.filterBounds(img.geometry())
    feature_ts =  img.reduceRegions(pts_fc_im,ee.Reducer.first(),proj_scale) ##.rename('bands')

    def get_ts_image(feat):
      feat2 = feat \
        .set('image_id',img.get(ic_property_id))
      return feat2

    ts_image = feature_ts.map(get_ts_image)
    return ts_image
  ts = ic.map(get_ts).flatten()
  
  # old_names = ts.first().propertyNames()
  # new_names = old_names.replace('first',ee.List(bands).get(0))

  return ts #.select(orig_names,new_names)

def extractQABits(qaBand, bitStart, bitEnd, RADIX = 2):
  """
  5. Extract QA bits using the position of the bit
  qaBand: image with a single QA band
  bitStart: index of the start of the QA bits, zero-referenced, from the right
  bitEnd:   index of the end   of the QA bits, zero-referenced, from the right

  rs.extractQABits(img.select('pixel_qa'),3,3).rename('shadow')
  rs.extractQABits(img.select('pixel_qa'),5,5).rename('clouds')
  """
  # RADIX = RADIX || 2
  numBits = bitEnd - bitStart + 1
  power = RADIX ** numBits
  qaBits = qaBand.rightShift(bitStart).mod(power)
  return qaBits

def get_qaband_clouds_shadows_func(qa_bandname, cloud_bit, shadow_bit, keep_orig_bands = False):
  """Get function to extract clouds and shadows mask from QA band
  The Landsat datasets has a pixel_qa band which includes a mask for clouds 
  and cloud shadows. 
  
  The function returns an image with three bands: 'clouds', 'shadows', and 'clouds_shadows'.

  Parameters
  ----------
  qa_bandname : str
  String containing the bandname with qa pixels
  cloud_bit : int
  shadow_bit : int
  Position of the binary digit containing the corresponding mask
  keep_orig_bands : bool
  If True, keep bands from original image

  Parameters for generated function:
  ----------------------------------
  im : ee.Image
  Must contain QA band with masks associated with bits
  
  Examples
  --------
  import rs
  import ee
  watershed_pt = ee.Geometry.Point([76, 12.9]) # Cauvery basin

  oli8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
    .filterBounds(watershed_pt)

  get_qaband_clouds_shadows = rs.get_qaband_clouds_shadows_func(
        qa_bandname = 'QA_PIXEL', 
        cloud_bit = 3, 
        shadow_bit = 4) 
  oli8_clouds = oli8.map(get_qaband_clouds_shadows)

  import geemap
  Map = geemap.Map()
  Map.addLayer(oli8_clouds.first(),{'bands':['clouds','clouds','clouds_shadows'],'min': 0, 'max': 1}, 'clouds')
  Map.addLayer(oli8.first(), {'bands':['SR_B5','SR_B3','SR_B4'], 'min': 0, 'max':30000},'oli8')
  Map.centerObject(oli8.first().geometry(), 8)
  Map.addLayerControl()
  """
  
  if keep_orig_bands:
      
      def get_qaband_clouds_shadows(im):
        shadows = extractQABits(im.select(qa_bandname),shadow_bit,shadow_bit).rename('shadows')
        clouds = extractQABits(im.select(qa_bandname),cloud_bit,cloud_bit).rename('clouds')
        cloud_shadow_mask = clouds.add(shadows).rename('clouds_shadows')
        return im.addBands(clouds).addBands(shadows).addBands(cloud_shadow_mask)
  else:
      
      def get_qaband_clouds_shadows(im):
        shadows = extractQABits(im.select(qa_bandname),shadow_bit,shadow_bit).rename('shadows')
        clouds = extractQABits(im.select(qa_bandname),cloud_bit,cloud_bit).rename('clouds')
        cloud_shadow_mask = clouds.add(shadows).rename('clouds_shadows')
        return im.select(None).addBands(clouds).addBands(shadows).addBands(cloud_shadow_mask)

  return get_qaband_clouds_shadows


def get_month_year_averages(imcollection, raster_bands, start_year, end_year):
  """Get average image for each month-year combination

  Parameters
  ----------
  imcollection : ee.ImageCollection
  - image collection over which to calculate monthly average
  raster_bands : list
  - list of the band names over which to take average
  start_year : int
  - start year of the range
  end_year : int
  - end year of the range
  """

  months = ee.List.sequence(1,12)
  year_seq = ee.List.sequence(start_year,end_year)

  # for each year and month, combine into a list (2 levels) of dates
  def get_all_dates(yr):
    def get_month_dates(mo):
      return(ee.Date.fromYMD(yr,mo,1))

    month_dates = months.map(get_month_dates)
    return(month_dates)

  all_dates = ee.List(year_seq.map(get_all_dates)).flatten()
  
  # Generate image collection where each image is the avg anomaly each month

  # Get Image for each month, map over all year-month combinations
  def get_31day_mean(date):
    ee_mo = ee.Date(date)
    raster_year_mo = imcollection \
      .select(raster_bands) \
      .filterDate(date,ee.Date(date).advance(1,'month')) \
      .mean() \
      .set('year',ee_mo.get('year')) \
      .set('month',ee_mo.get('month')) \
      .set('yr_mo',ee.Number(ee_mo.get('year')) \
        .add(ee.Number(ee_mo.get('month')).divide(100))) \
      .set('system:time_start', date)
    return(raster_year_mo)

  raster_all_months = all_dates.map(get_31day_mean)
  ic_all = ee.ImageCollection(raster_all_months)#.toBands()

  return ic_all

def get_monthly_average_ic(imcollection, raster_band, start_year, end_year, average_of_unit = 'months'):
  """
  6. Get monthly average values for a band of an image collection
  This function calculates the average pixel value for each month of an image collection
  in the range start_year, end_year. It returns an image with 12 bands containing the 
  average value for each month (band 0 - January, etc.).

  Parameters
  ----------
  imcollection : image collection over which to calculate
  - raster_band : name of the band in the imcollection 
  start_year : int
  - start year of the range
  end_year : int
  - end year of the range
  average_of_unit : str
  - Ether 'months' or 'images'. If 'months', first take average of images in
  each month-year, then take average over years. If 'images', directly take
  average (for each month) of all images over all years.
  """

  raster_bands = [raster_band]
  months = ee.List.sequence(1,12)

  if average_of_unit == 'months':
    year_seq = ee.List.sequence(start_year,end_year)
    # print('year_seq',year_seq)
    
    
    # for each year and month, combine into a list (2 levels) of dates
    def get_all_dates(yr):
      def get_month_dates(mo):
        return(ee.Date.fromYMD(yr,mo,1))

      month_dates = months.map(get_month_dates)
      return(month_dates)

    all_dates = ee.List(year_seq.map(get_all_dates)).flatten()
    
    # Generate image collection where each image is the avg anomaly each month

    # Get Image for each month, map over all year-month combinations
    def get_31day_mean(date):
      ee_mo = ee.Date(date)
      raster_year_mo = imcollection \
        .select(raster_bands) \
        .filterDate(date,ee.Date(date).advance(1,'month')) \
        .mean() \
        .set('year',ee_mo.get('year')) \
        .set('month',ee_mo.get('month')) \
        .set('yr_mo',ee.Number(ee_mo.get('year')) \
          .add(ee.Number(ee_mo.get('month')).divide(100))) \
        .set('system:time_start', date)
      return(raster_year_mo)

    raster_all_months = all_dates.map(get_31day_mean)
    ic_all = ee.ImageCollection(raster_all_months)#.toBands()

    # reduce each month to single image
    def reduce_ic_for_month(mo):
      ic_for_mo = ic_all \
        .filter(ee.Filter.eq('month', mo)) 
      ras_mo = ic_for_mo \
        .mean() \
        .set('month',mo) \
        .set('num_obs', ic_for_mo.size())
        
      return(ras_mo)

  elif average_of_unit == 'images':
    ic_all = imcollection.select(raster_bands)


    # reduce each month to single image
    def reduce_ic_for_month(mo):
      ic_for_mo = ic_all \
        .filter(ee.Filter.calendarRange(mo, mo, 'month'))
      ras_mo = ic_for_mo \
        .mean() \
        .set('month',mo) \
        .set('num_obs', ic_for_mo.size())
        
      return(ras_mo)

  else:
    raise ValueError('average_of_unit is \'' + average_of_unit + 
      '\' but must be either \'months\' or \'images\'')



  raster_monthly_ic = ee.ImageCollection(months.map(reduce_ic_for_month))

  return raster_monthly_ic


def get_monthly_average(imcollection, reduce_region, raster_band, start_year, end_year, scale_m):
  """
  6. Get monthly average values for a region
  This function calculates the average pixel value for each month of an image collection
  in the range start_year, end_year. It then uses reduceRegion to calculate the overall average
  within a region. It returns a 12x1 array of monthly values.
  imcollection - image collection over which to calculate
  reduce_region - feature over which to calculate the average
  raster_band - name of the band in the imcollection 
  start_year - start year of the range
  end_year - end year of the range
  scale_m - scale for reduceRegion call
  monthly_ts = rs.get_monthly_average(modisET, region, 'PET', 2001, 2005, 10000)
  """

  raster_bands = [raster_band]
  year_seq = ee.List.sequence(start_year,end_year)
  # print('year_seq',year_seq)
  
  months = ee.List.sequence(1,12)
  
  # for each year and month, combine into a list (2 levels) of dates
  def get_all_dates(yr):
    def get_month_dates(mo):
      return(ee.Date.fromYMD(yr,mo,1))

    month_dates = months.map(get_month_dates)
    return(month_dates)

  all_dates = ee.List(year_seq.map(get_all_dates)).flatten()
  
  # Generate image collection where each image is the avg anomaly each month

  # Get IC for each month
  def reduce_ic_toMonths(mo):
    ee_mo = ee.Date(mo)
    raster_mo = imcollection \
      .select(raster_bands) \
      .filterDate(mo,ee.Date(mo).advance(1,'month')) \
      .mean() \
      .set('year',ee_mo.get('year')) \
      .set('month',ee_mo.get('month')) \
      .set('yr_mo',ee.Number(ee_mo.get('year')) \
        .add(ee.Number(ee_mo.get('month')).divide(100)))
    return(raster_mo)
  
  raster_months = all_dates.map(reduce_ic_toMonths)
  
  raster_months_ic = ee.ImageCollection(raster_months)
  
  # reduce to region
  def reduce_ic_toRegion(mo):
    ras_mo = raster_months_ic.select(raster_bands) \
      .filter(ee.Filter.eq('month',mo)) \
      .reduce(ee.Reducer.mean()) \
      .set('month',mo) \
      .reduceRegion(ee.Reducer.mean(),reduce_region, scale_m) \
      .get(raster_band + '_mean')
      
    return(ras_mo)
  raster_monthly_ts = months.map(reduce_ic_toRegion)
  
  monthly_array = ee.Array(raster_monthly_ts)
  
  return(monthly_array)



def get_monthly_water_balance(reduce_region, start_year, end_year, scale_m):
  """
  7. Get approximate monthly water balance for a region using CHIRPS and MODIS ET
  This function calculates the average P, PET, and ET for a region
  in the range start_year, end_year. It uses reduceRegion to calculate the overall average
  within the region. It returns a 12x3 array of monthly values.
  reduce_region - feature over which to calculate the average
  start_year - start year of the range
  end_year - end year of the range
  scale_m - scale for reduceRegion call
  plot_vars = rs.get_monthly_water_balance(reduce_region, start_year, end_year, scale_m)
  """
  
  modisET = ee.ImageCollection("MODIS/006/MOD16A2")
  chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
  
  pet_monthly = rs.get_monthly_average(modisET, reduce_region, 'PET', start_year, end_year, scale_m) \
    .multiply(0.01).multiply(30)
  
  et_monthly = get_monthly_average(modisET, reduce_region, 'ET', start_year, end_year, scale_m) \
    .multiply(0.01).multiply(30)
  
  P_monthly = get_monthly_average(chirps, reduce_region, 'precipitation', start_year, end_year, scale_m) \
    .multiply(30)
  
  plot_vars = ee.Array.cat([P_monthly, pet_monthly, et_monthly],1)
  
  return(plot_vars)



def plot_monthly_water_balance(reduce_region, start_year, end_year, scale_m, title):
  """
  8. Get a timeseries plot of monthly water balance for a region using CHIRPS and MODIS ET
  This function calculates the average P, PET, and ET for a region
  in the range start_year, end_year. It uses reduceRegion to calculate the overall average
  within the region. It returns a 12x3 array of monthly values.
  reduce_region - feature over which to calculate the average
  start_year - start year of the range
  end_year - end year of the range
  scale_m - scale for reduceRegion call
  title - title for plot
  wb_plot = rs.get_monthly_water_balance(reduce_region, start_year, end_year, scale_m, 'Plot title')
  print(wb_plot)
  """

  # Monthly time series and chart
  plot_vars = get_monthly_water_balance(reduce_region, start_year, end_year, scale_m)
  
  # print('plot_vars',plot_vars)
  months = ee.List.sequence(1,12)
  ts_chart = ui.Chart.array.values(plot_vars, 0, months) \
    .setSeriesNames(['CHIRPS P', 'MODIS PET','MODIS ET']) \
    .setSeriesNames(['P', 'PET','ET']) \
    .setOptions({
      'title': title,
      'hAxis': {
        'title': 'Month'
      },
      'vAxis': {title: 'Value, mm'},
      'pointSize': 5,                # <--- hide the points.
      'lineSize': 3                  # <--- show the line.
    })
    
  return(ts_chart)



def set_feature_id_func(field_name = 'pt_id'):
  """
  This returns a function whic sets a point id for a feature in a collection
  call as FeatureCollection.map(set_feature_id_func('pt_id'))

  Parameters
  ----------
  field_name : str
  - Name of the field in which to add the id

  Examples
  --------
  listOfFeatures = [
    ee.Feature(ee.Geometry.Point(77,13), {'key' : 'pt1'}),
    ee.Feature(ee.Geometry.Point(77,13.1), {'key' : 'pt2'}),
    ee.Feature(ee.Geometry.Point(77,13.2), {'key' : 'pt3'})
  ]
  fc = ee.FeatureCollection(listOfFeatures)
  fc.map(rs.set_feature_id_func(field_name = 'pt_id')).getInfo()
  """
  def setPt(pt):
    return pt.set(field_name,ee.Number.parse(pt.id()).int()) 
  return setPt


def get_grid_pts_func(raster, square_px, scale):
  """Get grid points in raster
  This returns a function which converts a point feature to a square grid of 
  points associated with a raster. This approach allows the functiont to be mapped
  over a FeatureCollection.
  It uses eeGeometry.buffer().bounds() to generate a square in the vicinity
  of the location. It then uses rs.get_raster_pts() to generate the grid of points.
  It expects tha the original point contains a metadata feature 'loc_id' which is added
  to the metadata of the gridded points. The gridded points are labelled
  with 'pt_id'.

  Parameters
  ----------
  raster : ee.Image
  - raster from which to sample pixels
  square_px : numeric
  - The number of pixels along each dimension of the grid
  scale : numeric
  - The nominal size of the raster pixels

  Examples
  --------
  srtm = ee.Image('CGIAR/SRTM90_V4')
  pt = ee.Feature(ee.Geometry.Point([77,13])).set('loc_id',1)
  square_px = 4
  nominal_scale = 900 # this is the width of a pixel in the image

  # generate the function
  get_loc_grid_pts = rs.get_grid_pts_func(srtm, square_px, nominal_scale)

  # Apply the function to a the feature
  loc_grid_pts = get_loc_grid_pts(pt)

  # loc_grid_pts.getInfo()

  # # Example mapping over a feature collection:
  sample_rect_pts = loc_grid_pts.map(rs.get_grid_pts_func(srtm, square_px = 4, scale = 90)).flatten()

  # # Map the output
  Map = geemap.Map()
  Map.addLayer(sample_rect_pts, {}, 'Grid of grid')
  Map.addLayer(loc_grid_pts.draw(color = 'red', pointRadius = 5),{},'Main grid')
  Map.centerObject(loc_grid_pts, 14)
  Map.addLayerControl()
  Map
  """
  # Define the function to be returned
  def gen_grid_pts(feat_pt):
    buffer_dist = square_px / 2 * scale
    sample_rect = feat_pt.buffer(buffer_dist).bounds().geometry()
    grid_pts_prep = get_raster_pts(raster = raster, region = sample_rect, scale = scale)

    loc_id = feat_pt.get('loc_id')

    def add_loc_id(ft): 
      return ft.set('loc_id',loc_id)

    grid_pts = grid_pts_prep.map(set_feature_id_func(field_name = 'pt_id')).map(add_loc_id)
    
    return grid_pts

  return gen_grid_pts

import pandas as pd

def ee_array_to_df(arr, list_of_bands):
    """Transforms client-side ee.Image.getRegion array to pandas.DataFrame.
    From https://developers.google.com/earth-engine/tutorials/community/intro-to-python-api
    """
    df = pd.DataFrame(arr)

    # Rearrange the header.
    headers = df.iloc[0]
    df = pd.DataFrame(df.values[1:], columns=headers)

    # Remove rows without data inside.
    # df = df[['longitude', 'latitude', 'time', *list_of_bands]].dropna()

    numeric_bands = ['id','longitude','latitude'] + list_of_bands
    # Convert the data to numeric values.
    for band in numeric_bands:
        df[band] = pd.to_numeric(df[band], errors='coerce')

    # Convert the time field into a datetime.
    df['datetime'] = pd.to_datetime(df['time'], unit='ms')

    # Keep the columns of interest.
    # df = df[['time','datetime',  *list_of_bands]]

    return df

