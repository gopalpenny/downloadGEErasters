#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 22:27:11 2023

@author: gopal
"""
import ee

def get_s2_sr_cld_col(aoi, params):
    """Join Sentinel-2 Surface Reflectance and Cloud Probability

    This function retrieves and joins ee.ImageCollections:
    'COPERNICUS/S2_SR' and 'COPERNICUS/S2_CLOUD_PROBABILITY'

    Parameters
    ----------
    aoi : ee.Geometry or ee.FeatureCollection
      Area of interested used to filter Sentinel imagery
    params : dict
      Dictionary used to select and filter Sentinel images. Must contain
      START_DATE : str (YYYY-MM-DD)
      END_DATE : str (YYYY-MM-DD)
      CLOUD_FILTER : int
        Threshold percentage for filtering Sentinel images
    """


    # start_date = params['START_DATE']
    # end_date = params['END_DATE']
    # Import and filter S2 SR.
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(aoi))
        # .filterDate(start_date, end_date)
        # .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', params['CLOUD_FILTER'])))

    # Import and filter s2cloudless.
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY') 
        .filterBounds(aoi))
        # .filterDate(start_date, end_date))

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))


def add_cloud_bands(img, params):
    """Add cloud bands to Sentinel-2 image

    Parameters
    ----------
    img : ee.Image
      Sentinel 2 image including (cloud) 'probability' band
    params : dict
      Parameter dictionary including
      CLD_PRB_THRESH : int
        Threshold percentage to identify cloudy pixels

    """
    # Get s2cloudless image, subset the probability band.
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(params['CLD_PRB_THRESH']).rename('clouds')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def add_shadow_bands(img, params):
    """Add cloud shadow bands to Sentinel-2 image

    Parameters
    ----------
    img : ee.Image
      Sentinel 2 image including (cloud) 'probability' band
    params : dict
      Parameter dictionary including
      NIR_DRK_THRESH : int
        Threshold percentage to identify potential shadow pixels as dark pixels from NIR band
      CLD_PRJ_DIST : int
        Distance to project clouds along azimuth angle to detect potential cloud shadows
    """

    # Identify water pixels from the SCL band.
    not_water = img.select('SCL').neq(6)

    # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    SR_BAND_SCALE = 1e4
    dark_pixels = img.select('B8').lt(params['NIR_DRK_THRESH']*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, params['CLD_PRJ_DIST']*10)
        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    # Identify the intersection of dark pixels with cloud shadow projection.
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shadow_mask_func(params):
  """Add cloud/shadow mask to Sentinel-2 imagery

  This is a helper function to allow for additional parameters in a 'map' command.
  It returns a function which take a single ee.Image input. This function can be called
  as as ImCollection.map(add_cld_shadow_mask_func(params)). See
  https://gis.stackexchange.com/questions/302760/gee-imagecollection-map-with-multiple-input-function
  for how this works.

  The returned function (add_cld_shadow_mask) adds a mask for clouds and shadows to Sentinel-2 imagery. The 
  cloud mask is determined based on a threshold value for clouds

  Parameters
  ----------
  params : dict
    Parameter dictionary including
    CLD_PRB_THRESH : int
      Threshold percentage to identify cloudy pixels
    NIR_DRK_THRESH : int
      Threshold percentage to identify potential shadow pixels as dark pixels from NIR band
    CLD_PRJ_DIST : int
      Distance to project clouds along azimuth angle to detect potential cloud shadows.
      Can be 1 for low cumulus clouds, 2.5 or greater for mountainous regions or high altitude clouds.
    BUFFER : int
      Distance to buffer combined cloud and shadow mask

  Parameters for returned function (add_cld_shadow_mask)
  ------------------------------------------------------
    img : ee.Image()
  """

  def add_cld_shdw_mask(img):
    # Add cloud component bands.
    img_cloud = add_cloud_bands(img, params)

    # Add cloud shadow component bands.
    img_cloud_shadow = add_shadow_bands(img_cloud, params)

    # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(params['BUFFER']*2/20)
        .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)

  return add_cld_shdw_mask