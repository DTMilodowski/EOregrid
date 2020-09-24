"""
regridAGB_UKLCM_5km_driver.py
--------------------------------------------------------------------------------
Regridding script to regrid ESA-BIOMASS CCI to 5km, stratified by aggregated
land cover indicated by the CEH LCM dataset.
"""
import os
import sys
import glob
import numpy as np
import xarray as xr

import geospatial_tools as gst

sys.path.append('./data_io/')
import data_io as io

"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- ESA BIOMASS-CCI AGB AND SD
- CEH 2015 UK land cover map
"""

# Set up directory structure as required
path2agb = '/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/geotiff/'
path2dest =  '/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/UKLCM_5km/'

try:
    os.mkdir(path2dest)
    print('path2dest not found, creating new directory %s' % path2dest)
except:
    print('path2dest is %s' % path2dest)

landcover = ['broadleaf_forest', 'coniferous_forest', 'arable', 'improved_grass',
            'seminatural_grass', 'heath']

# open land use file
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_BIOMASS = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb_ESACCI_BIOMASS_aggregated.tif'
lc_25m_agg = xr.open_rasterio(lc_agg_file_25m).sel(band=1)

# open AGB file
agb = xr.open_rasterio('%sN80W020_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv1.0.tif' % path2agb)[0]
agb_sd = xr.open_rasterio('%sN80W020_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-2017-fv1.0.tif' % path2agb)[0]

# regrid aggregated LCM
dy_orig = agb.y.values[1]-agb.y.values[0]
dx_orig = agb.x.values[1]-agb.x.values[0]
UKlat = agb.y.values[(agb.y.values>49.75)*(agb.y.values<61.25)]
UKlon = agb.x.values[(agb.x.values>-10.5)*(agb.x.values<2.25)]
W = UKlon[0]-dx_orig/2.
E = UKlon[-1]+dx_orig/2.
N = UKlat[0]-dy_orig/2.
S = UKlat[-1]+dy_orig/2.

os.system('gdalwarp -te %f %f %f %f -r mode -s_srs EPSG:27700 -t_srs EPSG:4326 \
            -tr %f %f -dstnodata 0 -overwrite %s %s' %
            (W, S, E, N, dx_orig, dy_orig, lc_agg_file_25m,lc_agg_file_BIOMASS))

# load in the new raster
lc_ref = xr.open_rasterio(lc_agg_file_BIOMASS)[0]
lc_id = [1,2,3,4,5,6]
dx_target = 0.05
dy_target = 0.05

agb_ref=agb.sel(x=slice(W,E),y=slice(N,S))
agb_sd_ref=agb_sd.sel(x=slice(W,E),y=slice(N,S))

# Reset grid
Yorig = agb_ref.y.values
Xorig = agb_ref.x.values

# define scanning window size
Ysize = np.abs(np.round(dy_target/dy_orig).astype('i'))
Xsize = np.abs(np.round(dx_target/dx_orig).astype('i'))

# define destination lat / lon arrays
Ydest = np.arange(N-dy_target/2.,S,-dy_target)
Xdest = np.arange(W+dx_target/2.,E,dx_target)

# calculate grid cell area
areas = gst.calculate_cell_areas(Yorig,Xorig,projected=False)

for ii,lc in enumerate(landcover):
        print("Regridding variable AGB for %s...             " % lc)
        # clip variable grid using lat_mask and lon_mask
        mask = lc_ref.values==lc_id[ii]
        agb_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest, Ysize, Xsize,
                                        areas, agb_ref.values, mask=mask,
                                        aggregation_mode='mean')
        # write to file
        agb_out = xr.DataArray(data=agb_regrid, coords={'x':Xdest,'y':Ydest}, dims=['y','x'],
                                attrs={'details':'regridded AGB for %s, based on ESACCI-BIOMASS and CEH LCM2015' % lc,
                                        'name':'AGB for %s' % lc,
                                        'units':'Mg/ha'})
