"""
regridGFW_UKLCM_5km_driver.py
--------------------------------------------------------------------------------
regrid tree cover loss to 5km based on land-cover subsets defined by CEH LCM2015

"""
import os
import numpy as np
import xarray as xr
import geospatial_tools as gst

"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- Copernicus 300 m LAI data
- CEH 2015 UK land cover map
"""

# Set up directory structure as required
path2orig = '/disk/scratch/local.2/GFW/'
path2dest = '/disk/scratch/local.2/GFW/UKLCM_5km_GFW_2000_2019/'
try:
    os.mkdir(path2dest)
    print('path2dest not found, creating new directory %s' % path2dest)
except:
    print('path2dest is %s' % path2dest)

landcover = ['broadleaf_forest', 'coniferous_forest', 'arable', 'improved_grass',
            'seminatural_grass', 'heath', 'built_up', 'non_forest']

for lc in landcover:
    path2dest_sub = '%s%s/' % (path2dest, lc)
    try:
        os.mkdir(path2dest_sub)
        print('subdirectory not found, creating new directory %s' % path2dest_sub)
    except:
        print('subdirectory is %s' % path2dest_sub)

W = -10.5
E = 2.25
N = 61.25
S = 49.75

# merge GFW tiles for UK
os.system('gdal_merge.py -o %s/merged_files/Hansen_GFC-2019-v1.7_lossyear_UK.tif \
                -ul_lr %f %f %f %f %s/source_files/Hansen_GFC-2019-v1.7_lossyear_*N_*.tif'
                % (path2orig, W, N, E, S, path2orig))
gfw = xr.open_rasterio('%s/merged_files/Hansen_GFC-2019-v1.7_lossyear_UK.tif' % path2orig).sel(band=1)
dy_orig = gfw.y.values[1]-gfw.y.values[0]
dx_orig = gfw.x.values[1]-gfw.x.values[0]

# open land use file
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_GFW = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb_GFW_aggregated.tif'

# regrid aggregated LCM to GFW grid
os.system('gdalwarp -te %f %f %f %f -r near -s_srs EPSG:27700 -t_srs EPSG:4326 \
            -tr %f %f -dstnodata 0 -overwrite %s %s' %
            (W, S, E, N, dx_orig, dy_orig, lc_agg_file_25m, lc_agg_file_GFW))
lc_gfw = xr.open_rasterio(lc_agg_file_GFW).sel(band=1)
lc_id = [1,2,3,4,5,6,10,'_']

# set up target grid
dx_target = 0.05
dy_target = 0.05

# Reset grid
Yorig = gfw.y.values
Xorig = gfw.x.values

# define scanning window size
Ysize = np.abs(np.round(dy_target/dy_orig).astype('i'))
Xsize = np.abs(np.round(dx_target/dx_orig).astype('i'))

# define destination lat / lon arrays
Ydest = np.arange(N-dy_target/2.,S,-dy_target)
Xdest = np.arange(W+dx_target/2.,E,dx_target)

# calculate grid cell area
areas = gst.calculate_cell_areas(Yorig,Xorig,projected=False)

for year in np.arange(2001,2020):
    gfw_code = year-2000
    var = (gfw.values==gfw_code).astype('float')
    for ii,lc in enumerate(landcover):
        print("Regridding tree cover loss in %i for %s...             " % (year,lc))
        # clip variable grid using lat_mask and lon_mask
        if lc == 'non_forest':
            mask = lc_ref.values>2.5
        else:
            mask = lc_ref.values==lc_id[ii]
        # weights average
        var_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest, Ysize, Xsize,
                                        areas, var, mask=mask, aggregation_mode='mean')

        # write to file
        path2dest_sub = '%s%s/' % (path2dest, lc)
        nc_out = '%stree_cover_loss_fraction_%i_%s_5km.nc' % (path2dest_sub, year, lc)
        var_out = xr.DataArray(data={'tree_cover_loss' : var_regrid}, coords={'longitude':Xdest, 'latitude':Ydest},
                                    dims=['latitude', 'longitude'],
                                    attrs={'details':'regridded GFW tree cover loss for %s, based on GFW v1.7 and CEH LCM2015' % lc,
                                    'name':'tree cover loss fraction for %s' % lc,
                                    'units':'fraction (of land cover within pixel)'})
        var_out.to_netcdf(path=nc_out)
