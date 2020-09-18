"""
LAI_test_drive.py
--------------------------------------------------------------------------------
A test drive of the LAI regridding script. Will start by taking a land cover
mask from ... and regridding Copernicus 300m LAI based on the forest class.

Let's see how it goes!
"""
import os
import sys
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import regridLAI as regridLAI

sys.path.append('./data_io/')
import data_io as io
"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- Copernicus 300 m LAI data
- CEH 2015 UK land cover map
"""

# Set up directory structure as required
path2orig = glob.glob('/disk/scratch/local.2/copernicus/LAI_300m_2014-2018/*')
path2dest = '/disk/scratch/local.2/copernicus/LAI_300m_2014_2018_UK_tiles/'
try:
    os.mkdir(path2dest)
    print('path2dest not found, creating new directory %s' % path2dest)
except:
    print('path2dest is %s' % path2dest)


landcover = ['broadleaf_forest', 'coniferous_forest', 'arable', 'improved_grass',
            'seminatural_grass', 'heath']
for lc in landcover:
    path2dest_sub = '%s%s/' % (path2dest, lc)
    try:
        os.mkdir(path2dest_sub)
        print('subdirectory not found, creating new directory %s' % path2dest_sub)
    except:
        print('subdirectory is %s' % path2dest_sub)

# open land use file
lc_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m.tif'
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_300m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb300m_aggregated.tif'
lc_25m = xr.open_rasterio(lc_file_25m).sel(band=1)

lc_25m_agg = lc_25m.copy(deep=True)
lc_25m_agg.values[(lc_25m.values>=5)*(lc_25m.values<=8)]=5
lc_25m_agg.values[(lc_25m.values>=9)*(lc_25m.values<=12)]=6
lc_25m_agg.values[lc_25m.values==13]=7
lc_25m_agg.values[lc_25m.values==14]=8
lc_25m_agg.values[(lc_25m.values>=15)*(lc_25m.values<=19)]=9
lc_25m_agg.values[(lc_25m.values>=20)*(lc_25m.values<=21)]=10

io.write_xarray_to_GeoTiff(lc_25m_agg,lc_agg_file_25m)

# regrid aggregated LCM to 300m Copernicus grid
dx = 0.0029761904762040103#0.002777777777778
dy = 0.0029761904761897995#0.002777777777778
cop300lat = np.arange(-60.,80,dy)
cop300lon = np.arange(-180.,180,dx)
UK300lat = cop300lat[(cop300lat>49.75)*(cop300lat<61.25)]
UK300lon = cop300lon[(cop300lon>-10.5)*(cop300lon<2.25)]
os.system('gdalwarp -te %f %f %f %f -r mode -s_srs EPSG:27700 -t_srs EPSG:4326 \
            -tr %f %f -dstnodata 0 -overwrite %s %s' %
            (UK300lon[0]-dx/2., UK300lat[0]-dy/2., UK300lon[-1]+dx/2, UK300lat[-1]+dy/2, dx, dy,
            lc_agg_file_25m,lc_agg_file_300m))

# load in the 300m resolution raster
lc_300m = xr.open_rasterio(lc_agg_file_300m)[0]
lc_id = [1,2,3,4,5,6]
dx_target = 0.05
dy_target = 0.05
extent = [UK300lat[-1]+dy/2,UK300lat[0]-dy/2,UK300lon[-1]+dx/2,UK300lon[0]-dx/2] # N,S,E,W
for ii,lc in enumerate(landcover):
    path2dest_sub = '%s%s/' % (path2dest, lc)
    for jj,dir in enumerate(path2orig):
        print('processing layer %i of %i, for %s' % (jj+1,len(path2orig),lc))
        input_nc = glob.glob('%s/*nc' % dir)[0]
        output_nc = '%s%s_%s_5km.nc' % (path2dest_sub,input_nc.split('/')[-1][:-3],lc)
        regrid_ds = regridLAI.regridLAI(input_nc, output_nc, dx_target, dy_target,
                    mask=(lc_300m==lc_id[ii]), subset_label=lc, extent=extent,
                    variables=['LAI','RMSE'], aggregation_mode = ['mean','quadrature'],
                    projected=False)
