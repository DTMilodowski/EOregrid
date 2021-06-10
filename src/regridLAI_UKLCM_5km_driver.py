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
dx_target = 0.05
dy_target = 0.05
N = 61.25
S = 49.75
E = 2.25
W = -10.5
extent = [N,S,E,W]
# Set up directory structure as required
path2orig = glob.glob('/disk/scratch/local.2/copernicus/LAI_300m_2014-2018/*')
path2orig += glob.glob('/disk/scratch/local.2/copernicus/LAI_300m_2019-2020/ftp.copernicus.vgt.vito.be/M0083177/LAI300_2*PROBAV*')
path2dest = '/disk/scratch/local.2/copernicus/LAI_300m_2014_2020_UK_tiles/'
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

# open land use file
lc_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m.tif'
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
#lc_agg_file_300m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb300m_aggregated.tif'
lc_agg_file_300m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb300m_aggregated.nc'
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
dx = 0.00297619#04762040103#0.002777777777778
dy = -0.00297619#04761897995#0.002777777777778
W_ = W-dx/2.; E_ = E+dx/2.
N_ = N+dy/2.; S_ = S-dy/2.

template_nc = glob.glob('%s/*nc' % path2orig[0])[0]
template = xr.open_dataset(template_nc).sel(lat=slice(N_,S_),lon=slice(W_,E_))
rows,cols = template.LAI.shape
dx = np.diff(template.lon.values[:2])
dy = np.diff(template.lat.values[:2])
W_ = template.lon.values.min()-dx/2.
E_ = template.lon.values.max()+dx/2.
N_ = template.lat.values.max()+np.abs(dy)/2.
S_ = template.lat.values.min()-np.abs(dy)/2.
os.system('gdalwarp -te %f %f %f %f -r mode -s_srs EPSG:27700 -t_srs EPSG:4326 \
            -ts %i %i -dstnodata 0 -ot Float64 -overwrite %s %s' %
            (W_,S_,E_,N_, cols, rows,lc_agg_file_25m,lc_agg_file_300m))

# load in the 300m resolution raster
#lc_ref = xr.open_rasterio(lc_agg_file_300m)[0].sel(y=slice(N,S),x=slice(W,E))
lc_ref = xr.open_dataset(lc_agg_file_300m).Band1
lc_ref=lc_ref.assign_coords(lat = lc_ref.lat.values[::-1],lon = lc_ref.lon.values)
lc_ref.values = lc_ref.values[::-1,:]
lc_ref = lc_ref[:,:-1]
lc_id = [1,2,3,4,5,6,10,'_']
for ii,lc in enumerate(landcover):
    path2dest_sub = '%s%s/' % (path2dest, lc)
    if lc == 'non_forest':
        mask = lc_ref.values>2.5
    else:
        mask = lc_ref.values==lc_id[ii]
    for jj,dir in enumerate(path2orig):
        print('processing layer %i of %i, for %s' % (jj+1,len(path2orig),lc))
        input_nc = glob.glob('%s/*nc' % dir)[0]
        output_nc = '%s%s_%s_5km.nc' % (path2dest_sub,input_nc.split('/')[-1][:-3],lc)
        regrid_ds = regridLAI.regridLAI(input_nc, output_nc, dx_target, dy_target,
                    mask=mask, subset_label=lc, extent=extent,
                    variables=['LAI','RMSE'], aggregation_mode = ['mean','mean'],
                    projected=False)
