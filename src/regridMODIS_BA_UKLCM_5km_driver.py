"""
regdridMODIS_BA_UKLCM_5km_driver.py
--------------------------------------------------------------------------------
Regrid MODIS burned area to coarser resolution for assimilation, but tiling
according to CEH LCM land cover.
"""
import os
import sys
import glob
import numpy as np
import xarray as xr
import osgeo.gdal as gdal
import matplotlib.pyplot as plt
#import regridLAI as regridLAI
sys.path.append('./data_io/')
import data_io as io
"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- MODIS Burned Area data (MCD64A1)
- CEH 2015 UK land cover map
"""
# Bounding extent of desired domain in WGS84 and MODIS sinusoidal
W = -10.5; Ws = -754380.789747
E = 2.25; Es = 120353.6190577065
N = 61.25; Ns = 6810696.933569999
S = 49.75; Ss = 5532417.148668328

dx_target = 0.05
dy_target = 0.05

# Set up directory structure as required
path2orig = '/home/dmilodow/DataStore_GCEL/BurnedArea/MCD64A1/rawHDF/'
path2merged = '%s/merged_files/' % path2orig
path2dest = '/disk/scratch/local.2/copernicus/MODIS_BurnedArea_2014_2018_UK_tiles/'
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

# MODIS files
doy_1_to_month = {'001':'01', '032':'02', '060':'03', '091':'04', '121':'05',
                    '152':'06', '182':'07', '213':'08', '244':'09', '274':'10',
                    '305':'11', '335':'12'}
doy_2_to_month = {'001':'01', '032':'02', '061':'03', '092':'04', '122':'05',
                    '153':'06', '183':'07', '214':'08', '245':'09', '275':'10',
                    '306':'11', '336':'12'}
modis_temp = []
for tile in glob.glob('%s/*hdf' % path2orig):
    modis_temp.append(tile.split('/')[-1][9:16])

modis_unique = np.unique(modis_temp)

modis_months = []
modis_doy = []
for date in modis_unique:
    year = date[:4]
    doy = date[4:]
    modis_doy.append('MCD64A1.A%s%s' % (year,doy))
    try:
        modis_months.append('MCD64A1.A%s%s' % (year,doy_1_to_month[doy]))
    except:
        try:
            modis_months.append('MCD64A1.A%s%s' % (year,doy_2_to_month[doy]))
        except:
            print("can't deal with this DoY: %s" % doy)

# merge the MODIS BA for the extent of the UK
for jj,month in enumerate(modis_doy):
    input_tile_string = ''
    for tile in glob.glob('%s/%s*hdf' % (path2orig,month)):
        input_tile_string+='HDF4_EOS:EOS_GRID:"%s:MOD_Grid_Monthly_500m_DB_BA:Burn Date" ' % tile
    # merge modis tiles for region of interest
    os.system('gdal_merge.py -o %s/%s_BA_UK.tif \
                -ul_lr %f %f %f %f %s'
                % (path2merged,modis_months[jj], Ws, Ns, Es, Ss, input_tile_string))

# open land use file
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_MODIS = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gbMODIS_aggregated.tif'

# MODIS Sinusoidalprojection WKT
modis_wkt = "PROJCS[\"unnamed\",GEOGCS[\"Unknown datum based upon the custom spheroid\",DATUM[\"Not_specified_based_on_custom_spheroid\",SPHEROID[\"Custom spheroid\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"longitude_of_center\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]]]"


"PROJCS[\"unnamed\",
        GEOGCS[\"Unknown datum based upon the custom spheroid\",
                DATUM[\"Not_specified_based_on_custom_spheroid\",
                        SPHEROID[\"Custom spheroid\",6371007.181,0]],
                PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],
        PROJECTION[\"Sinusoidal\"],
        PARAMETER[\"longitude_of_center\",0],
        PARAMETER[\"false_easting\",0],
        PARAMETER[\"false_northing\",0],
        UNIT[\"metre\",1,
            AUTHORITY[\"EPSG\",\"9001\"]]]"

# load tempate MODIS raster
modis = xr.open_rasterio('%s/%s_BA_UK.tif' % (path2merged,modis_months[0])).sel(band=1)
dy_modis = modis.y.values[1]-modis.y.values[0]
dx_modis = modis.x.values[1]-modis.x.values[0]
modis=None

# Convert LCM2015 to sinusoidal projection at 250m resolution (to match the
# MODIS data being regridded). Approach is to sample land cover prior to
# regridding to the final coarse lat/lon grid
ds = gdal.Warp('%s' % (lc_agg_file_MODIS),'%s' % (lc_agg_file_25m),
                srcSRS='EPSG:27700', dstSRS=modis_wkt,resampleAlg='mode',
                outputBounds=(W,S,E,N) , outputBoundsSRS='EPSG:4326',
                xRes=dx_modis, yRes=dy_modis)
del ds

lc_modis = xr.open_rasterio(lc_agg_file_MODIS).sel(band=1)[:-1]
lc_id = [1,2,3,4,5,6]

# regridding BA with tiling approach
for jj,month in enumerate(modis_doy):
    BA = xr.open_rasterio('%s/%s_BA_UK.tif' % (path2merged,modis_months[jj])).sel(band=1)
    for ii,lc in enumerate(landcover):
        path2dest_sub = '%s%s/' % (path2dest, lc)
        print('processing layer %i of %i, for %s' % (jj+1,len(modis_doy),lc))

        BA_lc = BA.copy(deep=True)
        BA_lc.values=np.zeros(BA.values.shape,dtype='float')
        BA_lc.values[BA.values>0]=1.
        BA_lc.values[lc_modis.values!=lc_id[ii]]=np.nan # ignore areas outside land cover

        # write to temporary array
        sinusoidal_file = '%s%s_BA_UK_500m_sinusoidal_%s.tif' % (path2merged,modis_months[jj],lc)
        io.write_xarray_to_GeoTiff(BA_lc,sinusoidal_file,wkt=modis_wkt)
        # regrid to wgs84 template
        wgs84_file = '%s%s_BA_%s_5km.tif' % (path2dest_sub,modis_months[jj],lc)
        ds = gdal.Warp('%s' % wgs84_file,'%s' % sinusoidal_file,
                srcSRS=modis_wkt, dstSRS='EPSG:4326',resampleAlg='average',
                outputBounds=(W,S,E,N) , outputBoundsSRS='EPSG:4326',
                xRes=dx_target, yRes=dy_target)
        del ds
