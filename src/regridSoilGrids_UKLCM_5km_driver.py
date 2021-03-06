"""
regridSoilGrids_UKLCM_5km_driver.py
--------------------------------------------------------------------------------
Regridding script to regrid SoilGrids to 5km, stratified by aggregated
land cover indicated by the CEH LCM dataset.
"""
import os
import numpy as np
import xarray as xr

from osgeo import gdal

import geospatial_tools as gst

"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- SoilGrids2 soil organic carbon, sand, clay and silt %
- CEH 2015 UK land cover map
The output layers for CARDAMOM are:
- SOC within upper 100cm
- sand, silt & clay % for 0-30cm
- sand, silt & clay % for 30-100cm
"""
dx_target = 0.05
dy_target = 0.05
N = 61.25
S = 49.75
E = 2.25
W = -10.5

# Set up directory structure as required
path2sg = '/disk/scratch/local.2/SoilGrids2/files.isric.org/soilgrids/latest/data/'
path2wgs84 = '/disk/scratch/local.2/SoilGrids2/wgs84/'
path2dest =  '/disk/scratch/local.2/SoilGrids2/UKLCM_5km/'

try:
    os.mkdir(path2dest)
    print('path2dest not found, creating new directory %s' % path2dest)
except:
    print('path2dest is %s' % path2dest)

landcover = ['broadleaf_forest', 'coniferous_forest', 'arable', 'improved_grass',
            'seminatural_grass', 'heath', 'built_up', 'non_forest']
soilgrids_vars = ['sand','silt','clay','ocd']

for lc in landcover:
    path2dest_sub = '%s%s/' % (path2dest, lc)
    try:
        os.mkdir(path2dest_sub)
        print('subdirectory not found, creating new directory %s' % path2dest_sub)
    except:
        print('subdirectory is %s' % path2dest_sub)


igh = "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs" # proj string for Homolosine projection
sg_url = '/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/'

# open land use file
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_soilgrids = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb_soilgrids_aggregated.tif'
lc_25m_agg = xr.open_rasterio(lc_agg_file_25m).sel(band=1)

# process each variable in the soilgrids catalogue
# - note that the soilgrids data are provided in equal area projection
# - thus, first step will be to regrid to WGS84 for consistency with the other
#   datasets. Will use gdalwarp and the "average" mode
# - native units for soilgrids variables are:
#       - ocd: organic carbon density: hg/dm3. Divide by 100 -> kg/m3
#       - sand/silt/clay: g/kg of fines. Divide by 10 -> %

# download and regrid to wgs84
for sgvar in soilgrids_vars:
    for layer in ['0-5','5-15','15-30','30-60','60-100','100-200']:
        print('%s; soil layer %scm            ' % (sgvar,layer),end='\r')
        ds = gdal.Warp('%s/%s/%s_%scm_mean.tif' % (path2wgs84,sgvar,sgvar,layer),
                        '%s/%s/%s_%scm_mean.vrt' % (sg_url,sgvar,sgvar,layer),
                        dstSRS='EPSG:4326', srcSRS=igh, outputBounds = (W,S,E,N),
                        outputBoundsSRS='EPSG:4326', resampleAlg='average')
        del ds

        if sgvar=='ocd':
            for quantile in ['Q0.05','Q0.95']:
                ds = gdal.Warp('%s/%s/%s_%scm_%s.tif' % (path2wgs84,sgvar,sgvar,layer,quantile),
                                '%s/%s/%s_%scm_%s.vrt' % (sg_url,sgvar,sgvar,layer,quantile),
                                dstSRS='EPSG:4326', srcSRS=igh, outputBounds = (W,S,E,N),
                                outputBoundsSRS='EPSG:4326', resampleAlg='average')
                del ds

# Now regrid
for ii,sgvar in enumerate(soilgrids_vars):

    var_000_005 = xr.open_rasterio('%s%s/%s_0-5cm_mean.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    var_005_015 = xr.open_rasterio('%s%s/%s_5-15cm_mean.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    var_015_030 = xr.open_rasterio('%s%s/%s_15-30cm_mean.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    var_030_060 = xr.open_rasterio('%s%s/%s_30-60cm_mean.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    var_060_100 = xr.open_rasterio('%s%s/%s_60-100cm_mean.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')

    var_000_005.values[var_000_005.values<0]=np.nan
    var_005_015.values[var_005_015.values<0]=np.nan
    var_015_030.values[var_015_030.values<0]=np.nan
    var_030_060.values[var_030_060.values<0]=np.nan
    var_060_100.values[var_060_100.values<0]=np.nan

    if sgvar == 'ocd':
        var_000_005.values*=(10.**-4)
        var_005_015.values*=(10.**-4)
        var_015_030.values*=(10.**-4)
        var_030_060.values*=(10.**-4)
        var_060_100.values*=(10.**-4)

    else:
        var_000_005.values*=(0.1)
        var_005_015.values*=(0.1)
        var_015_030.values*=(0.1)
        var_030_060.values*=(0.1)
        var_060_100.values*=(0.1)

    # regrid aggregated LCM
    if ii == 0:
        dy_orig = var_000_005.y.values[1]-var_000_005.y.values[0]
        dx_orig = var_000_005.x.values[1]-var_000_005.x.values[0]

        os.system('gdalwarp -te %f %f %f %f -r mode -s_srs EPSG:27700 -t_srs EPSG:4326 \
                    -tr %f %f -dstnodata 0 -overwrite %s %s' %
                    (W, S, E, N, dx_orig, dy_orig, lc_agg_file_25m,lc_agg_file_soilgrids))

        # load in the new raster
        lc_ref = xr.open_rasterio(lc_agg_file_soilgrids)[0]
        lc_id = [1,2,3,4,5,6,10,'_']

        # Reset grid
        Yorig = var_000_005.y.values
        Xorig = var_000_005.x.values

        # define destination lat / lon arrays
        Ydest = np.arange(N-dy_target/2.,S,-dy_target)
        Xdest = np.arange(W+dx_target/2.,E,dx_target)

        # calculate grid cell area
        areas = gst.calculate_cell_areas(Yorig,Xorig,projected=False)

    coords={'longitude':Xdest,'latitude':Ydest}
    for ii,lc in enumerate(landcover):
        print("Regridding %s for %s...             " % (sgvar,lc))
        path2dest_sub = '%s%s/' % (path2dest, lc)
        # clip variable grid using lat_mask and lon_mask
        if lc == 'non_forest':
            mask = lc_ref.values>2.5
        else:
            mask = lc_ref.values==lc_id[ii]
        print('\t0-30cm')
        # weights average
        var = (var_000_005.values*0.05 + var_005_015.values*0.10 +
                var_015_030.values*0.15)
        var_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest,
                                        areas, var, mask=mask, aggregation_mode='mean')
        # write to file
        if sgvar == 'ocd':
            nc_out = '%sSOC_000-030cm_mean_%s_5km.nc' % (path2dest_sub,lc)
            attrs={'details':'regridded SOC for %s for soil depths of 0-30cm, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids Organic Carbon Density 250 m data' % lc,
                    'name':'SOC stock for %s' % lc,
                    'units':'g/m2'}
            data_vars = {'SOC':(['latitude','longitude'],var_regrid*1000.,attrs)}
            var_out = xr.Dataset(data_vars=data_vars,coords=coords)
        elif sgvar in ['sand','silt','clay']:
            nc_out = '%s%s_000-030cm_mean_%s_5km.nc' % (path2dest_sub,sgvar,lc)
            attrs={'details':'Percentage %s for %s for soil depths of 0-30cm, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids 250 m data' % (sgvar, lc),
                    'name':'Percentage %s for %s' % (sgvar, lc),
                    'units':'%'}
            data_vars = {'%s' % sgvar : (['latitude','longitude'],var_regrid,attrs)}
            var_out = xr.Dataset(data_vars=data_vars, coords=coords)

        var_out.to_netcdf(path=nc_out)

        print('\t30-100cm')
        # weights average
        var = (var_030_060.values*0.30 + var_060_100.values*0.40)
        var_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest,
                                                areas, var, mask=mask,
                                                aggregation_mode='mean')
        # write to file
        if sgvar == 'ocd':
            attrs={'details':'regridded SOC for %s for soil depths of 30-100cm, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids Organic Carbon Density 250 m data' % lc,
                    'name':'SOC stock for %s' % lc,
                    'units':'g/m2'}
            nc_out = '%sSOC_030-100cm_mean_%s_5km.nc' % (path2dest_sub,lc)
            data_vars = {'SOC':(['latitude','longitude'],var_regrid*1000.,attrs)}
            var_out = xr.Dataset(data_vars=data_vars,coords=coords)

        elif sgvar in ['sand','silt','clay']:
            nc_out = '%s%s_030-100cm_mean_%s_5km.nc' % (path2dest_sub,sgvar,lc)
            attrs={'details':'Percentage %s for %s for soil depths of 30-100cm, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids 250 m data' % (sgvar, lc),
                    'name':'Percentage %s for %s' % (sgvar, lc),
                    'units':'%'}
            data_vars = {'%s' % sgvar : (['latitude','longitude'],var_regrid,attrs)}
            var_out = xr.Dataset(data_vars=data_vars, coords=coords)

        var_out.to_netcdf(path=nc_out)

# Now do uncertainty in ocd
for sgvar in ['ocd']:

    Q0_05_000_005 = xr.open_rasterio('%s%s/%s_0-5cm_Q0.05.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_05_005_015 = xr.open_rasterio('%s%s/%s_5-15cm_Q0.05.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_05_015_030 = xr.open_rasterio('%s%s/%s_15-30cm_Q0.05.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_05_030_060 = xr.open_rasterio('%s%s/%s_30-60cm_Q0.05.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_05_060_100 = xr.open_rasterio('%s%s/%s_60-100cm_Q0.05.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')

    Q0_05_000_005.values[Q0_05_000_005.values<0]=np.nan
    Q0_05_005_015.values[Q0_05_005_015.values<0]=np.nan
    Q0_05_015_030.values[Q0_05_015_030.values<0]=np.nan
    Q0_05_030_060.values[Q0_05_030_060.values<0]=np.nan
    Q0_05_060_100.values[Q0_05_060_100.values<0]=np.nan

    Q0_95_000_005 = xr.open_rasterio('%s%s/%s_0-5cm_Q0.95.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_95_005_015 = xr.open_rasterio('%s%s/%s_5-15cm_Q0.95.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_95_015_030 = xr.open_rasterio('%s%s/%s_15-30cm_Q0.95.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_95_030_060 = xr.open_rasterio('%s%s/%s_30-60cm_Q0.95.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')
    Q0_95_060_100 = xr.open_rasterio('%s%s/%s_60-100cm_Q0.95.tif' % (path2wgs84,sgvar,sgvar))[0].astype('float')

    Q0_95_000_005.values[Q0_95_000_005.values<0]=np.nan
    Q0_95_005_015.values[Q0_95_005_015.values<0]=np.nan
    Q0_95_015_030.values[Q0_95_015_030.values<0]=np.nan
    Q0_95_030_060.values[Q0_95_030_060.values<0]=np.nan
    Q0_95_060_100.values[Q0_95_060_100.values<0]=np.nan

    for ii,lc in enumerate(landcover):
        print("Regridding uncertainty in %s for %s...             " % (sgvar,lc))
        path2dest_sub = '%s%s/' % (path2dest, lc)
        if lc == 'non_forest':
            mask = lc_ref.values>2.5
        else:
            mask = lc_ref.values==lc_id[ii]

        unc = ((Q0_95_000_005.values*0.05 + Q0_95_005_015.values*0.15 +
                Q0_95_015_030.values*0.15) -
                (Q0_05_000_005.values*0.05 + Q0_05_005_015.values*0.15 +
                Q0_05_015_030.values*0.15))/2.

        if sgvar == 'ocd':
            unc*=10.**-4

        unc_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest,
                                                areas, unc, mask=mask, aggregation_mode='mean')

        nc_out = '%sSOC_000-030cm_uncertainty_%s_5km.nc' % (path2dest_sub,lc)
        attrs={'details':'half-width of the 5th-95th interquantile range, regridded OCD uncertatinty for %s, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids 250 m data' % lc,
                'name':'Uncertainty in SOC %s' % lc,
                'units':'g/m2'}
        data_vars = {'SOC_uncertainty' : (['latitude','longitude'], unc_regrid*1000., attrs)}
        unc_out = xr.Dataset(data_vars=data_vars, coords=coords)
        unc_out.to_netcdf(path=nc_out)

        unc = ((Q0_95_030_060.values*0.3 + Q0_95_060_100.values*0.4) -
                (Q0_05_030_060.values*0.3 + Q0_05_060_100.values*0.4))/2.

        if sgvar == 'ocd':
            unc*=10.**-4

        unc_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest,
                                                areas, unc, mask=mask, aggregation_mode='mean')

        nc_out = '%sSOC_030-100cm_uncertainty_%s_5km.nc' % (path2dest_sub,lc)
        attrs={'details':'half-width of the 5th-95th interquantile range, regridded OCD uncertatinty for %s, based on SoilGrids2 and CEH LCM2015, aggregated from SoilGrids 250 m data' % lc,
                'name':'Uncertainty in SOC %s' % lc,
                'units':'g/m2'}
        data_vars = {'SOC_uncertainty' : (['latitude','longitude'], unc_regrid*1000., attrs)}
        unc_out = xr.Dataset(data_vars=data_vars, coords=coords)
        unc_out.to_netcdf(path=nc_out)

"""
Additional loops to download global datasets at different resolutions
"""
"""
resolutions = [0.5,0.125,0.05]
for res in resolutions:
    os.system('mkdir %s/global_%f_degree/' % (path2wgs84,res))
    for sgvar in soilgrids_vars:
        os.system('mkdir %s/global_%f_degree/%s/' % (path2wgs84,res,sgvar))
        for layer in ['0-5','5-15','15-30','30-60','60-100','100-200']:
            print('resolution: %f degree; variable: %s; soil layer %scm            ' % (res,sgvar,layer),end='\r')
            ds = gdal.Warp('%s/global_%f_degree/%s/%s_%scm_mean.tif' % (path2wgs84,res,sgvar,sgvar,layer),
                            '%s/%s/%s_%scm_mean.vrt' % (sg_url,sgvar,sgvar,layer),
                            xRes=res, yRes=-res, dstSRS='EPSG:4326', srcSRS=igh,
                            outputBounds = (-180,-90,180,90),
                            outputBoundsSRS='EPSG:4326', resampleAlg='average')
            del ds

            if sgvar=='ocd':
                for quantile in ['Q0.05','Q0.95']:
                    ds = gdal.Warp('%s/global_%f_degree/%s/%s_%scm_%s.tif' % (path2wgs84,res,sgvar,sgvar,layer,quantile),
                                    '%s/%s/%s_%scm_%s.vrt' % (sg_url,sgvar,sgvar,layer,quantile),
                                    xRes=res, yRes=-res, dstSRS='EPSG:4326', srcSRS=igh,
                                    outputBounds = (-180,-90,180,90),
                                    outputBoundsSRS='EPSG:4326', resampleAlg='average')
                    del ds
"""
