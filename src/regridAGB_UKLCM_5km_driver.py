"""
regridAGB_UKLCM_5km_driver.py
--------------------------------------------------------------------------------
Regridding script to regrid ESA-BIOMASS CCI to 5km, stratified by aggregated
land cover indicated by the CEH LCM dataset.
"""
import os
import numpy as np
import xarray as xr

import geospatial_tools as gst


"""
PART A: LOAD IN THE REQUIRED DATA
This program will use
- ESA BIOMASS-CCI AGB AND SD
- CEH 2015 UK land cover map
"""
dx_target = 0.05
dy_target = 0.05
N = 61.25
S = 49.75
E = 2.25
W = -10.5
# Set up directory structure as required
path2agb = '/home/dmilodow/DataStore_DTM/DARE_UK/EOregrid/data/ESACCI_BIOMASS_vrt/'
path2dest =  '/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/UKLCM_5km/'

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
lc_agg_file_25m = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-25m_3627507/lcm2015gb25m_aggregated.tif'
lc_agg_file_BIOMASS = '/home/dmilodow/DataStore_DTM/DARE_UK/Data/lcm-2015-300m-wgs84/lcm2015gb_ESACCI_BIOMASS_aggregated.tif'
lc_25m_agg = xr.open_rasterio(lc_agg_file_25m).sel(band=1)

# open AGB file

years = [2010,2017,2018]
for year in years:
    agb = xr.open_rasterio('%s/%i/ESACCI-BIOMASS-L4-AGB-MERGED-100m-%i-fv2.0.vrt' % (path2agb,year,year))[0]
    agb_sd = xr.open_rasterio('%s/%i/ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-%i-fv2.0.vrt' % (path2agb,year,year))[0]

    # regrid aggregated LCM
    dy_orig = agb.y.values[1]-agb.y.values[0]
    dx_orig = agb.x.values[1]-agb.x.values[0]

    W_ = W-dx_orig/2.; E_ = E+dx_orig/2.
    N_ = N-dy_orig/2.; S_ = S+dy_orig/2.
    if year == 2010:
        os.system('gdalwarp -r mode -s_srs EPSG:27700 -t_srs EPSG:4326 \
                -te %f %f %f %f -tr %f %f -dstnodata 0 -overwrite %s %s' %
                (W_, S_, E_, N_, dx_orig, dy_orig, lc_agg_file_25m,lc_agg_file_BIOMASS))

    # load in the new raster
    lc_ref = xr.open_rasterio(lc_agg_file_BIOMASS)[0]
    lc_id = [1,2,3,4,5,6,10,'_']


    agb_ref=agb.sel(x=slice(W,E),y=slice(N,S))
    agb_sd_ref=agb_sd.sel(x=slice(W,E),y=slice(N,S))

    # Reset grid
    Yorig = agb_ref.y.values
    Xorig = agb_ref.x.values

    # define destination lat / lon arrays
    Ydest = np.arange(N-dy_target/2.,S,-dy_target)
    Xdest = np.arange(W+dx_target/2.,E,dx_target)

    # calculate grid cell area
    areas = gst.calculate_cell_areas(Yorig,Xorig,projected=False)

    coords = {'latitude': (['latitude'],Ydest,{'units':'degrees_north','long_name':'latitude'}),
            'longitude': (['longitude'],Xdest,{'units':'degrees_east','long_name':'longitude'})}

    for ii,lc in enumerate(landcover):
        print("Regridding variable AGB for %s...             " % lc)
        # clip variable grid using lat_mask and lon_mask
        if lc == 'non_forest':
            mask = lc_ref.values>2.5
        else:
            mask = lc_ref.values==lc_id[ii]
        agb_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest,
                                            areas, agb_ref.values, mask=mask,
                                            aggregation_mode='mean')
        # repeat for sd
        sd_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest, 
                                            areas, agb_sd_ref.values, mask=mask,
                                            aggregation_mode='mean')

        # set up attributes
        attrs = agb.attrs.copy()
        attrs['title'] = 'ESACCI Aboveground Biomass data for %i; tiled at %.3f for %s subset' % (year,dx_target,lc)
        attrs['history'] = 'ESACCI Aboveground Biomass data for %i; subsetted for %s based on the UKLCM2015 land cover map, and regridded to %.3f degrees for ingestion into stratified CARDAMOM simulations for the UK C cycle' % (year,lc,dx_target)

        data_vars = {}
        agb_attrs = {'long_name' : 'Aboveground Biomass (Mg / ha) for %s in %i; tiled at %.3f' % (lc,year,dx_target),
                    'standard_name' : 'AGBiomass', 'units' : 'Mg ha-1'}
        data_vars['AGBiomass'] = (['latitude','longitude'],agb_regrid,agb_attrs.copy())
        unc_attrs = {'long_name' : 'Uncertainty in Aboveground Biomass (Mg / ha) for %s in %i; tiled at %.3f' % (lc,year,dx_target),
                    'standard_name' : 'AGBiomass_Uncertainty', 'units' : 'Mg ha-1'}
        data_vars['AGBiomass_Uncertainty'] = (['latitude','longitude'],sd_regrid,unc_attrs.copy())

        # write to file
        agb_ds = xr.Dataset(data_vars=data_vars,coords=coords)
        path2dest_sub = '%s%s/' % (path2dest, lc)
        agb_ds.to_netcdf(path='%sAGBiomass_stocks_%i_with_lat_lon_%s.nc' % (path2dest_sub,year,lc))
