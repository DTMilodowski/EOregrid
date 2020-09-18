"""
================================================================================
RegridLAI.py
--------------------------------------------------------------------------------
This set of functions contains the architecture of the regridding script whereby
area weighted averages are estimated for regridding fine scale resolution raster
datasets to coarse grids. There is scope to regrid both geographic (lat/lon) and
projected grids using the relevant flag.

This code was produced as part of the DARE-UK project, based on earlier scripts
produced by Jean-Francois Exbrayat.
--------------------------------------------------------------------------------
Author: David T. Milodowski
Date: 28/04/2020
================================================================================
"""
import numpy as np
import xarray as xr
import glob as glob
import geospatial_tools as gst


"""
--------------------------------------------------------------------------------
regridLAI
--------------------------------------------------------------------------------
    This function regrids a netcdf file onto a regular grid with same
    extent and resolution provided in the arguments.
    - path2orig : path to original file holding the data
    - path2dest : path to new file with regridded info
    - latres    : latitude resolution of target grid
    - lonres    : longitude resolution of target grid
    - variables : which variables to regrid
    - projected : boolean flag; True = projected (regular) grid;
                                False (default) = lat/lon grid
    - mask      : to only consider pixels according to a specified mask, for
                  example a land cover class. This must be a boolean array with
                  same shape as variables in source data
    - extent    : list with bounding extent, provided in form: [N,S,E,W]. If
                  unspecified, then apply original extent
    - coords_cell_centre : boolean flag to specify whether coordinates
                           correspond to cell centre
--------------------------------------------------------------------------------
"""
def regridLAI(path2orig, path2dest, Xres, Yres, mask=None, extent=None,
            variables=['LAI','LAI_ERR'], aggregation_mode = ['mean','quadrature'],
            projected=False, subset_label='?'):
    # First do a quick check to ensure we don't overwrite an existing file
    if len(glob.glob(path2dest)) > 0:
        print('Destination file "%s" already exists, please remove before proceeding' % (glob.glob(path2dest)[0]))
        return 1
    else:
        print("Regridding file " + path2orig.split('/')[-1])

    # Open the reference dataset
    ref = xr.open_dataset(path2orig)
    dim_names = gst.check_dim_names(ref)
    try:
        Xvar,Yvar=dim_names
    except:
        return 2

    # get resolution
    Yorig_res = np.abs(ref[Yvar].values[1]-ref[Yvar].values[0])
    Xorig_res = np.abs(ref[Xvar].values[1]-ref[Xvar].values[0])

    if extent is not None:
        N,S,E,W = extent[:]
        if 'lat' in dim_names and 'lon' in dim_names:
            if ref[Yvar].values[1]<ref[Yvar].values[0]:
                ref=ref.sel(lon=slice(W,E),lat=slice(N,S))
            else:
                ref=ref.sel(lon=slice(W,E),lat=slice(S,N))
        elif 'y' in dim_names and 'x' in dim_names:
            if ref[Yvar].values[1]<ref[Yvar].values[0]:
                ref=ref.sel(x=slice(W,E),y=slice(N,S))
            else:
                ref=ref.sel(x=slice(W,E),y=slice(S,N))
        else:
            print('Dimension names not currently accepted', dim_names)
            return 3

    # Reset grid
    Yorig = ref.coords[Yvar].values
    Xorig = ref.coords[Xvar].values

    # define scanning window size
    Ysize = np.abs(np.round(Yres/Yorig_res).astype('i'))
    Xsize = np.abs(np.round(Xres/Xorig_res).astype('i'))

    # define destination lat / lon arrays
    Ydest = np.arange(N-Yres/2.,S,-Yres)
    Xdest = np.arange(W+Xres/2.,E,Xres)

    # calculate grid cell area
    areas = gst.calculate_cell_areas(Yorig,Xorig,projected=projected)

    # regrid to new resolution
    regridded = {}
    for iv,varname in enumerate(variables):
        print("Regridding variable %s ... " % (varname))
        # clip variable grid using lat_mask and lon_mask
        variable = ref.variables[varname].values
        var_regrid,fraction=gst.regrid_single(Yorig, Xorig, Ydest, Xdest, Ysize, Xsize,
                                        areas, variable, mask=mask,
                                        aggregation_mode=aggregation_mode[iv])
        regridded[varname]=var_regrid.copy()
        if iv==0:
            regridded['fraction']=fraction.copy()

    # compile new netcdf file
    if projected:
        coords = {Yvar: ([Yvar],Ydest,{'units':'m','long_name':'y coordinate'}),
                Xvar: ([Xvar],Xdest,{'units':'m','long_name':'x coordinate'})}
    else:
        coords = {Yvar: ([Yvar],Ydest,{'units':'degrees_north','long_name':'latitude'}),
                Xvar: ([Xvar],Xdest,{'units':'degrees_east','long_name':'longitude'})}

    attrs = ref.attrs.copy()
    attrs['title'] += '; tiled at %.3f for %s subset' % (Xres,subset_label)
    attrs['history'] += '; subsetted and regridded'

    data_vars = {}
    for iv,varname in enumerate(list(regridded.keys())):
        try:
            var_attrs = ref[varname].attrs.copy()
            var_attrs['long_name'] += '; tiled at %.3f for %s subset' % (Xres,subset_label)
            var_attrs['standard_name'] += ' for %s subset' % (Xres,subset_label)
            data_vars[varname] = (['lat','lon'],regridded[varname],var_attrs.copy())
        except:
            if varname=='fraction':
                var_attrs = {}
                var_attrs['long_name'] = 'Fraction of grid cell occupied by %s' % (subset_label)
                var_attrs['standard_name'] = 'fraction %s' % (subset_label)
                var_attrs['grid_mapping'] = 'crs'
                var_attrs['units'] = ''
                data_vars[varname] = (['lat','lon'],regridded[varname],var_attrs.copy())
            else:
                print('variable not found')

    regrid_ds = xr.Dataset(data_vars=data_vars,coords=coords)
    regrid_ds.to_netcdf(path=path2dest)
    return regrid_ds
