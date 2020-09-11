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
import glob
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
def regridLAI(path2orig,path2dest,Xres,Yres,variables=['LAI','LAI_ERR'],
                projected=False,mask=None, extent=None, coords_cell_centre = False):
    # First do a quick check to ensure we don't overwrite an existing file
    if len(glob.glob(path2dest)) > 0:
        print('Destination file "%s" already exists, please remove before proceeding' % (glob.glob(path2dest)[0]))
        return 1
    else:
        print("Regridding file " + path2orig.split('/')[-1])

    # Open the reference dataset
    ref = xr.open_dataset(path2orig)
    dim_names = gst.check_dim_names(ref)
    if len(dim_names)<2:
        return 2
    else:
        Xvar,Yvar=dim_names

    # Set up grids
    Yorig = ref.coords[Yvar].values
    Xorig = ref.coords[Xvar].values
    temporal = False
    if 'time' in ref.dims:
        temporal=True

    # get resolution
    Yorig_res = np.abs(Yorig[1]-Yorig[0])
    Xorig_res = np.abs(Xorig[1]-Xorig[0])

    # adjust the lat/lon to centre of cell
    if not coords_cell_centre:
        Xorig = Xorig+Xorig_res/2.
        Yorig = Yorig-Yorig_res/2.

    # Apply mask to bbox extent if required
    if extent is not None:
        N,S,E,W = extent
        Ymask = np.all((Yorig<=N+Yres,Yorig>=S-Yres),axis=0); Yorig = Yorig[Ymask]
	    Xmask = np.all((Xorig<=E+Xres,Xorig>=W-Xres),axis=0); Xorig = Xorig[Xmask]


    # define scanning window size
    Ysize = np.abs(np.round(Yres/Yorig_res).astype('i'))
    Xsize = np.abs(np.round(Xres/Xorig_res).astype('i'))

    # define destination lat / lon arrays
    Ydest = np.arange(N-Yres/2.,S,-Yres)
	Xdest = np.arange(W+Xres/2.,E,Xres)

    # calculate grid cell area
    areas = gst.calculate_cell_areas(Yorig,Xorig,projected=projected)

    # regrid to new resolution
    for iv,varname in enumerate(variables):
        print("Regridding variable %s ... " % (varname))
        # clip variable grid using lat_mask and lon_mask
        if temporal:
        	variable = ref.variables[varname].values[:,Ymask,:]
	    	variable = variable[:,:,Xmask]
            target,fraction=gst.regrid_single(Yorig,Xorig,Ydest,Xdest,Ysize,Xsize,
            areas,variable,mask=mask,temporal=temporal)
	    else:
            variable = ref.variables[varname][Ymask,:]
	    	variable = variable[:,Xmask]
            target,fraction=gst.regrid_single(Yorig,Xorig,Ydest,Xdest,Ysize,Xsize,
            areas,variable,mask=mask,temporal=temporal)
