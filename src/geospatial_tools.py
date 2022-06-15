import numpy as np

"""
--------------------------------------------------------------------------------
check_dim_names
--------------------------------------------------------------------------------
"""
def check_dim_names(ref):

    # check possible name of x
    if 'x' in ref.dims:
        Xvar = 'x'
    elif 'X' in ref.dims:
        Xvar = 'X'
    elif 'lon' in ref.dims:
        Xvar = 'lon'
    elif 'Lon' in ref.dims:
        Xvar = 'Lon'
    elif 'longitude' in ref.dims:
        Xvar = 'longitude'
    elif 'Longitude' in ref.dims:
        Xvar = 'Longitude'
    else:
        print('Cannot find x variable in coords')
        return 2

    # check possible name of y
    if 'y' in ref.dims:
        Yvar = 'y'
    elif 'Y' in ref.dims:
        Yvar = 'Y'
    elif 'lat' in ref.dims:
        Yvar = 'lat'
    elif 'Lat' in ref.dims:
        Yvar = 'Lat'
    elif 'latitude' in ref.dims:
        Yvar = 'latitude'
    elif 'Latitude' in ref.dims:
        Yvar = 'Latitude'
    else:
        print('Cannot find y variable in coords')
        return 2

    return Xvar,Yvar


"""
--------------------------------------------------------------------------------
calculate_cell_areas
--------------------------------------------------------------------------------
"""
def calculate_cell_areas(Y,X,projected=False):
    #calculate grid cell area
    Yres = Y[1]-Y[0]; Xres = X[1]-X[0]
    areas = np.zeros((Y.size,X.size))
    print(areas.shape)
    if projected:
        areas+=np.abs(Xres*Yres)
    else:
        for iy,y in enumerate(Y):
            areas[iy]= (6371e3)**2 * (np.deg2rad(0+Xres/2.)-np.deg2rad(0-Xres/2.)) * \
                        (np.sin(np.deg2rad(y+Yres/2.))-np.sin(np.deg2rad(y-Yres/2.)))
    return areas


def regrid_single(Yorig,Xorig,Ydest,Xdest,
            areas,variable,mask=None,aggregation_mode='mean'):
    # Set up array
    dy_target = Ydest[1]-Ydest[0]
    dx_target = Xdest[1]-Xdest[0]
    target = np.zeros([Ydest.size,Xdest.size])*np.nan
    fraction = np.zeros(target.shape)*np.nan
    counter = 0

    # regridding, here we go!
    for iY, y in enumerate(Ydest):
        slcYsub = abs(Yorig-y)<abs(dy_target/2.)
        for iX, x in enumerate(Xdest):
            slcXsub = abs(Xorig-x)<abs(dx_target/2.)
            slcsub = np.ix_(slcYsub,slcXsub)
            counter+=1
            if counter%100 == 0:
                print('Regridding pixel %i / %i' % (counter, len(Ydest)*len(Xdest)),end='\r')

            #if there's a mask, extract the data
            if mask is not None:
                #slcmask = mask[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
                slcmask = mask[slcsub]
            else:
                #slcmask = np.ones([Ysize,Xsize],dtype='bool')
                slcmask = np.ones([slcYsub.sum(),slcXsub,sum()],dtype='bool')

            #check that there's data within mask at this pixel
            if slcmask.sum() != 0:
                # extract contributing pixel areas
                #slcarea = areas[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
                slcarea = areas[slcsub]
                # extract the high res pixels inside the destinatino pixel
                #slcdata = variable[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
                slcdata = variable[slcsub]

                fraction[iY,iX]=(slcmask*slcarea).sum()/slcarea.sum()

                # check whether there is any data
                slcmask_filter_nodata = np.isfinite(slcdata)*slcmask
                if slcmask_filter_nodata.sum() != 0:
                    if aggregation_mode == 'mean':
                        target[iY,iX] = np.nansum(slcmask_filter_nodata*slcdata*slcarea)/np.nansum(slcmask_filter_nodata*slcarea)
                    elif aggregation_mode == 'quadrature': # for uncertainties
                        target[iY,iX] = np.sqrt(np.nansum(((slcmask_filter_nodata*slcdata*slcarea)/np.nansum(slcmask_filter_nodata*slcarea))**2))
                    elif aggregation_mode == 'sum':
                        target[iY,iX] = np.nansum(slcmask*slcdata)
                    else: # assume mean if not quadrature
                        target[iY,iX] = np.nansum(slcmask*slcdata*slcarea)/np.nansum(slcmask*slcarea)

    return target,fraction

"""
extract nearest neighbour for a point source
"""
def find_nearest_row_col(Xloc,Yloc,Xref,Yref):
    row = np.argmin(np.abs(Yref-Yloc))
    col = np.argmin(np.abs(Xref-Xloc))
    return row,col

"""
def regrid_temporal(Yorig,Xorig,Ydest,Xdest,Ysize,Xsize,time,
            areas,variable,mask=None,mask_fraction=False):
    # Set up arrays
    target = np.zeros([time.size,Ydest.size,Xdest.size])*np.nan
    if mask_fraction:
        fraction = np.zeros(target.shape[1],target.shape[2])*np.nan
    counter = 0

    # regridding, here we go!
    for iT, t in time:
        for iY, y in enumerate(Ydest):
            for iX, x in enumerate(Xdest):
            counter+=1
            print '\rRegridding %.5f %' % (counter/float(time.size,Ydest.size*Xdest.size))
		    slcarea = areas[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]

            #if there's a mask, extract the data
            if mask is not None:
                slcmask = mask[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
            else:
                slcmask = np.ones([Ysize,Xsize],dtype='bool')

            #check that there's data within mask at this pixel
            if slcmask.sum() != 0:
                # extract the high res pixels inside the destinatino pixel
                slcdata = variable[iT,(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]

                # if there is a built in mask to the dataset, update (this might be redundant in new version)
                if 'mask' in dir(slcdata):
                    #replace the mask where the land cover is not taken into account
			        print(variable.shape,slcdata.shape)
                    slcdata.mask[~slcmask] = True
                    if slcdata.mask.sum() != slcdata.size:
				        print(iY,iX,target.shape,slcdata.shape,slcarea.shape)
                    target[iT,iY,iX] = (slcdata*slcarea).sum()/(~slcdata.mask*slcarea).sum()

                    # first pass, extract the fraction of pixel with valid data
                    if mask_fraction:
                        if iT==0:
                            fraction[iY,iX] = (~slcdata.mask*slcarea).sum()/(slcarea.sum())
                else:
                    target[iT,iY,iX] = (slcdata*slcarea).sum()/slcarea.sum()
                    if mask_fraction:
                        if iT==0:
                            fraction[iY,iX]=1.

    if mask_fraction:
        return target,fraction
    else:
        return target
"""
