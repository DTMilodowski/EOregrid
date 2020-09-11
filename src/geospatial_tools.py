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
    areas = np.zeros([Y.size,X.size])
	print(areas.shape)
    if projected:
        areas+=np.abs(Xres*Yres)
    else:
        for iy,y in enumerate(Y):
            areas[iy]= (6371e3)**2 * (np.deg2rad(0+Xres/2.)-np.deg2rad(0-Xres/2.)) * \
                        (np.sin(np.deg2rad(y+Yres/2.))-np.sin(np.deg2rad(y-Yres/2.)))
    return areas


def regrid_single(Yorig,Xorig,Ydest,Xdest,Ysize,Xsize,
            areas,variable,mask=None,temporal=False):
    # Set up arrays
    if temporal = True:
        target = np.zeros([1,Ydest.size,Xdest.size])*np.nan
    else:
        target = np.zeros([Ydest.size,Xdest.size])*np.nan
    fraction = np.zeros(target.shape)*np.nan
    counter = 0

    # regridding, here we go!
    for iY, y in enumerate(Ydest):
        for iX, x in enumerate(Xdest):
            counter+=1
            print '\rRegridding pixel %i / %i' % (counter, len(Ydest)*len(Xdest))
		    slcarea = areas[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]

            #if there's a mask, extract the data
            if mask is not None:
                slcmask = mask[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
            else:
                slcmask = np.ones([Ysize,Xsize],dtype='bool')

            #check that there's data within mask at this pixel
            if slcmask.sum() != 0:
                # extract the high res pixels inside the destinatino pixel
                if temporal:
                    slcdata = variable[0,(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]
                else:
                    slcdata = variable[(iY*Ysize):((iY+1)*Ysize),(iX*Xsize):((iX+1)*Xsize)]

                # if there is a built in mask to the dataset, update (this might be redundant in new version)
                if 'mask' in dir(slcdata):
                    #replace the mask where the land cover is not taken into account
			        print(variable.shape,slcdata.shape)
                    slcdata.mask[~slcmask] = True
                    if slcdata.mask.sum() != slcdata.size:
				        print(iY,iX,target.shape,slcdata.shape,slcarea.shape)
                    if temporal:
                        target[0,iY,iX] = (slcdata*slcarea).sum()/(~slcdata.mask*slcarea).sum()
				    else:
                        target[iY,iX] = (slcdata*slcarea).sum()/(~slcdata.mask*slcarea).sum()

                    # first pass, extract the fraction of pixel with valid data
                    if temporal:
                        fraction[0,iY,iX] = (~slcdata.mask*slcarea).sum()/(slcarea.sum())
                    else:
                        fraction[iY,iX] = (~slcdata.mask*slcarea).sum()/(slcarea.sum())
                else:
                    if temporal:
                        target[0,iY,iX] = (slcdata*slcarea).sum()/slcarea.sum()
                        fraction[0,iY,iX]=1.
    			    else:
                        target[iY,iX] = (slcdata*slcarea).sum()/slcarea.sum()
                        fraction[iY,iX]=1.

    return target,fraction

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
