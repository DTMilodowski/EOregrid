import glob
import numpy as np
import xarray as xr
import pandas as pd

# quick function to find nearest row and column of the grid
def find_nearest_row_col(Xloc,Yloc,Xref,Yref):
    row = np.argmin(np.abs(Yref-Yloc))
    col = np.argmin(np.abs(Xref-Xloc))
    return row,col

# define path to lai data
path2lai = '/disk/scratch/local.2/copernicus/LAI_300m/'

# define a path to the output
path2output = '~/' #'/path/to/output/'

# list of lai files to be processed
lai_files = glob.glob('%s*PROBAV*.nc' % path2lai)
lai_files.sort()
n_lai=len(lai_files)
query_variables = ['LAI','RMSE','QFLAG']

# read in csv file with site locations
sites_file = '~/TG_site_coords.csv'
sites = pd.read_csv(sites_file)
n_sites = sites.shape[0]

# get the grid coordinates for the LAI data
ds = xr.open_dataset(lai_files[0])
lon_ref = ds.lon.values
lat_ref = ds.lat.values
ds=None
# now loop through the sites and extract all the data
for ii in range(0,n_sites):
    lai_info = {}
    for var in query_variables:
        lai_info[var] = np.zeros(n_lai)*np.nan
    lai_info['date'] = np.zeros(n_lai,dtype='U8')
    site = sites.iloc[ii]
    print('processing lai time series for %s' % site.Site_code)
    row,col = find_nearest_row_col(site.Lon,site.Lat,lon_ref,lat_ref)
    # loop through the datasets and extract the desired information for the 
    # nearest pixel to the target location
    for jj,lai_file in enumerate(lai_files):
        print('\tprocessing... %.1f percent' % (float(jj)/float(n_lai)*100), end='\r')
        ds = xr.open_dataset(lai_file)
        for var in query_variables:
            lai_info[var][jj] = float(ds[var][row,col].values)
        lai_info['date'][jj] = lai_file.split('/')[-1][10:18]
        ds=None
    print('\tprocessing... done!                 ')
    # combine into pandas dataframe and write to a csv file
    df = pd.DataFrame(lai_info)
    outfile = '%s/copernicus_lai_300m_%s_nearest.csv' % (path2output,site.Site_code)
    df.to_csv(outfile)
