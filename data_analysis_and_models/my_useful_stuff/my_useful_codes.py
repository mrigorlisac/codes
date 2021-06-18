def array_to_tif(array, x_pixels = 8, y_pixels = 20, PIXEL_SIZE = 0.25, x_min = -74, y_max = -37,
                 dst_filename = '/Users/lisac/Documents/Data_analysis/TRMM_test.tif',
                 wkt_projection = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'):
    """Transform ordinary array format to .tiff geopositioned format
    inputs: array to transform, size, latitude, longitude
    outputs: array with embedded geographic location - .tiff"""
    from osgeo import gdal
    driver = gdal.GetDriverByName('GTiff')

    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float32, )

    dataset.SetGeoTransform((
        x_min,       # 0
        PIXEL_SIZE,  # 1
        0,           # 2
        y_max,       # 3
        0,           # 4
        -PIXEL_SIZE))  

    dataset.SetProjection(wkt_projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.
    return dataset, dataset.GetRasterBand(1)  # If you need to return, remember to return also the dataset because 
                                              # the band don`t live without dataset



def append_nans_to_timeseries(dataaaray, start_date, end_date, sampling_frequency, mode='beginning'):
    
    '''Function to add Nans to the beginning or the end of timeseries because incomplete series mess up some statistics. 
    Output is the dataarray with the same dimensions as input, but with added nans.
    Inputs: start_date, end_date, sampling_frequency --> strings
    mode: append to beginning or the end of timeseries'''

    import numpy as np
    import pandas as pd
    
    time_coords=pd.date_range(start=start_date, end=end_date, freq=sampling_frequency)
    nan_array = np.zeros((time_coords.shape[0], dataaaray.latitude.shape[0], dataaaray.longitude.shape[0]))
    nan_array[:] = np.nan


    nan_da = xr.DataArray(nan_array, coords={'time': time_coords, 'latitude': dataaaray.latitude, 
                                             'longitude': dataaaray.longitude}, dims=['time', 'latitude', 'longitude'])

    if mode == 'beginning':
        precip_appended = xr.concat([nan_da, dataaaray], dim='time')
        
    elif mode == 'end':
        precip_appended = xr.concat([dataaaray, nan_da], dim='time')
        
    else:
        print('not valid mode')

    return precip_appended





def spatial_pdf_params(dataarray, season='full', threshold=0.1, mode='gamma_pdf'):
    
    '''Function used for fitting the PDFs and getting the PDF parameters in space, for each grid cell.
    Takes xarray 3-dim timeseries, as input.
    season = full, wet, dry, everything
    mode = gamma_pdf, exponential_pdf, genpareto_pdf
    Seasons are divided into wet and dry based on calendar months.
    Output is dataarray of PDF parameters'''
    
    from scipy import stats
    import numpy as np
    import xarray as xr
    import pandas as pd

    if mode == 'gamma_pdf':
    
        m = len(dataarray.latitude) # number of rows (latitude)
        n = len(dataarray.longitude) # number of columns (longitude)
        season = season

        # threshold for rainy days
        rainfall_data = dataarray.where(dataarray>threshold).copy()


        # wet and dry season masks
        month = rainfall_data.groupby('time.month').apply(lambda x: x).month
        wet = (month > 4) & (month <= 10) # April-October
        dry = (month <= 4) | (month > 10)


        if season == 'full':
            gamma_kappa_full_ar = np.zeros((m,n))
            gamma_betta_full_ar = np.zeros((m,n))

            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    kappa_full, loc, theta_full = stats.gamma.fit(rainfall_data.isel(latitude=i, longitude=j).dropna('time').values, floc=0)
                    betta_full=1/theta_full
                    gamma_kappa_full_ar[i,j] = kappa_full
                    gamma_betta_full_ar[i,j] = betta_full

            gamma_kappa_full_da = xr.DataArray(gamma_kappa_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_full_da = xr.DataArray(gamma_betta_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            return gamma_kappa_full_da, gamma_betta_full_da



        if season == 'wet':
            gamma_kappa_wet_ar = np.zeros((m,n))
            gamma_betta_wet_ar = np.zeros((m,n))
            wet_data = rainfall_data.where(wet)


            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    kappa_wet, loc, theta_wet = stats.gamma.fit(wet_data.isel(latitude=i, longitude=j).dropna('time').values, floc=0)
                    betta_wet=1/theta_wet
                    gamma_kappa_wet_ar[i,j] = kappa_wet
                    gamma_betta_wet_ar[i,j] = betta_wet

            gamma_kappa_wet_da = xr.DataArray(gamma_kappa_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_wet_da = xr.DataArray(gamma_betta_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            return gamma_kappa_wet_da, gamma_betta_wet_da       





        if season == 'dry':
            gamma_kappa_dry_ar = np.zeros((m,n))
            gamma_betta_dry_ar = np.zeros((m,n))
            dry_data = rainfall_data.where(dry)


            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    kappa_dry, loc, theta_dry = stats.gamma.fit(dry_data.isel(latitude=i, longitude=j).dropna('time').values, floc=0)
                    betta_dry=1/theta_dry
                    gamma_kappa_dry_ar[i,j] = kappa_dry
                    gamma_betta_dry_ar[i,j] = betta_dry

            gamma_kappa_dry_da = xr.DataArray(gamma_kappa_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_dry_da = xr.DataArray(gamma_betta_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            return gamma_kappa_dry_da, gamma_betta_dry_da



        if season == 'everything':
            gamma_kappa_wet_ar = np.zeros((m,n))
            gamma_betta_wet_ar = np.zeros((m,n))
            gamma_kappa_dry_ar = np.zeros((m,n))
            gamma_betta_dry_ar = np.zeros((m,n))
            gamma_kappa_full_ar = np.zeros((m,n))
            gamma_betta_full_ar = np.zeros((m,n))
            wet_data = rainfall_data.where(wet)
            dry_data = rainfall_data.where(dry)


            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    kappa_wet, loc, theta_wet = stats.gamma.fit(wet_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    betta_wet=1/theta_wet
                    kappa_dry, loc, theta_dry = stats.gamma.fit(dry_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    betta_dry=1/theta_dry
                    kappa_full, loc, theta_full = stats.gamma.fit(rainfall_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    betta_full=1/theta_full

                    gamma_kappa_wet_ar[i,j] = kappa_wet
                    gamma_betta_wet_ar[i,j] = betta_wet
                    gamma_kappa_dry_ar[i,j] = kappa_dry
                    gamma_betta_dry_ar[i,j] = betta_dry
                    gamma_kappa_full_ar[i,j] = kappa_full
                    gamma_betta_full_ar[i,j] = betta_full 


            gamma_kappa_full_da = xr.DataArray(gamma_kappa_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_full_da = xr.DataArray(gamma_betta_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            gamma_kappa_wet_da = xr.DataArray(gamma_kappa_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_wet_da = xr.DataArray(gamma_betta_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            gamma_kappa_dry_da = xr.DataArray(gamma_kappa_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            gamma_betta_dry_da = xr.DataArray(gamma_betta_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])


            gamma_pdf_full_da = xr.concat([gamma_kappa_full_da, gamma_betta_full_da], dim=pd.Index(['shape', 'rate'], name='gamma_pdf_param'))
            gamma_pdf_wet_da = xr.concat([gamma_kappa_wet_da, gamma_betta_wet_da], dim=pd.Index(['shape', 'rate'], name='gamma_pdf_param'))
            gamma_pdf_dry_da = xr.concat([gamma_kappa_dry_da, gamma_betta_dry_da], dim=pd.Index(['shape', 'rate'], name='gamma_pdf_param'))


            dim_name = 'season'
            coords_names = ['full', 'wet', 'dry']

            gamma_pdf_all_da = xr.concat([gamma_pdf_full_da, gamma_pdf_wet_da, gamma_pdf_dry_da], dim=pd.Index(coords_names, name=dim_name))

            return gamma_pdf_all_da
        
        
        
        
    if mode == 'exponential_pdf':
    
        m = len(dataarray.latitude) # number of rows (latitude)
        n = len(dataarray.longitude) # number of columns (longitude)
        season = season

        # threshold for rainy days
        rainfall_data = dataarray.where(dataarray>threshold).copy()


        # wet and dry season masks
        month = rainfall_data.groupby('time.month').apply(lambda x: x).month
        wet = (month > 4) & (month <= 10) # April-October
        dry = (month <= 4) | (month > 10)
        
        
        
        if season == 'everything':
            expon_alpha_full_ar = np.zeros((m,n))
            expon_alpha_wet_ar = np.zeros((m,n))
            expon_alpha_dry_ar = np.zeros((m,n))
            
            wet_data = rainfall_data.where(wet)
            dry_data = rainfall_data.where(dry)


            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    loc, alpha_wet = stats.expon.fit(wet_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    
                    loc, alpha_dry = stats.expon.fit(dry_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    
                    loc, alpha_full = stats.expon.fit(rainfall_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    

                    expon_alpha_full_ar[i,j] = alpha_full
                    expon_alpha_wet_ar[i,j] = alpha_wet
                    expon_alpha_dry_ar[i,j] = alpha_dry


            expon_alpha_full_da = xr.DataArray(expon_alpha_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            

            expon_alpha_wet_da = xr.DataArray(expon_alpha_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            

            expon_alpha_dry_da = xr.DataArray(expon_alpha_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            



            dim_name = 'season'
            coords_names = ['full', 'wet', 'dry']

            expon_pdf_all_da = xr.concat([expon_alpha_full_da, expon_alpha_wet_da, expon_alpha_dry_da], dim=pd.Index(coords_names, name=dim_name))
            expon_pdf_all_da.name = 'scale'
            
            
            return expon_pdf_all_da
        
        
        
        
    if mode == 'genpareto_pdf':
        m = len(dataarray.latitude) # number of rows (latitude)
        n = len(dataarray.longitude) # number of columns (longitude)
        season = season

        # threshold for rainy days
        rainfall_data = dataarray.where(dataarray>threshold).copy()


        # wet and dry season masks
        month = rainfall_data.groupby('time.month').apply(lambda x: x).month
        wet = (month > 4) & (month <= 10) # May-October
        dry = (month <= 4) | (month > 10)
        
        

        if season == 'everything':
            genpareto_shape_wet_ar = np.zeros((m,n))
            genpareto_scale_wet_ar = np.zeros((m,n))
            genpareto_shape_dry_ar = np.zeros((m,n))
            genpareto_scale_dry_ar = np.zeros((m,n))
            genpareto_shape_full_ar = np.zeros((m,n))
            genpareto_scale_full_ar = np.zeros((m,n))
            
            wet_data = rainfall_data.where(wet)
            dry_data = rainfall_data.where(dry)


            for i in range(len(rainfall_data.latitude)):
                for j in range(len(rainfall_data.longitude)):
                    shape_wet, loc, scale_wet = stats.genpareto.fit(wet_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    
                    shape_dry, loc, scale_dry = stats.genpareto.fit(dry_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    
                    shape_full, loc, scale_full = stats.genpareto.fit(rainfall_data.isel(latitude=i, longitude=j).dropna(dim='time').values, floc=0)
                    

                    genpareto_shape_wet_ar[i,j] = shape_wet
                    genpareto_scale_wet_ar[i,j] = scale_wet
                    genpareto_shape_dry_ar[i,j] = shape_dry
                    genpareto_scale_dry_ar[i,j] = scale_dry
                    genpareto_shape_full_ar[i,j] = shape_full
                    genpareto_scale_full_ar[i,j] = scale_full 


            genpareto_shape_full_da = xr.DataArray(genpareto_shape_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            genpareto_scale_full_da = xr.DataArray(genpareto_scale_full_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            genpareto_shape_wet_da = xr.DataArray(genpareto_shape_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            genpareto_scale_wet_da = xr.DataArray(genpareto_scale_wet_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])

            genpareto_shape_dry_da = xr.DataArray(genpareto_shape_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])
            genpareto_scale_dry_da = xr.DataArray(genpareto_scale_dry_ar, coords=[rainfall_data.coords['latitude'], rainfall_data.coords['longitude']], 
                 dims=['latitude', 'longitude'])


            genpareto_pdf_full_da = xr.concat([genpareto_shape_full_da, genpareto_scale_full_da], dim=pd.Index(['shape', 'scale'], name='genpareto_pdf_param'))
            genpareto_pdf_wet_da = xr.concat([genpareto_shape_wet_da, genpareto_scale_wet_da], dim=pd.Index(['shape', 'scale'], name='genpareto_pdf_param'))
            genpareto_pdf_dry_da = xr.concat([genpareto_shape_dry_da, genpareto_scale_dry_da], dim=pd.Index(['shape', 'scale'], name='genpareto_pdf_param'))


            dim_name = 'season'
            coords_names = ['full', 'wet', 'dry']

            genpareto_pdf_all_da = xr.concat([genpareto_pdf_full_da, genpareto_pdf_wet_da, genpareto_pdf_dry_da], dim=pd.Index(coords_names, name=dim_name))
            
            
            return genpareto_pdf_all_da
        
        

def my_to_datetime(date_str):
    ''' Function to translate strings into datetime format with special care of "24" and 
    not "00" for the day start/end'''
    if date_str[8:10] != '24':
        return pd.to_datetime(date_str, format='%Y%m%d%H')

    date_str = date_str[0:8] + '00' + date_str[10:]
    return pd.to_datetime(date_str, format='%Y%m%d%H') + dt.timedelta(days=1)
        



def average_stations_to_grid_cells(station_data, latbins=None, lonbins=None):
    '''Function used to average stations that fall in a particular grid cell, in order to compare them to gridded rainfall datasets.
    Input is: - processed station data with calculated shape parameter of some PDF.
              - bounding box in latitude and longitude to which to translate the station data.
              - e.g. latbins = np.r_[21.5:25.55:0.1]
                     lonbins=numpy.r_[120:122.1:0.1]
              
              <xarray.Dataset>
Dimensions:                (station: 467, time: 42609)
Coordinates:
  * time                   (time) datetime64[ns] 2003-01-01 ... 2017-08-01
  * station                (station) object 'C0A510' 'C0A520' ... 'C1Z090'
Data variables:
    daily_precipitation    (time, station) float64 ...
    start_of_observation   (station) object ...
    elevation              (station) float64 ...
    latitude               (station) float64 24.89 24.97 24.94 ... 23.44 24.12
    longitude              (station) float64 121.4 121.4 121.7 ... 121.3 121.6
    mask_5years            (station) float64 ...
    genpareto_shape_param  (station) float64 0.6229 0.5113 ... 0.6114 0.5293'''
    
    import numpy as np
    
    station_gpd_3hourly = station_data
    latbins = latbins
    lonbins = lonbins
    lon = station_gpd_3hourly.longitude.values
    lat = station_gpd_3hourly.latitude.values
    gpd_shape_params = station_gpd_3hourly.genpareto_shape_param.values

    nsamples, xx, yy = np.histogram2d(lat, lon, bins=(latbins, lonbins))
    gpd_shape_params_sum, xx, yy = np.histogram2d(lat, lon, bins=(latbins, lonbins), weights=gpd_shape_params)
    gpd_shape_params_avg = gpd_shape_params_sum / nsamples
    gpd_shape_params_avg_nonans = np.ma.masked_invalid(gpd_shape_params_avg)
    
    return gpd_shape_params_avg_nonans




def plot_seasonal_pdf_parameters_as_heatmap(data, colorbar_min=None, colorbar_max=None, mode='genpareto_pdf_param=shape', 
                                           fmt='.2f' , annot_kws={"size":4}, annot=False):
    '''Function to plot parameters of particular PDF or anything else that has the right structure.
    <xarray.DataArray '__xarray_dataarray_variable__' (season: 3, genpareto_pdf_param: 2, latitude: 40, longitude: 20)>
             kwargs --> fmt='.2f' , annot_kws={"size":4}):'''
    
    
    import seaborn as sns
    from matplotlib.ticker import FormatStrFormatter
    import matplotlib.ticker as mtick
    import matplotlib.style as style
    style.use('seaborn-poster')
    
    
    gsmap_taiwan_params_3hourly = data.copy()
    
    shape_df = pd.DataFrame(gsmap_taiwan_params_3hourly.sel(genpareto_pdf_param='shape').to_series().reset_index())


    df1 = shape_df[shape_df.season == 'full'].pivot('latitude', 'longitude', '__xarray_dataarray_variable__')[::-1]
    df2 = shape_df[shape_df.season == 'wet'].pivot('latitude', 'longitude', '__xarray_dataarray_variable__')[::-1]
    df3 = shape_df[shape_df.season == 'dry'].pivot('latitude', 'longitude', '__xarray_dataarray_variable__')[::-1]



    df1.columns = np.around(df1.columns, decimals=2)
    df1.index = np.round(df1.index, decimals=2)
    df2.columns = np.around(df2.columns, decimals=2)
    df2.index = np.round(df2.index, decimals=2)
    df3.columns = np.around(df3.columns, decimals=2)
    df3.index = np.round(df3.index, decimals=2)

    fig, axn = plt.subplots(1, 3, dpi=200, sharey=True)
    df_all = pd.concat([df1, df2, df3], axis=1, keys=['full', 'wet', 'dry'])

    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    seasons = ['full', 'wet', 'dry']
    for i, ax in enumerate(axn.flat):

        sns.heatmap(df_all[seasons[i]], ax=ax,
                    cbar=i == 0,
                    vmin=colorbar_min, vmax=colorbar_max, square=True, linewidth=0., annot=annot, fmt=fmt , annot_kws=annot_kws, 
                    cbar_ax=None if i else cbar_ax).figure.axes[-1].yaxis.label.set_size(8) #cbar_kws={'label': 'gamma shape parameter'}



        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)  
            tick.label.set_rotation(s=80)

        for tick1 in ax.yaxis.get_major_ticks():
            tick1.label.set_fontsize(14)  
            tick1.label.set_rotation(s=10)
        ax.set_title('season = '+ seasons[i], fontdict={'fontsize':17})    
        if i!=0:
            ax.get_yaxis().set_visible(False)
        else:
            pass

        ax.set_ylabel('latitude', fontsize=17)
        ax.set_xlabel('longitude', fontsize=17)

    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=16)



    #fig.suptitle('Gamma PDF shape parameter', x=0.5, y=0.8)
    fig.tight_layout(rect=[0, 0, 0.9, 1])
    # fig.savefig(base+'/gsmap_genpareto_shape_param_seasons_taiwan_3hourly.pdf', transparent=True, dpi=200, bbox_inches='tight',
    #             format='pdf')
    return




def pdf_params_diff_temporal_resolutions(dataset, resolutions=['D', '12H', '9H', '6H', '3H', '1H'], 
                                         thresholds=[0.1, 0.1/2, 9*0.1/24, 0.1/4, 0.1/8, 0.1/24]):
    '''Function to loop over different temporal scales (with spatial_pdf_params()) in order to estimate spatial patterns of PDF parameters. It
    appends the arrays into a list and then concatenates everything together. Input is a dataset, where precipitation is stored in dataset.precipitation. The output is DataArray of seasonal PDF
    parameters, BUT with added new dimension "temporal_resolution".
                Coordinates:
              * latitude             (latitude) float64 25.45 25.35 25.25 ... 21.65 21.55
              * longitude            (longitude) float64 120.1 120.2 120.2 ... 121.8 121.9
              * genpareto_pdf_param  (genpareto_pdf_param) object 'shape' 'scale'
              * season               (season) object 'full' 'wet' 'dry'
              * temporal_resolution  (temporal_resolution) <U3 'D' '12H' '9H' '6H' '3H' '1H' '''
    
    import xarray as xr
    import sys
    sys.path.insert(0, "/Users/lisac/Documents/Data_analysis/my_useful_codes/")
    from useful_stuff import spatial_pdf_params
    
    dataarrays = []
    for idx, one_timescale in enumerate(resolutions):
        
        gsmap_resampled = dataset.resample({'time': one_timescale}).sum('time', skipna=False)
        spatial_pdf_params_temp = spatial_pdf_params(gsmap_resampled.precipitation, season='everything', 
                                                     threshold=thresholds[idx], mode='genpareto_pdf')
        spatial_pdf_params_temp = spatial_pdf_params_temp.assign_coords(temporal_resolution=one_timescale).expand_dims('temporal_resolution')

        dataarrays.append(spatial_pdf_params_temp)
    
    concated_da = xr.concat(dataarrays, dim='temporal_resolution')
    concated_da = concated_da.rename('genpareto_pdf_parameters')
    
    return concated_da


def from_shp_to_json(data_dir):     
    '''Function to take all of the .shp files in the data directory, stick them together into one file, and make it 
    .json format'''
    
    import glob
    import geopandas as gpd
    import pandas as pd

    shapefilesList = []
    for shp_file in sorted(glob.glob(data_dir)):
        shapefilesList.append(gpd.read_file(shp_file))

    return gpd.GeoDataFrame(pd.concat(shapefilesList, ignore_index=True))




def transform_from_latlon(lat, lon):
    '''function used in rasterize function'''
    import numpy as np
    from affine import Affine
    
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.
    """
    from rasterio import features
    import xarray as xr
    import numpy as np
    fill = np.nan
    
    import sys
    sys.path.insert(0, "/Users/lisac/Documents/Data_analysis/my_useful_codes/")
    from useful_stuff import transform_from_latlon
    
    transform = transform_from_latlon(coords['latitude'], coords['longitude'])
    out_shape = (len(coords['latitude']), len(coords['longitude']))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    return xr.DataArray(raster, coords=coords, dims=('latitude', 'longitude'))




def rasterize_shapefiles(path_to_shapefiles, coords_for_burning):
    '''Function to rasterize .json format shapefiles, and put each one into separate xarray DataArray layer.
    It uses custom rasterize and transform_from_latlon functions. Input is filepath to .json shapefiles, and 
    DataArray to use as a blueprint for new DataArray's latitude and longitude into which the catchments will be burned.
    Output is a xarray.DataArray: 
    Coordinates:
              * longitude         (longitude) float64 78.05 78.15 78.25 ... 89.85 89.95
              * latitude          (latitude) float64 31.45 31.35 31.25 ... 25.75 25.65 25.55
              * catchment_number  (catchment_number) int64 0 1 2 3 4 5 ... 25 26 27 28 29 30'''
    
    import geopandas as gpd
    import xarray as xr
    import sys
    sys.path.insert(0, "/Users/lisac/Documents/Data_analysis/my_useful_codes/")
    from useful_stuff import rasterize
    
    
    shp_file_large = gpd.read_file(path_to_shapefiles)
    ds = xr.Dataset(coords={'longitude': coords_for_burning.longitude,   # np.linspace(78.05, 89.95, num=120)
                              'latitude': coords_for_burning.latitude})  # np.linspace(31.45, 25.55, num=60)
    
    shapes = [(shape, n) for n, shape in enumerate(shp_file_large.geometry)]

    one_da = rasterize([shapes[0]], ds.coords).assign_coords(catchment_number=shapes[0][1]).expand_dims('catchment_number')

    for geometry_tuple in shapes[1:]:
        one_da = xr.concat([one_da, rasterize([geometry_tuple], ds.coords).assign_coords(catchment_number=geometry_tuple[1]).expand_dims('catchment_number')], dim='catchment_number')

    return one_da





def erosion_model_fscape(mi, ni, b=2, tmax=5*10**6, dt=1*10**5, threshold_option = 'no_threshold', gama = 1, xl = 100000, 
                        yl = 100000, nx=101, ny=101, U=2*10**(-3), K = 10**(-5), n = 1, ):
    
    '''Semimplicit/explicit solution to advective equation in a numerically most efficient way, order O(n). Fastscape + Deal et al. (2018)'''

    m = 0.4*n
    nn = nx*ny
    dx = xl / (nx-1)
    dy = yl / (ny-1)
    h = np.random.rand(nn) # put random seed somewhere
    
    def find_stack(ij, don, ndon, nn, stack, nstack):
        for k in range(ndon[ij]):
            ijk = don[ij, k]
            stack[nstack] = ijk
            nstack += 1
            nstack = find_stack(ijk, don, ndon, nn, stack, nstack)
        return nstack


    
    for t in range(0, tmax, dt):

        print(f'time = {str(t)[:-3]} * 10^3 years')

        rec = np.zeros(nn, dtype=np.int)

        for ij in range(nn):
            rec[ij] = ij

        length = np.zeros(nn)

        for j in range(1, ny - 1):
                for i in range(1, nx - 1):
                    ij = i + nx * j
                    smax=0
                    for jj in range(-1, 2):
                        for ii in range(-1, 2):
                            iii = ii + i
                            jjj = j + jj
                            ijk = iii + jjj*nx
                            if ijk != ij:
                                l = ((dx * ii)**2 + (dy * jj)**2)**0.5
                                slope = (h[ij] - h[ijk]) / l
                                if slope > smax:
                                    smax = slope
                                    rec[ij] = ijk
                                    length[ij] = l

        ndon = np.zeros(nn, dtype=np.int)
        don = np.empty((nn, 8), dtype=np.int)
        stack = np.empty_like(ndon)
        ndon[:] = 0

        for ij in range(nn):
                if rec[ij] != ij:
                    ijk = rec[ij]
                    ndon[ijk] += 1
                    don[ijk, ndon[ijk] - 1] = ij



        nstack = 0

        for ij in range(nn):
                if rec[ij] == ij:
                    stack[nstack] = ij
                    nstack += 1
                    nstack = find_stack(ij, don, ndon, nn, stack, nstack)


        area = np.empty(nn)
        area[:] = dx * dy                        


        for ijk in stack[-1::-1]:
                if rec[ijk] != ijk:
                    area[rec[ijk]] += area[ijk]


        for j in range(1, ny-1):
                for i in range(1, nx-1):
                    ij = i + j*nx
                    h[ij] = h[ij] + U*dt

        q = np.zeros(nn)
        mi_e = np.ones(nn)
        lambda_e = np.zeros(nn)
        Ke = np.zeros(nn)
        slopes = np.zeros(nn)
        C1 = np.zeros(nn)
        C2 = np.zeros(nn)
        q_c = np.zeros(nn)

        for ij in range(nn):
                ijk = stack[ij]
                ijr = rec[ijk]
                if ijr == ijk:
                    continue
                if threshold_option == 'no_threshold':
                    h_0 = h[ijk]

                    if b==1:
                        mi_e[ijk] = sp.special.gamma(1/ni + gama) / (sp.special.gamma(1/ni)*ni**(-gama))
                    elif b==2:
                        mi_e[ijk] = sp.special.gamma(1/ni + 1 - gama) / (sp.special.gamma(1/ni)*ni**(gama-1))
                    else:
                        pass

                    Ke[ijk] = mi_e[ijk] * mi**m * K
                    F = Ke[ijk] * dt * area[ijk]**m / length[ijk]**n

                    h[ijk] = (h_0 + F*h[ijr]) / (1. + F)

                elif threshold_option == 'variable_Shields_number':
                    h_0 = h[ijk]
                    slopes[ijk] = (h[ijk] - h[ijr]) / length[ijk]
                    q_c[ijk] = q_c_func(slope=slopes[ijk], area=area[ijk])
                    mi_e[ijk] = mi_e_func(ni, q_c[ijk], b=1)
                    lambda_e[ijk] = lambda_e_func(ni, q_c[ijk], b=1)
                    C1[ijk] = (mi_e[ijk] * K * mi**m * dt * area[ijk]**m) / length[ijk]**n
                    C2[ijk] = (k_e * lambda_e[ijk] * D**a_t * 0.15**a_t * dt) / length[ijk]**(0.25*a_t)
                    #continue
                    for _ in range(10):
                        num = h[ijk] - h_0 + C1[ijk]*(h[ijk] - h[ijr])**n - C2[ijk]*(h[ijk] - h[ijr])**(0.25*a_t)
                        denum = 1 + n*C1[ijk]*(h[ijk] - h[ijr])**(n-1) - C2[ijk]*0.25*a_t*(h[ijk] - h[ijr])**(0.25*a_t - 1)
                        h[ijk] = h[ijk] - num/denum



                else:

                    pass
                
    return h


def spatial_chi_statistic(dataarray, number_of_k_intervals=5, threshold=0.1, mode='genpareto'):
      '''  Function for estimation of Chi square statistic for Maximum likelihood estimator (MLE) per grid cell, 
            using the method with bootstrapping. 
            Input is 3-dim (time, latitude, longitude) xarray dataarray.
            Threshold is to cut the rainfall data. 
            Number of k intervals is the number of bins in histogram.  '''
    
    
    import numpy as np
    from scipy import stats
    
    rainfall_data = dataarray.where(dataarray>threshold).copy()
    chi_statistic_ar = np.empty((len(rainfall_data.latitude), len(rainfall_data.longitude)))
    chi_statistic_ar[:] = np.nan


    expected_counts = []
    observed_counts = []

    
    if mode == 'gamma':
        for i in range(len(rainfall_data.latitude)):
            for j in range(len(rainfall_data.longitude)):

                gamma_data = rainfall_data.isel(latitude=i, longitude=j).dropna('time')

                size_of_data = len(gamma_data.time)
                size_of_samples = size_of_data
                resample_gamma = np.random.choice(gamma_data.values, size=size_of_samples, replace=True)

                alpha_fit, loc_fit, beta_fit = stats.gamma.fit(resample_gamma, floc=0)

                k_bins = np.empty((number_of_k_intervals, len(gamma_data.time)))
                k_bins[:] = np.nan

                k_intervals = np.linspace(0, 1, number_of_k_intervals+1)



                for idx_element, element in enumerate(gamma_data.values):
                    for bin_ in range(len(k_intervals)-1):
                        if (stats.gamma.cdf(element, alpha_fit, loc=0, scale=beta_fit) >= k_intervals[bin_]) & (stats.gamma.cdf(element, alpha_fit, loc=0, scale=beta_fit) < k_intervals[bin_+1]):
                            k_bins[bin_][idx_element] = element
                            break

                del expected_counts
                expected_counts = []           

                for _bin in range(len(k_intervals)-1):
                    expected_counts.append((k_intervals[1] - k_intervals[0]) * len(gamma_data.values))


                del observed_counts
                observed_counts = []
                observed_counts.append(np.sum(~np.isnan(k_bins), axis=1))


                chi_statistic_ar[i,j] = np.sum((np.array(observed_counts) - np.array(expected_counts))**2 / np.array(expected_counts))
                


        return (chi_statistic_ar, observed_counts, expected_counts)
    
    
    
    if mode == 'expon':
        for i in range(len(rainfall_data.latitude)):
            for j in range(len(rainfall_data.longitude)):

                gamma_data = rainfall_data.isel(latitude=i, longitude=j).dropna('time')

                size_of_data = len(gamma_data.time)
                size_of_samples = size_of_data
                resample_gamma = np.random.choice(gamma_data.values, size=size_of_samples, replace=True)

                loc_fit, alpha_fit = stats.expon.fit(resample_gamma, floc=0)

                k_bins = np.empty((number_of_k_intervals, len(gamma_data.time)))
                k_bins[:] = np.nan

                k_intervals = np.linspace(0, 1, number_of_k_intervals+1)



                for idx_element, element in enumerate(gamma_data.values):
                    for bin_ in range(len(k_intervals)-1):
                        if (stats.expon.cdf(element, loc=0, scale=alpha_fit) >= k_intervals[bin_]) & (stats.expon.cdf(element, loc=0, scale=alpha_fit) < k_intervals[bin_+1]):
                            k_bins[bin_][idx_element] = element
                            break

                del expected_counts
                expected_counts = []           

                for _bin in range(len(k_intervals)-1):
                    expected_counts.append((k_intervals[1] - k_intervals[0]) * len(gamma_data.values))


                del observed_counts
                observed_counts = []
                observed_counts.append(np.sum(~np.isnan(k_bins), axis=1))


                chi_statistic_ar[i,j] = np.sum((np.array(observed_counts) - np.array(expected_counts))**2 / np.array(expected_counts))



        return (chi_statistic_ar, observed_counts, expected_counts)







    if mode == 'genpareto':
        for i in range(len(rainfall_data.latitude)):
            for j in range(len(rainfall_data.longitude)):

                gamma_data = rainfall_data.isel(latitude=i, longitude=j).dropna('time')

                size_of_data = len(gamma_data.time)
                size_of_samples = size_of_data
                resample_gamma = np.random.choice(gamma_data.values, size=size_of_samples, replace=True)

                shape_fit, loc_fit, scale_fit = stats.genpareto.fit(resample_gamma, floc=0)

                k_bins = np.empty((number_of_k_intervals, len(gamma_data.time)))
                k_bins[:] = np.nan

                k_intervals = np.linspace(0, 1, number_of_k_intervals+1)



                for idx_element, element in enumerate(gamma_data.values):
                    for bin_ in range(len(k_intervals)-1):
                        if (stats.genpareto.cdf(element, shape_fit, loc=0, scale=scale_fit) >= k_intervals[bin_]) & (stats.genpareto.cdf(element, shape_fit, loc=0, scale=scale_fit) < k_intervals[bin_+1]):
                            k_bins[bin_][idx_element] = element
                            break

                del expected_counts
                expected_counts = []           

                for _bin in range(len(k_intervals)-1):
                    expected_counts.append((k_intervals[1] - k_intervals[0]) * len(gamma_data.values))


                del observed_counts
                observed_counts = []
                observed_counts.append(np.sum(~np.isnan(k_bins), axis=1))


                chi_statistic_ar[i,j] = np.sum((np.array(observed_counts) - np.array(expected_counts))**2 / np.array(expected_counts))



        #return (chi_statistic_ar, observed_counts, expected_counts, (np.array(observed_counts) - np.array(expected_counts))**2 / np.array(expected_counts), shape_fit, scale_fit)
        
        return chi_statistic_ar #, observed_counts, expected_counts, (np.array(observed_counts) - np.array(expected_counts))**2 / np.array(expected_counts), shape_fit, scale_fit)




def concat_thur_stations_with_metadata(path_to_hydrocode_output_files, mode='my_areas'):
    '''input: path to folder where you have .nc hydrocode output station files, in order to concat them to one file 
       with added dimension station, and also to add some metadata, e.g. drainage areas...'''
    
    
    import xarray as xr
    
    list_of_stations = []
    names = []
    nc_paths = []
    attrs = dict()
    
    areas = {'Thur-Andelfingen': 1700, 'Sitter-Appenzell': 74, 'Murg-Waengi': 80, 'Thur-Halden': 1085,
            'Thur-Jonschwil': 493, 'Glatt-Herisau': 17}

    my_areas = {'Thur-Andelfingen': 1710.98, 'Sitter-Appenzell': 88.02, 'Murg-Waengi': 76.82, 'Thur-Halden': 1103.56, 
                'Thur-Jonschwil': 493.13, 'Glatt-Herisau': 16.51}

    station_areas = xr.DataArray([(areas[i]) for i in areas], dims=('station',), coords={'station': [i for i in areas]}, name='areas')
    station_reported_areas = xr.DataArray([(my_areas[i]) for i in my_areas], dims=('station',), coords={'station': [i for i in my_areas]}, name='my_areas')



    for nc_file in sorted(glob.glob(path_to_hydrocode_output_files)):
        dis_station_da = xr.open_dataset(nc_file)
        list_of_stations.append(dis_station_da)
        names.append(nc_file.split('/')[-1].split('_')[1])
        nc_paths.append(nc_file)
        attrs[nc_file.split('/')[-1].split('_')[1]] = dis_station_da.attrs

    hydro_results_calendar_seasons = xr.concat(list_of_stations, dim=pd.Index(names, name='station'))
    hydro_results_calendar_seasons = xr.merge([hydro_results_calendar_seasons, station_areas, station_reported_areas])

    hydro_results_calendar_seasons.seasonId.values = [i.decode('ascii') for i in hydro_results_calendar_seasons.seasonId.values]
    hydro_results_calendar_seasons.stationId.values = [i.decode('ascii') for i in hydro_results_calendar_seasons.stationId.values]

    hydro_results_calendar_seasons.to_netcdf(base + '_'.join(nc_paths[0].split('_')[-3:]), mode='w')
    
    attrs_df = pd.DataFrame(attrs)
    attrs_df.to_csv(base+f"attrs_{nc_file.split('_')[-1][:-3]}.csv")

    
    return print('Done!')




def P_and_lambda_time_rasters(path_to_dems_for_masks, path_to_rainfall_interpolations_npy, start_year = 1970
                             ):
    '''Function to transform values that fall inside particular catchment area, bounded by the DEM.
    station variable: /Users/lisac/Documents/Data_analysis/Thur/DEMs/nested_catchments_utm/waengi_dem25_utm.txt,
     path_to_dems_for_masks: base + /*.txt  
     path_to_rainfall_interpolations_npy: base + /*.npy''' 
    
    
    import glob
    import numpy as np
    import pandas as pd
    
    names = []
    P_sums_stations = []
    
    for station in sorted(glob.glob(path_to_dems_for_masks)): # DEMs in .txt, of stations from folder for masks
        print(f'Starting with: {station}')
        dem = np.loadtxt(station)
        h = dem.flatten()
        masknan = h==-9999
        P_sums = [] # "station" variable: '/Users/lisac/Documents/Data_analysis/Thur/DEMs/nested_catchments_utm/waengi_dem25_utm.txt'
        names.append(station.split('/')[-1].split('_')[0].capitalize())
        for year in sorted(glob.glob(path_to_rainfall_interpolations_npy)): # year range
            rainfall_year = np.load(year)
            print(f'Done loading: {year}')
            for day in rainfall_year:
                rain_sum = np.sum(day.flatten() * ~masknan) / (~masknan).sum() / 100.
                P_sums.append(rain_sum)

            del rainfall_year


        P_sums_stations.extend([P_sums])
        print(' ')

    end_year = int(year[-8:-4])
    columns = pd.date_range(start=f'01-01-{start_year}', end=f'31-12-{end_year}', freq='D')
    
    
    P_sums_stations_df = pd.DataFrame(P_sums_stations, index=names, columns=columns).transpose()
    
    
    
    print('Done with one function run!')
    return P_sums_stations_df



def calculate_alpha_from_stations(path_to_rainfall_station_data, rain_threshold=0.1):
    '''Function to calculate observed seasonal values of mean rainfall intensity from rainfall station data.
     Input is .txt of rainfall station data, the same format used to interpolate on daily basis. '''
    
    import pandas as pd
    import numpy as np

    all_stations_rain = pd.read_csv(path_to_rainfall_station_data, index_col=0, parse_dates=[0])

    summer = pd.date_range('21-06-1970', '21-09-1970', freq='D', closed='left').dayofyear.values
    fall = pd.date_range('21-09-1970', '21-12-1970', freq='D', closed='left').dayofyear.values
    winter = pd.date_range('21-12-1970', '21-03-1971', freq='D', closed='left').dayofyear.values
    spring = pd.date_range('21-03-1971', '21-06-1971', freq='D', closed='left').dayofyear.values
    seasons = [spring, summer, fall, winter]
    season_names = ['spring', 'summer', 'fall', 'winter']
    catchments = {'Waengi': ['Eschlikon'], 
                  'Appenzell': ['Saentis'], 
                  'Herisau': ['Herisau'], 
                  'Jonschwil': ['Ricken', 'Sulgen', 'Starkenbach'],
                  'Halden': ['Saentis', 'Herisau', 'Ricken', 'Sulgen', 'Starkenbach', 'Urnaesch', 'Appenzell', 'Teufen',
                             'Flawil', 'Bischofszell'],
                  'Andelfingen': ['Saentis', 'Herisau', 'Ricken', 'Sulgen', 'Starkenbach', 'Urnaesch', 'Appenzell', 
                                  'Teufen', 'Flawil', 'Bischofszell', 'Eschlikon', 'Affeltrangen', 'Andelfingen', 
                                  'Frauenfeld', 'Illhart', 'Kalchrain', 'Niederneunforn', 'St.Peterzell', 'Weinfelden']}

    for idx_catch,catchment in enumerate(catchments):
        all_stations_rain_temp = all_stations_rain.loc[:,catchments[catchment]].copy()
        all_stations_rain_thresh = all_stations_rain_temp[all_stations_rain_temp>rain_threshold].copy()

        one_catchment = []
        for idx,season in enumerate(seasons):
            one_catchment.append([season_names[idx],
                                  all_stations_rain_thresh[np.in1d(all_stations_rain_thresh.index.dayofyear, season)].mean(axis=1).mean()])
                                  

        if idx_catch==0:
            all_catchments_df = pd.DataFrame(one_catchment, columns=['season', catchment])

        else:
            one_catchment_df = pd.DataFrame(one_catchment, columns=['season', catchment])
            all_catchments_df = pd.merge(all_catchments_df, one_catchment_df)


    all_catchments_df['rain threshold'] = None
    all_catchments_df['rain threshold'] = str(rain_threshold) + ' mm'
    all_catchments_df = pd.melt(all_catchments_df, id_vars=['season', 'rain threshold'], value_vars=None, 
                                var_name='station', value_name='observed mean rainfall intensity [mm] (E)')


    return all_catchments_df



def calculate_alpha_from_interpolated_P_rasters(P_timeseries_rasters_datetimeindex_dataframe, rain_threshold=0.1, year_start=1970, year_end=1979):
    '''Calculate mean rainfall intensity from interpolated station data, and sum up in the contributing 6 catchments
        for Thur region.'''
    
    import pandas as pd
    import numpy as np
    
    P_rasters = P_timeseries_rasters_datetimeindex_dataframe.copy()
    P_rasters = P_rasters.loc[str(year_start):str(year_end)]
    lambda_p_rasters = (P_rasters > rain_threshold).copy()


    seasonal_rasters = []
    seasonal_rasters_lambda = []


    summer = pd.date_range('21-06-1970', '21-09-1970', freq='D', closed='left').dayofyear.values
    fall = pd.date_range('21-09-1970', '21-12-1970', freq='D', closed='left').dayofyear.values
    winter = pd.date_range('21-12-1970', '21-03-1971', freq='D', closed='left').dayofyear.values
    spring = pd.date_range('21-03-1971', '21-06-1971', freq='D', closed='left').dayofyear.values
    dry = np.append(np.arange(220,366), np.arange(0, 26))
    wet = np.arange(70,176)
    annual = np.arange(0,366)

    seasons = [spring, summer, fall, winter]
    season_names = ['spring', 'summer', 'fall', 'winter']


    for idx, season in enumerate(seasons):
        P_seasons_df = pd.DataFrame(P_rasters[np.in1d(P_rasters.index.dayofyear, season)].mean(axis=0),
                               columns=[season_names[idx]])
        lambda_p_seasons_df = pd.DataFrame(lambda_p_rasters[np.in1d(lambda_p_rasters.index.dayofyear, season)].mean(axis=0),
                               columns=[season_names[idx]])


        P_seasons_df.index.name = 'station'
        lambda_p_seasons_df.index.name = 'station'

        seasonal_rasters.append(P_seasons_df)
        seasonal_rasters_lambda.append(lambda_p_seasons_df)

        P_seasons_df = pd.concat(seasonal_rasters, axis=1)
        lambda_p_seasons_df = pd.concat(seasonal_rasters_lambda, axis=1)


    alpha_seasons_df = P_seasons_df / lambda_p_seasons_df
    alpha_seasons_df['rain threshold'] = None
    alpha_seasons_df['rain threshold'] = str(rain_threshold) + ' mm'
    alpha_seasons_df = alpha_seasons_df.reset_index()
    alpha_seasons_df = pd.melt(alpha_seasons_df, id_vars=['station', 'rain threshold'], value_vars=None, 
                                var_name='season', value_name='modeled mean rainfall intensity [mm] (E)')
    
    P_seasons_df['rain threshold'] = None
    P_seasons_df['rain threshold'] = str(0.1) + ' mm'
    P_seasons_df = P_seasons_df.reset_index()
    P_seasons_df = pd.melt(P_seasons_df, id_vars=['station', 'rain threshold'], value_vars=None, 
                            var_name='season', value_name='modeled mean rainfall rate [mm/day] (E)')

    lambda_p_seasons_df['rain threshold'] = None
    lambda_p_seasons_df['rain threshold'] = str(0.1) + ' mm'
    lambda_p_seasons_df = lambda_p_seasons_df.reset_index()
    lambda_p_seasons_df = pd.melt(lambda_p_seasons_df, id_vars=['station', 'rain threshold'], value_vars=None, 
                            var_name='season', value_name='modeled mean rainfall frequency [1/day] (E)')

    
    return (P_seasons_df, lambda_p_seasons_df, alpha_seasons_df)



def compare_plots_of_station_and_interpolated_data_at_the_position_of_stations(station_name, year_start=1970, year_end=1970, 
                                                                              mode='mm'):
    '''Function to compare visual representations of station and interpolated data'''
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    # row and colum position of the weather station in the interpolated raster
    weather_stations_dict = {'Saentis': [1794, 1996], 'Illhart': [142, 1075], 'Kalchrain': [190, 623], 
                             'Andelfingen': [229, 13], 'Niederneunforn': [257, 294]}
    
    base = '/Users/lisac/Documents/Data_analysis/Thur/'
    all_years = []
    
    weather_station = pd.read_csv(base + f'weather_data/processed/{station_name}.csv', index_col=0, parse_dates=[0], 
                                  usecols=[0,2])
    if mode=='cm':
        weather_station = 0.1*weather_station

    
    
    for year in list(range(year_start, year_end+1)):
        if mode=='cm':
            fp_interp = base + f'interpolations/with_cm/interpolation_idw_dem25_utm_{year}.npy'
        else:  
            fp_interp = base + f'interpolations/interpolation_idw_dem25_utm_{year}.npy'
        
        year_interp = np.load(fp_interp)
        year_interp = year_interp/100.

        interp_station_lt = list(year_interp[:, weather_stations_dict[station_name][0], 
                                     weather_stations_dict[station_name][1]])
        all_years.extend(interp_station_lt)
        
    fig = plt.figure(figsize=(15,7))   
    plt.plot(weather_station.loc[str(year_start):str(year_end)].values, color='red', label='station data')

    plt.plot(all_years, marker='o', linestyle=':', color='k', 
             label='interpolated data at the position of the station')

    plt.legend()
    
    
    
    return print(f'Position of weather station {station_name} for the time span 01.01.{year_start}.-31.12.{year_end}.')





def estimate_lambda_from_discharge_timeseries(path_to_discharge_stations_csv_files, parametar='modeled mean rainfall intensity [mm] (E)',
                                             parametar1='observed mean rainfall intensity [mm] (E)'):
    '''Calculate seasonal mean discharge and discharge frequency from discharge and climate station data.
       (base+'hydrologic_data/processed/specific_discharge_my_areas/*.csv')'''
    
    
    
    import glob
    import pandas as pd
    import numpy as np

    
    base = '/Users/lisac/Documents/Data_analysis/Thur/' 
    alpha_data = pd.read_csv(base + 'weather_data/working_data/thur_observed_modeled_my_and_paper_values.csv')
    alpha_modeled = alpha_data.where(alpha_data['rain threshold'] == '0.1 mm').dropna()


    summer = pd.date_range('21-06-1970', '21-09-1970', freq='D', closed='left').dayofyear.values
    fall = pd.date_range('21-09-1970', '21-12-1970', freq='D', closed='left').dayofyear.values
    winter = pd.date_range('21-12-1970', '21-03-1971', freq='D', closed='left').dayofyear.values
    spring = pd.date_range('21-03-1971', '21-06-1971', freq='D', closed='left').dayofyear.values
    seasons = [spring, summer, fall, winter]
    season_names = ['spring', 'summer', 'fall', 'winter']

    

    one_station_lt = []

    for csv_station_file in sorted(glob.glob(path_to_discharge_stations_csv_files)):
        station_temp = pd.read_csv(csv_station_file, index_col=[0], parse_dates=[[0,1,2]])
        #break

        for idx,season in enumerate(seasons):
            station_season = station_temp[np.in1d(station_temp.index.dayofyear, season)]        
            
            station_offset = np.append(station_season.values[1:], np.nan)
            station_offset = np.expand_dims(station_offset, axis=1)
            mean_lambda = (station_offset > station_season.values)[:-1].mean()
            
            season_mean = station_season.mean().values
            alpha_season = alpha_modeled[(alpha_modeled['station'] == csv_station_file.split('-')[1].split('_')[0]) & (alpha_modeled['season'] == season_names[idx])][parametar].values[0]
            alpha_season1 = alpha_modeled[(alpha_modeled['station'] == csv_station_file.split('-')[1].split('_')[0]) & (alpha_modeled['season'] == season_names[idx])][parametar1].values[0]

            season_lambda = season_mean[0]/alpha_season
            season_lambda1 = season_mean[0]/alpha_season1
            
            
            
            one_station_lt.append([csv_station_file.split('-')[1].split('_')[0], season_names[idx], season_mean[0], mean_lambda, season_lambda, season_lambda1])
            





    all_lambda_E = pd.DataFrame(one_station_lt, columns=['station', 'season', 'observed mean discharge [mm/day]', 'mean streamflow frequency [1/d] (E1)', 
                                'mean streamflow frequency [1/d] (E2)', 'mean streamflow frequency [1/d] (E3)'])

    

    lambda_paper = pd.read_csv(base + 'weather_data/working_data/thur_paper_mean_discharge_frequency_values.csv')
    lambda_merged = pd.merge(all_lambda_E, lambda_paper)
    
    
    
    
    
    return (lambda_merged)


def calculate_catchment_scale_seasonal_pet():
    '''Estimating catchment-scale Potential evapotranspiration (PET) from CGIAR data.'''
    
    import numpy as np
    import glob
    import xarray as xr
    from calendar import monthrange
    import pandas as pd

    base = '/Users/lisac/Documents/Data_analysis/Thur/'
    fp_dem_txt = base + 'DEMs/nested_catchments_utm/*.txt'
    cgiar_months_all = xr.open_dataset(base+'cgiar_pet_monthly/cgiar_pet_monthly_averages_thur_utm.nc').pet
    station_names = []
    PET_sums_stations = []
    cgiar_seasons = {'winter': [1, 2, 3], 'spring': [4, 5, 6], 'summer': [7, 8, 9], 'fall': [10, 11, 12]}

    for station in sorted(glob.glob(fp_dem_txt)): # DEMs in .txt, of stations from folder for masks
        print(f"Starting with: {station.split('/')[-1]}")
        dem = np.loadtxt(station)
        h = dem.flatten()
        masknan = h==-9999
        station_names.append(station.split('/')[-1].split('_')[0].capitalize())
        PET_sums = []
        season_names = []
        for season in cgiar_seasons:
            #print(f'current season is {season}')

            cgiar_season_sum = cgiar_months_all.sel(month=cgiar_seasons[season]).sum('month')

            if season == 'winter':
                number_of_days = monthrange(2011, cgiar_seasons[season][0])[1] + monthrange(
                    2011, cgiar_seasons[season][1])[1] + monthrange(2011, cgiar_seasons[season][2])[1] + 0.25
            else:
                number_of_days = monthrange(2011, cgiar_seasons[season][0])[1] + monthrange(
                    2011, cgiar_seasons[season][1])[1] + monthrange(2011, cgiar_seasons[season][2])[1]


            #print(f'Number of days {number_of_days}')
            seasonal_pet = cgiar_season_sum/number_of_days
            pet_sum = np.sum(seasonal_pet.values.flatten() * ~masknan) / (~masknan).sum()
            PET_sums.append(pet_sum)
            season_names.append(season)


        PET_sums_stations.append(PET_sums)


    PET_sums_stations_df = pd.DataFrame(PET_sums_stations, columns=season_names, index=station_names)
    PET_sums_stations_df.index.names = ['station']
    PET_sums_stations_df = PET_sums_stations_df.reset_index()
    PET_sums_stations_df = pd.melt(PET_sums_stations_df, id_vars='station', value_vars=None, var_name='season', value_name='catchment-scale PET [mm/day]')
    
    
    return PET_sums_stations_df




def create_correlated_timeseries_M1(uncorrelated_rainfall_daily_data, g=1):
    '''Method 1 to create artificial (fake) correlated rainfall timeseries on a landscape. It averages the values in a window, 
    i.e. moving square. g is the number of layers around a central square, which value is being calculated by averaging.'''
    
    bottom_boundary = -g
    top_boundary = 1 + g

    smoothed_averages = []
    

    for j in range(g, ny-g):
        for i in range(g, nx-g):
            ij = i + j*nx
            
            square_window = []
            for jj in range(bottom_boundary, top_boundary):
                for ii in range(bottom_boundary, top_boundary):
                    iii = i + ii
                    jjj = j + jj
                    ijk = iii + jjj*nx
                    fake_data_flatten = uncorrelated_rainfall_daily_data.flatten()
                    square_window.append(fake_data_flatten[ijk])

            smoothed_averages.append(np.mean(square_window))

    fake_correlated_data = np.reshape(smoothed_averages, (nx-2*g , ny-2*g))
    
    return fake_correlated_data




def create_correlated_timeseries_M2(nx=100, ny = 100, n = 18, timeseries_length = 2):
    '''Method 2 for creating correlated fake timeseries on a landscape. It's the moving window method.
    n is the size of the moving square given in number of pixels.'''
    
    timeseries = np.zeros((timeseries_length, ny, nx))


    for one_day in range(timeseries_length):
        starty, endy = 0, n
        for column in range(int(np.ceil(ny/n))):
            startx = 0
            endx = n
            for row in range(int(np.ceil(nx/n))):
                timeseries[one_day, starty:endy, startx:endx] = stats.expon.rvs(loc=0, scale=scale_param, size=1)
                startx += n
                endx += n


            starty += n
            endy += n
            
    return timeseries      