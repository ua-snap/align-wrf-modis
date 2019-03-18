def melt_year( df, metric, station_name ):
    month_days = {'Oct':31, 'Nov':30, 'Dec':31, 'Jan':31, 'Feb':28,
                 'Mar':31, 'Apr':30, 'May':31, 'Jun':30,'Jul':31, 
                 'Aug':31, 'Sep':30}
    month_ids = {'Oct':10, 'Nov':11, 'Dec':12, 'Jan':1, 'Feb':2, 
                'Mar':3, 'Apr':4, 'May':5, 'Jun':6,'Jul':7, 
                'Aug':8, 'Sep':9}

    nrows,ncols = df.shape
    year = df.Year.astype(int).iloc[0]  
    months = [i for i in df.columns if 'Year' not in i and 'Day' not in i]
    out = []
    for month in months:
        if month == 'Feb' and calendar.isleap(year):
            ndays = 29
        else:
            ndays = month_days[month]

        if month in df.columns:
            cur_month = df[['Year','Day',month]]
            cur_month = cur_month[ cur_month.Day <= ndays ]
            cur_month['Month'] = month_ids[month]
            # update colname
            cur_month.columns = [station_name if month in i else i for i in cur_month.columns]
            if month in [10,11,12]:
                # make due to weird formatting the first 3 months are the preceding year
                cur_month['Year'] = cur_month['Year'] -1
            out = out + [cur_month]

    out_df = pd.concat(out)
    # make a DatetimeIndex
    index = pd.DatetimeIndex([ pd.Timestamp('{}-{}-{}'.format(*x[['Year','Month','Day']].astype(int))) for row,x in out_df.iterrows()])
    out_df.index = index
    return out_df[station_name].asfreq('1D')


if __name__ == '__main__':
    import os, calendar
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point
    import numpy as np

    base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
    fn = os.path.join(base_path, 'STATION_DATA', 'snotel', 'SNOTEL_points_all.csv')
    output_path = os.path.join(base_path,'STATION_DATA','snotel','prepped')
    
    if not os.path.exists(output_path):
        _ = os.makedirs(output_path)

    downloaded_path = os.path.join(base_path, 'STATION_DATA', 'snotel', 'raw')
    if not os.path.exists(downloaded_path):
        os.makedirs(downloaded_path)

    df = pd.read_csv( fn )
    df['geometry'] = df.apply( lambda x: Point(x.Longitude, x.Latitude), axis=1 )
    df['ID'] = df.ID.apply(lambda x: x.replace('\t', '').replace(',', ''))
    df['Name'] = [''.join([j for j in i if j.isalnum()]) for i in df.Name ]

    # subset to remove some junk
    df = df[['Name','State' ,'ID', 'Longitude', 'Latitude', 'geometry']]

    # make a shapefile with all of the data...
    pts = gpd.GeoDataFrame(df, geometry='geometry')
    pts.to_file(os.path.join(base_path, 'STATION_DATA', 'snotel', 'SNOTEL_points_all.shp'))
    
    # read in the huc polygon we are using for the test-case.
    huc_fn = os.path.join(base_path, 'shapefiles','chena_river_huc_19040506.shp')
    huc = gpd.read_file( huc_fn )
    huc = huc.to_crs( epsg=4326 )

    # now find which ones intersect
    ak_pts = pts[pts.geometry.apply(lambda x: x.intersects(huc.geometry[0]))]
    ak_pts.to_file(os.path.join(base_path, 'STATION_DATA', 'snotel', 'SNOTEL_points_chena_river_huc.shp'))
    
    # get the data
    import requests
    station_files = []
    for station_id in ak_pts.ID:
        # get the station name from the ID
        station_name = ak_pts[ak_pts.ID == station_id].Name
        station_name = ''.join([i for i in station_name if i.isalnum()])

        # download the data from the server
        url = 'https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customGroupByMonthReport/daily/{0}:AK:SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/TMIN::value,TMAX::value,TAVG::value'.format(station_id)
        response = requests.get(url)
        out_fn = os.path.join(downloaded_path, '{}.csv'.format(str(station_id)))
        with open(out_fn, 'w') as f:
            f.write( response.text )
        station_files = station_files + [out_fn]

    # NOW DO STUFF TO THE DOWNLOADED DATA TO MAKE IT CLEANER
    # separate variables and restructure for processing similar to the ACIS data
    # initial column-fu to deal with the nonsensical data format they are passing.
    colnames = ['Year','Day','Oct-TMIN','Oct-TMAX','Oct-TAVG','Nov-TMIN','Nov-TMAX','Nov-TAVG','Dec-TMIN','Dec-TMAX','Dec-TAVG','Jan-TMIN','Jan-TMAX',
                'Jan-TAVG','Feb-TMIN','Feb-TMAX','Feb-TAVG','Mar-TMIN','Mar-TMAX','Mar-TAVG','Apr-TMIN','Apr-TMAX','Apr-TAVG','May-TMIN','May-TMAX',
                'May-TAVG','Jun-TMIN','Jun-TMAX','Jun-TAVG','Jul-TMIN','Jul-TMAX','Jul-TAVG','Aug-TMIN','Aug-TMAX','Aug-TAVG','Sep-TMIN','Sep-TMAX','Sep-TAVG']

    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    for metric in ['TMAX', 'TMIN', 'TAVG']:
        out_data = []
        for station_fn in station_files:

            station_id = os.path.basename(station_fn).split('.')[0]

            station_name = ak_pts[ak_pts.ID == station_id].Name
            station_name = ''.join([i for i in station_name if i.isalnum()])
            print(station_name)
            df = pd.read_csv(station_fn, skiprows=56, header=None, names=colnames )
            
            cols = ['Year','Day'] + [ col for col in colnames if metric in col ]
            sub_df = df[ cols ]
            sub_df.columns = [i.split('-')[0] for i in sub_df.columns]
            sub_df = sub_df[['Year','Day'] + months]
            
            # melt it and make the index a datetime index
            month_ids = {'Oct':10, 'Nov':11, 'Dec':12, 'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6,'Jul':7, 'Aug':8, 'Sep':9}
            month_days = {'Oct':31, 'Nov':30, 'Dec':31, 'Jan':31, 'Feb':28,
                     'Mar':31, 'Apr':30, 'May':31, 'Jun':30,'Jul':31, 
                     'Aug':31, 'Sep':30}

            melted = pd.concat([ melt_year( df, metric, station_name ) for year, df in sub_df.groupby('Year')])
            out_data = out_data + [melted]

        # write to csv
        out_df = pd.concat(out_data, axis=1)
        # rescale the temps to celcius from fahrenheit
        out_df = ((out_df - 32.0) * 5.0 / 9.0).round(1) # make celcius
        out_fn = os.path.join(output_path,'{}_snotel_stationdata.csv'.format(metric.lower()) )
        out_df.to_csv( out_fn )
