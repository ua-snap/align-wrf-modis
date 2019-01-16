def make_mask( shp_fn, meta ):
	from rasterio.features import rasterize
	shp = gpd.read_file(shp_fn)
	arr = rasterize(list(shp.geometry), out_shape=(meta['height'], meta['width']), fill=0, transform=meta['transform'], 
				all_touched=False, default_value=1 )
	return arr

def extract_mean_aoi( fn, mask ):
	with rasterio.open(fn) as rst:
		arr = rst.read(1)
	# return arr[(mask == 1)].mean()
	return arr[(mask == 1) & (arr != rst.nodata)].mean()

def make_datetime( fn ):
	timestr = os.path.basename(fn).split('_')[-2]
	return datetime.datetime.strptime( timestr,'%Y%j' )

def make_datetime_lst( fn ):
	timestr = os.path.basename(fn).split('_')[1]
	return datetime.datetime.strptime( timestr,'A%Y%j' )


if __name__ == '__main__':
	import matplotlib
	matplotlib.use('agg')
	from matplotlib import pyplot as plt
	import os, glob, itertools, functools, datetime
	import xarray as xr
	import rasterio
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	import multiprocessing as mp

	# input args
	ncpus = 32
	shp_fn = '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_19040506.shp'
	lst_prepped = '/workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos'
	tsk_prepped = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_20km_3338_lanczos'

	# groups
	models = ['GFDL-CM3','NCAR-CCSM4','ERA-Interim']
	scenarios = ['rcp85']
	sensors = ['MOD11A2']
	metrics = ['mean','min','max']
	
	df_list = []
	for model, scenario, sensor, metric in itertools.product(models,scenarios,sensors,metrics):
		if model == 'ERA-Interim':
			scenario = 'historical'

		files = sorted(glob.glob(os.path.join(tsk_prepped,'*{}*{}*{}*{}*.tif'.format(model,scenario,metric,sensor))))

		with rasterio.open(files[0]) as tmp:
			meta = tmp.meta.copy()

		# get the dates from the filenames 
		datetimes = pd.DatetimeIndex([ make_datetime(fn) for fn in files ])

		# make mask
		mask = make_mask( shp_fn, meta )
		f = functools.partial( extract_mean_aoi, mask=mask )
		pool = mp.Pool( ncpus )
		out = pool.map( f, files )
		pool.close()
		pool.join()

		df_list	= df_list + [pd.DataFrame( {'tsk_{}_{}_{}'.format(model,scenario,metric):out}, index=datetimes )]
		
	# now get the lst values...
	for sensor in sensors:
		files = sorted(glob.glob(os.path.join(lst_prepped,'*{}*.tif'.format(sensor))))
		# get the dates from the filenames 
		datetimes = [ make_datetime_lst(fn) for fn in files ]

		# make mask
		mask = make_mask( shp_fn, meta )
		f = functools.partial( extract_mean_aoi, mask=mask )
		pool = mp.Pool( ncpus )
		out = pool.map( f, files )
		pool.close()
		pool.join()

		df = pd.DataFrame( {'lst_{}'.format(sensor):out}, index=datetimes )

		# bad_date = datetime.datetime.strptime('2000-08-04', '%Y-%m-%d').strftime('%Y%j')
		# fn, = [fn for fn in files if bad_date in fn]

		df_list = df_list + [df]

	# glue the df's together
	df = df_list.pop(0)
	df = df.join( df_list )
	
	# make celcius
	df = df - 273.15

	# for year in range(2001,2017):
	# 	df_sel = df.copy().loc[slice('{}-01-01'.format(str(year)),'{}-12-31'.format(str(year)))].dropna(axis=1, how='all').copy()
	# 	ax = df_sel.plot(kind='line', title='Compare MODIS LST and WRF TSK\nyear:{}'.format(str(year)))
	# 	ax.set_ylabel('Temperature ($^\circ$C)')
	# 	plt.savefig('/workspace/Shared/Users/malindgren/MODIS_DATA/comparison_plots/COMPARE_LST_TSK_ChenaRiver_HUC_{}.png'.format(year))
	# 	plt.close()
	# 	plt.cla()



# Do A Series-length 2 week climatology -- this could be smarter by only looking at futures or something at a single time.
groups = ['rcp85', 'historical']
timelen = { 'rcp85':[2006,2018], 'historical':[2000,2015] }
for group in groups:
	begin,end = timelen[group]
	sub_df = df.loc[slice(str(begin),str(end))][[ i for i in df.columns if group in i or 'lst_' in i  ]].copy()
	twoweek_mean = sub_df.resample('2W').mean()
	weeknums = twoweek_mean.index.map( lambda x: int(x.strftime('%W')) )
	twoweek_clim = twoweek_mean.groupby(weeknums).mean()
	
	# get just the mean columns
	cols = [ i for i in twoweek_clim.columns if 'mean' in i or 'lst_' in i ]
	ax = twoweek_clim[cols].plot(kind='line', title='Compare MODIS LST and WRF TSK\nAveraged Across 2-week Intervals and All Years ({}-{})'.format( str(begin), str(end) ) )

	# # # # pull apart some stuff here for use in fill_between semantics
	models = {'rcp85':['GFDL-CM3', 'NCAR-CCSM4'], 'historical':['ERA-Interim']}

	for model in models[group]:
		# some column name -fu 
		cols = [ i for i in twoweek_clim.columns if model in i ]
		mincol, = [ i for i in cols if 'min' in i ]
		maxcol, = [ i for i in cols if 'max' in i ]
		ax.fill_between(np.array(twoweek_clim.index), twoweek_clim[mincol], twoweek_clim[maxcol], alpha=0.3)

	# # # # END fill between step

	ax.set_ylabel('Temperature ($^\circ$C)')
	ax.set_xlabel('Week of the Year')
	plt.savefig('/workspace/Shared/Users/malindgren/MODIS_DATA/comparison_plots/COMPARE_LST_TSK_ChenaRiver_HUC_twoweek_climatology_{}_{}-{}.png'.format(group, str(begin), str(end)))
	plt.close()
	plt.cla()

