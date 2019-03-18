# make a plot dataframe for comparing the snotel, acis, lst, tsk, t2 data
# across the chena river huc polygon


if __name__ == '__main__':
	import os, glob
	import pandas as pd

	# setup pathing
	groups = ['snotel','acis']
	for group in groups:
		base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
		output_path = os.path.join( base_path,'outputs','plots' )

		# list the data
		extracted_path = os.path.join(base_path, 'data_extractions','chena_river_huc_stations_extracted',group,'chena_river_huc_station_avg')
		tsk_mod_era = glob.glob( os.path.join(extracted_path, '*{}*{}*{}*{}*'.format( 'tsk','ERA-Interim_historical','MOD11A2', group))) # get the longest series
		t2_mod_era = glob.glob( os.path.join(extracted_path, '*{}*{}*{}*{}*'.format( 't2','ERA-Interim_historical','MOD11A2', group))) # get the longest series
		station_path = os.path.join( base_path,'STATION_DATA',group,'chena_river_huc_station_avg' )
		station = glob.glob( os.path.join(station_path, '*{}*'.format( '8day')) )
		lst_mod, = glob.glob( os.path.join(extracted_path, '*{}*{}*'.format('lst','MOD11A2')))
		lst_myd, = glob.glob( os.path.join(extracted_path, '*{}*{}*'.format('lst','MYD11A2')))

		# read in the data to dict then stack to a single dataframe with the same dates
		# load the min/max/mean stations
		out = []
		for metric in ['tmin','tmax','tavg']:
			dat, = [pd.read_csv(fn, index_col=0) for fn in station if metric in fn]
			new_metric = {'tmin':'min','tmax':'max','tavg':'mean'}[metric]
			dat.columns = [group+' '+new_metric]
			out = out + [dat]
		observed = pd.concat(out, axis=1)
		del out

		# load the min/max/mean tsk
		out = []
		for metric in ['min','max','mean']:
			tsk, = [pd.read_csv(fn, index_col=0) for fn in tsk_mod_era if metric in fn]
			tsk.columns = ['tsk'+' '+metric]
			out = out + [tsk.copy(deep=True)]
		tsk = pd.concat(out, axis=1)
		del out
		
		# load the min/max/mean t2
		out = []
		for metric in ['min','max','mean']:
			t2, = [pd.read_csv(fn, index_col=0) for fn in t2_mod_era if metric in fn]
			t2.columns = ['t2'+' '+metric]
			out = out + [t2.copy(deep=True)]
		t2 = pd.concat(out, axis=1)
		del out

		# load in the mod/myd
		mod = pd.read_csv(lst_mod, index_col=0)
		mod.columns = ['MOD11A2']
		myd = pd.read_csv(lst_myd, index_col=0)
		myd.columns = ['MYD11A2']

		final = pd.concat([observed,tsk,t2,mod,myd], axis=1)
		final.index = pd.DatetimeIndex(final.index)
		out_fn = os.path.join( output_path, 'modis_wrf_{}_comparison_chena_river_huc.csv'.format(group) )
		final.to_csv( out_fn )

		# # PLIOTTING TO REMOVE
		# monfinal = final.resample('m').mean()
		# monfinal.plot()
		# plt.savefig('/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/TEST_ALLL_PLOT_month.png')
		# plt.close()

