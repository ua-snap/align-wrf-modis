# plot the timeseries...
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import os, glob, itertools
import pandas as pd

base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
path = os.path.join(base_path,'STATION_DATA')
out_path = os.path.join(base_path,'outputs','plots')
scenarios = ['rcp85','historical']
variables = ['tmax', 'tavg', 'tmin']

for scenario,variable in itertools.product(scenarios, variables):
	snotel
	acis_files = glob.glob(os.path.join(path, '{}*{}*.csv'.format(variable, 'acis')))
	lst_files = glob.glob(os.path.join(path, 'MOD11A2_lst_*.csv'))
	files = glob.glob( os.path.join(path, '{}*{}*.csv'.format(variable, scenario)))
	files = files + acis_files + lst_files
	for station in ['AKFAIRBANKSAP2','AKFTKNOXMINE','AKUNIVERSITYEXPSTN','AKFAIRBANKSINTLAP']:
		out_fn = os.path.join( out_path, '{}_{}_{}.png'.format(variable, scenario, station) )
		hold = []
		for fn in files:
			colname = os.path.basename(fn).split('_')[1]
			df = pd.read_csv(fn, index_col=0)
			if 'acis' in os.path.basename(fn):
				df.index = pd.DatetimeIndex(df.index)
			else:
				df.index = [pd.Timestamp.strptime(str(i),'%Y%j') for i in df.index ]
			tmp = df[station]
			tmp.name = colname
			hold = hold + [df[station].to_frame()]

		bdf = hold.pop(0)
		for i in hold:
			bdf = bdf.join(i)

		out_df = bdf.resample('M').mean()
		title = 'Compare WRF / MOD11A2 / Chena River HUC Stations\n8 Day {} averages\n{} - {}'.format(variable.upper(),station,scenario)
		colors = {'lst': (197,57,50), 'GFDL-CM3':(57,118,175), 'NCAR-CCSM4':(239,133,54),'acis':(82,157,62), 'ERA-Interim':(239,133,54) }
		colors_list = ['#{:02x}{:02x}{:02x}'.format( *colors[col] ) for col in out_df.columns]
		ax = out_df.plot(kind='line', title=title, color=colors_list )
		ax.set_ylabel('Temperatures (\u00B0C)')
		ax.set_xlabel('Time\n(monthly avg)')
		plt.tight_layout()
		plt.savefig( out_fn, figsize=(16,9), dpi=300 )
		plt.close()
		plt.cla()

