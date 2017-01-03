import os,glob,sys,gc
import numpy as np 
import pandas as pd
from netCDF4 import Dataset,netcdftime,num2date

# checked whether the rx5 compute with the r script is correct

#all_files=glob.glob('data/raw/pr_bced_1960_1999_hadgem2-es_rcp2p6_2011-2020.*')
models=["GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC-ESM-CHEM","NorESM1-M"]

job_id=1
print 'job_id=',job_id

for scenario in ["rcp2p6",'historical']:
	all_files=glob.glob('/p/projects/isimip/isimip/inputdata_bced/'+models[job_id]+'/pr_bced_*'+models[job_id].lower()+'*'+scenario+'*')
	all_files.sort()
	first=True
	index=0

	print all_files

	for file in all_files:
		print file
		nc_in=Dataset(file,"r")

		# extract time information
		time=nc_in.variables['time'][:]
		time_unit=nc_in.variables['time'].units
		datevar = []
		try:	# check if there is calendar information
			cal_temps = nc_in.variables['time'].calendar
			datevar.append(num2date(time,units = time_unit,calendar = cal_temps))
		except:
			datevar.append(num2date(time,units = time_unit))

		years=np.array([int(str(date).split("-")[0]) for date in datevar[0][:]])
		months=np.array([int(str(date).split("-")[1]) for date in datevar[0][:]])

		# create or extend output array
		if first:
			rx5=np.zeros([len(set(years))*12,360,720])
			out_time=np.zeros([len(set(years))*12])
			first=False
		else:
			rx5=np.concatenate((rx5,np.zeros([len(set(years))*12,360,720])),axis=0)
			out_time=np.concatenate((out_time,np.zeros([len(set(years))*12])),axis=0)

		pr=nc_in.variables['pr'][:,:,:]

		# get moving sum
		movsum=np.cumsum(pr,axis=0)
		movsum[5:]=movsum[5:]-movsum[:-5]

		# find monthly maximum
		for yr in sorted(set(years)):
			for mth in sorted(set(months)):
				days_in_month=np.where((years==yr) & (months==mth))[0]
				# the 4 first values of each month are ignored 
				# thats what happened in the R version
				rx5[index,:,:]=np.max(movsum[days_in_month[5:],:,:],axis=0)
				out_time[index]=time[max(days_in_month)]
				index+=1
		
		# problem with nearly ocean cells
		rx5[np.where(rx5>999999)]=np.nan

		# free memory... working?
		del pr
		gc.collect()



	out_file='/p/projects/tumble/carls/shared_folder/rx5/raw/mon_rx5_1960-1999_'+models[job_id]+'_'+scenario+'_1950-2099_reproduce.nc4'
	os.system("rm "+out_file)
	nc_out=Dataset(out_file,"w")
	for dname, the_dim in nc_in.dimensions.iteritems():
		if dname=='time':nc_out.createDimension(dname, len(out_time))
		else:nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

	# Copy variables
	for v_name, varin in nc_in.variables.iteritems():
		if v_name=='pr':	
			outVar = nc_out.createVariable('rx5', varin.datatype, varin.dimensions,fill_value='NaN')					    
			outVar.setncattr('long_name','monthly maximal 5 day cumulative precipitation')
			outVar.setncattr('units','mm')
			outVar.setncattr('comment','includes all types (rain, snow, large-scale, convective, etc.)')
			outVar[:] = rx5[:,:,:]*86400
		elif v_name=='time':	
			outVar = nc_out.createVariable('time', varin.datatype, varin.dimensions)					    
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = out_time[:]
		else:	
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)					    
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]

	# close the output file
	nc_out.close()
	nc_in.close()


# merge historic rcp
import os,glob,sys,gc
model="HadGEM2-ES"
for scenario in ['rcp2p6']:
	os.system('cdo mergetime /p/projects/tumble/carls/shared_folder/rx5/raw/mon_rx5_1960-1999_'+model+'_historical_1950-2099_reproduce.nc4 '+'/p/projects/tumble/carls/shared_folder/rx5/raw/mon_rx5_1960-1999_'+model+'_'+scenario+'_1950-2099_reproduce.nc4 '+'/p/projects/tumble/carls/shared_folder/rx5/mon_rx5_1960-1999_'+model.lower()+'_'+scenario+'_1950-2099_reproduce.nc4')




