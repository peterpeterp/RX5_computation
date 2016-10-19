import os,glob,sys,gc,time
import numpy as np 
import pandas as pd
from netCDF4 import Dataset,netcdftime,num2date
import collections

def dry_period_identifier(ind):
	# not straight forward but faster
	su=np.cumsum(ind)
	counter=collections.Counter(su)

	dry=ind*0
	index=0
	if ind[0]==0 and ind[1]==1:
		dry[0]=1
	for count,val in zip(counter.values(),counter.keys()):
		index+=count
		if count>1:
			dry[index-(count-1)/2-1]=count-1

	return(dry)

def basic_and_understandable(ind):
	dry=ind*0
	state,count=0,0
	for i in range(len(ind)):
		if ind[i]==state:
			count+=1
		if ind[i]!=state:
			dry[i-count/2-1]=((-1)**state)*count
			count=0
			if ind[i]!=99:
				state+=1
				if state==2:
					state=0
				count=1

	# still an issue with last spell??
	if state==1:	dry[i-count/2]=((-1)**1)*count
	if state==0:	dry[i-count/2]=((-1)**0)*count

	return(dry)

def test_dry_spell_identifier(N):
	ind=np.random.random(N)
	ind[ind<0.5]=0
	ind[ind>=0.5]=1
	ind=np.array(ind,'i')
	print(ind)

	start_time = time.time()
	print(basic_and_understandable(ind))
	print("--- basic_and_understandable %s seconds ---" % (time.time() - start_time))

	start_time = time.time()
	print(dry_period_identifier(ind))
	print("--- dry_period_identifier %s seconds ---" % (time.time() - start_time))




nc_in=Dataset('data/raw/pr_bced_1960_1999_hadgem2-es_rcp2p6_2011-2020.nc4',"r")
tmp=nc_in.variables['pr'][:,:,:]
pr=np.ma.getdata(tmp)*86400
mask=np.ma.getmask(tmp)
del tmp
gc.collect()

pr=np.array(pr,'i')
pr[pr<1]=0
pr[pr>=1]=1
pr[mask]=99

del mask
gc.collect()

per=np.zeros(pr.shape)

for y in range(360):
	print y
	for x in range(720):
		ind=pr[:,y,x]
		if len(np.where(ind==99)[0])<1000:
			state,count=0,0
			for i in range(pr.shape[0]):
				if ind[i]==state:
					count+=1
				if ind[i]!=state:
					per[i-count/2-1,y,x]=((-1)**state)*count
					count=0
					if ind!=99:
						state+=1
						if state==2:
							state=0
						count=1

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

first,index=True,0
# create or extend output array
if first:
	dry_spell=np.zeros([len(set(years))*12,360,720])
	out_time=np.zeros([len(set(years))*12])
	first=False
else:
	dry_spell=np.append(dry_spell,np.zeros([len(set(years))*12,360,720]))
	out_time=np.append(out_time,np.zeros([len(set(years))*12]))

# find monthly maximum
for yr in sorted(set(years)):
	for mth in sorted(set(months)):
		days_in_month=np.where((years==yr) & (months==mth))[0]
		dry_spell[index,:,:]=np.max(per[days_in_month,:,:],axis=0)
		out_time[index]=time[max(days_in_month)]
		index+=1


# copy netcdf and write zoomed file
out_file='data/raw/dry_spell_bced_1960_1999_hadgem2-es_rcp2p6_2011-2020.nc4'
os.system("rm "+out_file)
nc_out=Dataset(out_file,"w")
for dname, the_dim in nc_in.dimensions.iteritems():
	if dname=='time':nc_out.createDimension(dname, len(out_time))
	else:nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy variables
for v_name, varin in nc_in.variables.iteritems():
	if v_name=='pr':	
		outVar = nc_out.createVariable('dry_spell', varin.datatype, varin.dimensions)					    
		outVar.setncattr('long_name','maximal dry spell length with midpoint in month')
		outVar.setncattr('units','days')
		outVar[:] = dry_spell[:,:,:]
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

#np.save('data/raw/persistence_test',per)



