
library("RNetCDF")

movsum <- function(x,n=5){filter(x,rep(1,n), sides=1)}


compute_rx5_cmip5 <- function(model,scenario,path){
    print(paste(model,scenario,path))
	files<-list.files(path)
	rx5<-array(NA,c(720,360,1800))
	time_out<-array(NA,1800)
	for (file in files){
		if (strsplit(file,"_")[[1]][1]=="pr" && strsplit(file,"_")[[1]][6]==scenario){
			print(file)
			nc_in<-open.nc(paste(path,file,sep=""))
			pr<-var.get.nc(nc_in,"pr")
			time<-var.get.nc(nc_in,"time")
            years<-as.numeric(strftime(as.Date(time, origin = "1860-01-01"),"%Y"))
            months<-as.numeric(strftime(as.Date(time, origin = "1860-01-01"),"%m"))
			yrs<-unique(years)
            mths<-unique(months)
			for (yr in yrs){
                for (mth in mths){
                    index<-(yr-1950)*12+mth
                    selected_year<-years==yr
                    selected_month<-months==mth
                    selected_days<-which(selected_year+selected_month==2)
                    time_out[index]<-time[max(selected_days)]
                    rx5[,,index]<-apply(pr[,,selected_days],c(1,2),function(x) max(movsum(x,5),na.rm=TRUE))
                }
			}
			rm(pr)
			rm(time)
			gc()
		}
	}
	nc_out<-create.nc(paste("data/mon_rx5/mon_rx5_",model,"_",scenario,"_1950-2099.nc4",sep=""))

    dim.def.nc(nc_out,"lon",dimlength=720,unlim=FALSE)
    dim.def.nc(nc_out,"lat",dimlength=360,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=1800,unlim=FALSE)

    var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
    att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lon",var.get.nc(nc_in,"lon"))

    var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lat",var.get.nc(nc_in,"lat"))    

    var.def.nc(nc_out,"time","NC_DOUBLE",c(2))
    att.put.nc(nc_out,"time", "missing_value", "NC_DOUBLE", -999)
    att.put.nc(nc_out, "time", "units", "NC_CHAR","days since 1860-1-1 00:00:00")
    var.put.nc(nc_out,"time",time_out) 

    var.def.nc(nc_out,"mon_rx5","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out,"mon_rx5", "missing_value", "NC_DOUBLE", -99.9)
    att.put.nc(nc_out,"mon_rx5", "units", "NC_CHAR","mm")
    att.put.nc(nc_out,"mon_rx5", "long_name", "NC_CHAR","RX5")
    var.put.nc(nc_out,"mon_rx5",rx5*86400)  #unit conversion 1 kg/m2/s = 86400 mm/day

    close.nc(nc_out)
}

compute_rx5_ncep <- function(path="/home/pepflei/ClimateAnalytics/data/ncep_day_pr/"){
	files<-list.files(path)
	rx5<-array(NA,c(192,94,804))
	time_out<-array(NA,804)
	for (file in files){
		if (strsplit(file,".sfc.")[[1]][1]=="prate"){
            print(file)
            nc_in<-open.nc(paste(path,file,sep=""))
            pr<-var.get.nc(nc_in,"prate")
            time<-var.get.nc(nc_in,"time")
            years<-as.numeric(strftime(as.Date(time/24, origin = "1800-01-01 00:00"),"%Y"))
            months<-as.numeric(strftime(as.Date(time/24, origin = "1800-01-01 00:00"),"%m"))
            yrs<-unique(years)
            mths<-unique(months)
            for (yr in yrs){
                for (mth in mths){
                    index<-(yr-1948)*12+mth
                    selected_year<-years==yr
                    selected_month<-months==mth
                    selected_days<-which(selected_year+selected_month==2)
                    time_out[index]<-time[max(selected_days)]
                    rx5[,,index]<-apply(pr[,,selected_days],c(1,2),function(x) max(movsum(x,5),na.rm=TRUE))
                }
            }
			rm(pr)
			rm(time)
			gc()
		}
	}
	nc_out<-create.nc(paste("data/mon_rx5/NCEP/mon_rx5_ncep_1948-2014.nc4",sep=""))

    dim.def.nc(nc_out,"lon",dimlength=192,unlim=FALSE)
    dim.def.nc(nc_out,"lat",dimlength=94,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=804,unlim=FALSE)

    var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
    att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lon",var.get.nc(nc_in,"lon"))

    var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lat",var.get.nc(nc_in,"lat"))    

    print(time_out)

    var.def.nc(nc_out,"time","NC_DOUBLE",c(2))
    att.put.nc(nc_out,"time", "missing_value", "NC_DOUBLE", -999)
    att.put.nc(nc_out, "time", "units", "NC_CHAR","hours since 1800-1-1 00:00:00")
    var.put.nc(nc_out,"time",time_out) 

    print(dim(rx5))

    var.def.nc(nc_out,"mon_rx5","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out,"mon_rx5", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out,"mon_rx5", "units", "NC_CHAR","mm")
    att.put.nc(nc_out,"mon_rx5", "long_name", "NC_CHAR","RX5")
    var.put.nc(nc_out,"mon_rx5",rx5*86400)  #unit conversion 1 kg/m2/s = 86400 mm/day

    close.nc(nc_out)
}

#-------------CMIP5

cmip5_merge_hist_rcp <- function(){
    MODELS<-c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC-ESM-CHEM","NorESM1-M")
    RCPS<-c("2p6","8p5")
    models<-c("gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr","miroc-esm-chem","noresm1-m")
    rcps<-c(2.6,8.5)

    for (mod in 1:5){
        MODEL<-MODELS[mod]
        model<-models[mod]
        filename<-paste("data/mon_rx5/CMIP5/mon_rx5_",MODEL,"_","historical","_1950-2099.nc4",sep="") ; print(filename)
        nc_hst<-open.nc(filename)
        rx5_hst<-var.get.nc(nc_hst,"mon_rx5")
        yr_hst<-var.get.nc(nc_hst,"time")
        for (rc in 1:2){
            RCP<-RCPS[rc]
            rcp<-rcps[rc]
            filename<-paste("data/mon_rx5/CMIP5/mon_rx5_",MODEL,"_rcp",RCP,"_1950-2099.nc4",sep="")   ;   print(filename)
            nc_rcp<-open.nc(filename)
            rx5_rcp<-var.get.nc(nc_rcp,"mon_rx5")
            yr_rcp<-var.get.nc(nc_rcp,"time")

            filename<-paste("data/mon_rx5/CMIP5/RCP",rcp,"/mon_rx5_",model,"_rcp",rcp,"_1950-2099.nc4",sep="")  ;   print(filename)
            nc_out<-create.nc(filename)

            dim.def.nc(nc_out,"lon",dimlength=720,unlim=FALSE)
            dim.def.nc(nc_out,"lat",dimlength=360,unlim=FALSE)
            dim.def.nc(nc_out,"time",dimlength=1800,unlim=FALSE)

            var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
            att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
            var.put.nc(nc_out,"lon",var.get.nc(nc_hst,"lon"))

            var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
            att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
            var.put.nc(nc_out,"lat",var.get.nc(nc_hst,"lat"))    

            time_hst<-var.get.nc(nc_hst,"time")
            time_rcp<-var.get.nc(nc_rcp,"time")
            time_out<-time_hst
            time_out[!is.na(time_rcp)]=time_rcp[!is.na(time_rcp)]

            rx5<-array(NA,c(720,360,1800))

            rx5[,,which(!is.na(time_hst))]=rx5_hst[,,which(!is.na(time_hst))]
            rx5[,,which(!is.na(time_rcp))]=rx5_rcp[,,which(!is.na(time_rcp))]

            var.def.nc(nc_out,"time","NC_DOUBLE",c(2))
            att.put.nc(nc_out,"time", "missing_value", "NC_DOUBLE", -999)
            att.put.nc(nc_out, "time", "units", "NC_CHAR","days since 1860-1-1 00:00:00")
            var.put.nc(nc_out,"time",time_out) 

            var.def.nc(nc_out,"mon_rx5","NC_DOUBLE",c(0,1,2))
            att.put.nc(nc_out,"mon_rx5", "missing_value", "NC_DOUBLE", -99.9)
            att.put.nc(nc_out,"mon_rx5", "units", "NC_CHAR","mm")
            att.put.nc(nc_out,"mon_rx5", "long_name", "NC_CHAR","mon_rx5")
            var.put.nc(nc_out,"mon_rx5",rx5)

            close.nc(nc_out)

        }
    }
   
}

cmip5_parallel <- function(){
    id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    model_list=c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC-ESM-CHEM","NorESM1-M")  
    if (id<6){
        compute_rx5_cmip5(model=model_list[(id)],scenario="historical",path=paste("/p/projects/isimip/isimip/inputdata_bced/",model_list[(id)],"/",sep=""))
    }
    if (id>5 && id<11){
        asdas
        #compute_rx5_cmip5(model=model_list[(id-5)],scenario="rcp2p6",path=paste("/p/projects/isimip/isimip/inputdata_bced/",model_list[(id-5)],"/",sep=""))
    }
    if (id>10){
        adasd
        #compute_rx5_cmip5(model=model_list[(id-10)],scenario="rcp8p5",path=paste("/p/projects/isimip/isimip/inputdata_bced/",model_list[(id-10)],"/",sep=""))
    }
}

#compute_rx5_cmip5(model="HadGEM2-ES",scenario="rcp2p6",path=paste("/p/projects/isimip/isimip/inputdata_bced/","HadGEM2-ES","/",sep=""))
#cmip5_parallel()
cmip5_merge_hist_rcp()
#-------------CMIP5


#-------------NCEP
ncep_parallel <- function(){
    id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    print(id)
    if (id==1){compute_rx5_ncep()}
}
#ncep_parallel()
#-------------NCEP

#compute_rx5_ncep()






