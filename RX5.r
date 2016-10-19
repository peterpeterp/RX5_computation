
library("RNetCDF")

movsum <- function(x,n=5){filter(x,rep(1,n), sides=1)}


compute_rx5_cmip5 <- function(model,scenario,path){
	files<-list.files(path)
	rx5<-array(NA,c(720,360,150))
	year<-array(NA,150)
	for (file in files){
		if (strsplit(file,"_")[[1]][1]=="pr" && strsplit(file,"_")[[1]][6]==scenario){
			print(file)
			nc_in<-open.nc(paste(path,file,sep=""))
			pr<-var.get.nc(nc_in,"pr")
			time<-var.get.nc(nc_in,"time")
			time<-as.numeric(strftime(as.Date(time, origin = "1860-01-01"),"%Y"))
			yrs<-unique(time)
			print(yrs)
			for (yr in yrs){
				year[yr-1949]<-yr
				rx5[,,(yr-1949)]<-apply(pr[,,which(time==yr)],c(1,2),function(x) max(movsum(x,5),na.rm=TRUE))
			}
			rm(pr)
			rm(time)
			gc()
		}
	}
	nc_out<-create.nc(paste("data/RX5/RX5_",model,"_",scenario,"_1950-2099.nc4",sep=""))

    dim.def.nc(nc_out,"lon",dimlength=720,unlim=FALSE)
    dim.def.nc(nc_out,"lat",dimlength=360,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=150,unlim=FALSE)

    var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
    att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lon",var.get.nc(nc_in,"lon"))

    var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lat",var.get.nc(nc_in,"lat"))    

    var.def.nc(nc_out,"time","NC_SHORT",c(2))
    att.put.nc(nc_out,"time", "missing_value", "NC_SHORT", -999)
    att.put.nc(nc_out, "time", "units", "NC_CHAR","year")
    var.put.nc(nc_out,"time",year) 

    var.def.nc(nc_out,"RX5","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out,"RX5", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out,"RX5", "units", "NC_CHAR","mm")
    att.put.nc(nc_out,"RX5", "long_name", "NC_CHAR","RX5")
    var.put.nc(nc_out,"RX5",rx5*86400)

    close.nc(nc_out)
}

compute_rx5_ncep <- function(path="/home/pepflei/ClimateAnalytics/data/ncep_day_pr/"){
	files<-list.files(path)
	rx5<-array(NA,c(192,94,70))
	year<-array(NA,70)
	for (file in files){
		if (strsplit(file,".sfc.")[[1]][1]=="prate"){
			print(file)
			nc_in<-open.nc(paste(path,file,sep=""))
			pr<-var.get.nc(nc_in,"prate")
			time<-var.get.nc(nc_in,"time")
			time<-as.numeric(strftime(as.POSIXct(time*3600,origin='1800-01-01 00:00'),"%Y"))
			yrs<-unique(time)
			print(yrs)
			for (yr in yrs){
				year[yr-1947]<-yr                
				rx5[,,(yr-1947)]<-apply(pr[,,which(time==yr)],c(1,2),function(x) max(movsum(x,5),na.rm=TRUE))
			}
			rm(pr)
			rm(time)
			gc()
		}
	}
	nc_out<-create.nc(paste("data/RX5/RX5_ncep_1948-2014.nc4",sep=""))

    dim.def.nc(nc_out,"lon",dimlength=192,unlim=FALSE)
    dim.def.nc(nc_out,"lat",dimlength=94,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=70,unlim=FALSE)

    var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
    att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lon",var.get.nc(nc_in,"lon"))

    var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
    att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
    var.put.nc(nc_out,"lat",var.get.nc(nc_in,"lat"))    

    var.def.nc(nc_out,"time","NC_SHORT",c(2))
    att.put.nc(nc_out,"time", "missing_value", "NC_SHORT", -999)
    att.put.nc(nc_out, "time", "units", "NC_CHAR","year")
    var.put.nc(nc_out,"time",year) 

    var.def.nc(nc_out,"RX5","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out,"RX5", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out,"RX5", "units", "NC_CHAR","mm")
    att.put.nc(nc_out,"RX5", "long_name", "NC_CHAR","RX5")
    var.put.nc(nc_out,"RX5",rx5*86400)

    close.nc(nc_out)
}

#-------------CMIP5
cmip5_parallel <- function(){
    id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    model_list=c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC-ESM-CHEM","NorESM1-M")
    path=paste("/p/projects/isimip/isimip/inputdata_bced/",model,"/",sep="")
    if (id<6){
        compute_rx5(model=model_list[(id)],scenario="historical",path=path)
    }
    if (id>5 && id<11){
        compute_rx5(model=model_list[(id-5)],scenario="rcp2p6",path=path)
    }
    if (id>10){
        compute_rx5(model=model_list[(id-10)],scenario="rcp8p5",path=path)
    }
}

cmip5_merge_hist_rcp <- function(){
    MODELS<-c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC-ESM-CHEM","NorESM1-M")
    RCPS<-c("2p6","8p5")
    models<-c("gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr","miroc-esm-chem","noresm1-m")
    rcps<-c(2.6,8.5)

    for (mod in 1:5){
        MODEL<-MODELS[mod]
        model<-models[mod]
        filename<-paste("data/RX5/CMIP5/RX5_",MODEL,"_","historical","_1950-2099.nc4",sep="") ; print(filename)
        nc_hst<-open.nc(filename)
        rx5_hst<-var.get.nc(nc_hst,"RX5")
        yr_hst<-var.get.nc(nc_hst,"time")
        for (rc in 1:2){
            RCP<-RCPS[rc]
            rcp<-rcps[rc]
            filename<-paste("data/RX5/CMIP5/RX5_",MODEL,"_rcp",RCP,"_1950-2099.nc4",sep="")   ;   print(filename)
            nc_rcp<-open.nc(filename)
            rx5_rcp<-var.get.nc(nc_rcp,"RX5")
            yr_rcp<-var.get.nc(nc_rcp,"time")

            filename<-paste("data/RX5/CMIP5/RCP",rcp,"/RX5_",model,"_rcp",rcp,"_1950-2099.nc4",sep="")  ;   print(filename)
            nc_out<-create.nc(filename)

            dim.def.nc(nc_out,"lon",dimlength=720,unlim=FALSE)
            dim.def.nc(nc_out,"lat",dimlength=360,unlim=FALSE)
            dim.def.nc(nc_out,"time",dimlength=150,unlim=FALSE)

            var.def.nc(nc_out,"lon","NC_DOUBLE",c(0))
            att.put.nc(nc_out, "lon", "units", "NC_CHAR","deg")
            var.put.nc(nc_out,"lon",var.get.nc(nc_hst,"lon"))

            var.def.nc(nc_out,"lat","NC_DOUBLE",c(1))
            att.put.nc(nc_out, "lat", "units", "NC_CHAR","deg")
            var.put.nc(nc_out,"lat",var.get.nc(nc_hst,"lat"))    

            time<-1950:2099

            rx5<-array(NA,c(720,360,150))

            rx5[,,which(time %in% yr_hst)]=rx5_hst[,,which(time %in% yr_hst)]
            rx5[,,which(time %in% yr_rcp)]=rx5_rcp[,,which(time %in% yr_rcp)]

            var.def.nc(nc_out,"time","NC_SHORT",c(2))
            att.put.nc(nc_out,"time", "missing_value", "NC_SHORT", -999)
            att.put.nc(nc_out, "time", "units", "NC_CHAR","year")
            var.put.nc(nc_out,"time",time) 

            var.def.nc(nc_out,"RX5","NC_DOUBLE",c(0,1,2))
            att.put.nc(nc_out,"RX5", "missing_value", "NC_DOUBLE", -99999.9)
            att.put.nc(nc_out,"RX5", "units", "NC_CHAR","mm")
            att.put.nc(nc_out,"RX5", "long_name", "NC_CHAR","RX5")
            var.put.nc(nc_out,"RX5",rx5)

            close.nc(nc_out)

        }
    }
   
}

#compute_rx5(model="HadGEM2-ES",scenario="historical")
#cluster_parallel()
#cmip5_merge_hist_rcp()
#-------------CMIP5


#-------------NCEP
ncep_parallel <- function(){
    id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    print(id)
    if (id==1){compute_rx5_ncep()}
}
ncep_parallel()
#-------------NCEP

#compute_rx5_ncep()






