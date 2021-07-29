#------------------------------------------------------------------------------------------------
# CO2 concentration simulation based on WRF-STILT/HYSPLIT-NAM12 footprints and
# SiB3/4 biosphere fluxes, CarbonTracker background fluxes and boundary conditions
# 
# created by W.He on Feb 20, 2015 
#------------------------------------------------------------------------------------------------ 
source("/Users/wei/co2simu/Rcode/lpdm/rsource/load.ncdf4.r")
source("/Users/wei/co2simu/Rcode/empirical_boundary/src/id2info.r")
source("/Users/wei/co2simu/Rcode/empirical_boundary/src/julian.r")
source("/Users/wei/co2simu/Rcode/rsource/time.r")
source("/Users/wei/co2simu/Rcode/rsource/ddate.r")
source("/Users/wei/co2simu/Rcode/rsource/datelist.r")
source("/Users/wei/co2simu/Rcode/stilt-cvs-20090417/stiltR/sourceall.r")

source("/Users/wei/co2simu/Rcode/stilt-cvs-20090417/stiltR/Trajecfoot.r")
source("/Users/wei/co2simu/Rcode/stilt-cvs-20090417/stiltR/getgridp.r")
source("/Users/wei/co2simu/Rcode/rsource/assignr.r")

source("/Users/wei/co2simu/Rcode/output.ncdf4.v2.2.r")

args <- commandArgs(trailingOnly = TRUE)
procs= as.integer(args[1])


#inputs and outputs
stiltfile = paste("stilt_",procs,".rc", sep="")

curdir="/Storage/CO2/wei/ctdas-stilt-parallel-f-only-newcov-bugfix-updatedfoot-sumwinmdm-bcCTE/exec/da/stilt/"
fn=paste(curdir, stiltfile,sep="")
conf=as.matrix(read.table(fn,header=FALSE))

sibdir  = conf[1,3]
bgfdir  = conf[2,3]
footpdir= conf[3,3]
bounddir= conf[4,3]
samdir  = conf[5,3]
outdir  = conf[6,3]

# data choice flag
bioflux_flag = "SiBCASA"  #"SiB3", "SiBCASA", "CT_OPT"
boundary_flag= "CTE"   #"CTNA", "CTE", "EMP"

# optimization option (flux only, or flux + boundary)
opt_para_flag = "FLUX"    # "FLUX_BOU","FLUX"

#----------------------------------------------------------------------------------------------
#Convolve footprints with sib hourly fluxes, for observations from both Aircraft and towers
#----------------------------------------------------------------------------------------------
endtime=240
foottimes=seq(0,endtime,1)     #vector of times (backtimes) in hours between which footprint is computed
zbot=0                         #lower vertical bound for influence projection, in meters agl
ztop=0                         #upper vertical bound for influence projection, in meters agl
#if ztop set to zero, *surface* influence will be calculated

#set up an equivalent domain among footprints,fluxes, but not CO2 boundary 
ncol2=66                       #number of pixels in x directions in grid 
nrow2=40                       #number of pixels in y directions in grid 
LLLon2=(-129)                  #lower left corner of grid 
LLLat2=22                      #lower left corner of grid 
ResLon2=1                      #resolution in degrees longitude 
ResLat2=1                      #resolution in degrees latitude 

# constants for level-height transfrom
levelhgt=c(34.5,111.9,256.9,490.4,826.4,1274.1,1839.0,2524.0,3329.9,4255.6,5298.5,6453.8,7715.4,9076.6,
           + 10533.3,12108.3,13874.2,15860.1,18093.2,20590.0,24247.3,29859.6,35695.0,42551.5,80000.0)
#data from http://www.esrl.noaa.gov/gmd/ccgg/carbontracker-ch4/documentation_tm5.html

TBegin <- proc.time()
library(ncdf4)
library(futile.logger)

logfile = paste("stilt_",procs,".log", sep="")
flog.appender(appender.file(logfile), name='logger.b')
flog.info("Launch WRF-STILT-based forward simulations %s", "file", name='logger.b')

nsam=as.numeric(conf[7,3])
ndayscycle=as.numeric(conf[8,3])

startdate=conf[9,3]
year=as.numeric(substring(startdate,1,4))
month=as.numeric(substring(startdate,6,7))
day=as.numeric(substring(startdate,9,10))
b0 = ISOdatetime(year,month,day,0,0,0,tz="UTC")

startdate=conf[10,3]
year=as.numeric(substring(startdate,1,4))
month=as.numeric(substring(startdate,6,7))
day=as.numeric(substring(startdate,9,10))
if(ndayscycle!=10)
  flog.fatal("Sorry, only support the case nndayscycle equals 10", name='logger.b') 

b2 = ISOdatetime(year,month,day,0,0,0,tz="UTC")
b1 = b2 - ndayscycle*24*3600      # for footprint
b3 = b1 + 2*ndayscycle*24*3600    # for flux

tb1 = paste(substring(b1,1,4),substring(b1,6,7),substring(b1,9,10),sep="")
tb2 = paste(substring(b2,1,4),substring(b2,6,7),substring(b2,9,10),sep="")
tb3 = paste(substring(b3,1,4),substring(b3,6,7),substring(b3,9,10),sep="")

scalefacarr1=array(NA, dim=c(ncol2,nrow2,nsam))   
scalefacarr2=array(NA, dim=c(ncol2,nrow2,nsam))   
# make CTDAS to generate domain for North America, read data partly?

flog.info("Reading scaling factor files", name='logger.b')
for(i in 0:(nsam-1))    #parameters.000.2010010100_2010011100.nc
{
  if (i<10)
    ii = paste("00",i,sep="")
  if (i<100 && i>=10)
    ii = paste("0",i,sep="")
  if (i>=100)
    ii = i
  
  if (b1<b0)   #b1==b0, revise on Aug 5,2015
    scalefacarr1[,,i+1] = 1
  
  if (b1>=b0)  #b1>b0
  {
    ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb1,"00","_",tb2,"00",".nc",sep="")) 
    scalefac <- ncvar_get(ncf,"parametermap",start=c(52,113),count=c(ncol2,nrow2))   #real52:117,113:152,start=c(52,113),count=c(66,40)
    scalefacarr1[,,i+1] = scalefac
    nc_close(ncf)
  }
  
  ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb2,"00","_",tb3,"00",".nc",sep="")) 
  scalefac <- ncvar_get(ncf,"parametermap",start=c(52,113),count=c(ncol2,nrow2))  
  scalefacarr2[,,i+1] = scalefac  
  nc_close(ncf)
}

#---------------------------------------------------------------------------------
# Centering at the state vector especially the ready-optimized cycle, look for data
# from corresponding period for optimzation
#---------------------------------------------------------------------------------
# according to b1, b2, decide which months
y1=substring(b2,1,4)
y2=substring(b3,1,4)
m1=substring(b2,6,7)
m2=substring(b3,6,7)
uyr=unique(c(y1,y2))
umon=unique(c(m1,m2))

flog.info("Determining ready-use footprint files", name='logger.b')
fns=NULL
for(yy in 1:length(uyr))
{
  flog.info("Footprints in year %s",uyr[yy],name='logger.b')
  for(mm in 1:length(umon))  
  {
      pfbpath=paste(footpdir,uyr[yy],"/",umon[mm],"/",sep="")  #"stilt2010x01x01x18x59x33.4057Nx081.8334Wx00285.nc"
      tmp=list.files(pfbpath,pattern="nc", all=T)
      fns=c(fns,tmp)
  }
}
# read needed event ids from sample_coordinates files
ncf <- nc_open(paste(samdir,"sample_coordinates_",tb2,"00","_",tb3,"00","-",procs,".nc",sep="")) 
eventidarr <- ncvar_get(ncf,"obs_num")

# filter before loop
pfbfns=NULL
fn=NULL
for(mm in 1:length(fns))  
{  
  
  eid=as.numeric(substring(fns[mm],48,53))
  for(i in 1:length(eventidarr))  
  { 
    if(eid==eventidarr[i])
    {
      pfbfns = c(pfbfns,fns[mm])
      break
    }
  }
  
}
nc_close(ncf)

#evid matching result
n_obs=length(eventidarr)
n_foot=length(fns)
n_used=length(pfbfns)
log=paste("Number of event ids from ctdas: ",n_obs,", Number of matched: ",n_used)
flog.info(log, name='logger.b')  #June 23

log=paste("For state vectors from ",b1,"to",b2,",",length(pfbfns),"observations have been found")
flog.info(log, name='logger.b')

#flog.info("Start convolving for all footprint files", name='logger.b')
newfile=T #output results into a new file

#
#read EMP boundary monthly into memory
#

if(boundary_flag == "EMP")
{
    bouf=paste(bounddir,y1,"/",m1,"/","CO2.v201209_v1.init.stilt.",y1,m1,".txt",sep="")
    print(bouf)
    val<-as.matrix(read.table(bouf,header=TRUE))
}

for(mm in 1:length(pfbfns))  
{
  #--------------------------------------------
  # STEP 1: Read footprints
  #--------------------------------------------
  yr=substring(pfbfns[mm],6,9)
  mo=substring(pfbfns[mm],11,12)
  fn=paste(footpdir,yr,"/",mo,"/",pfbfns[mm],sep="")
  footp=load.ncdf(fn)  
  
  footnc=nc_open(fn)
  foot=ncvar_get(footnc,"foot1",start=c(41,12,1),count=c(ncol2,nrow2,-1))
  
  # get particles
  #part=footp$partfoot 
  #partna=footp$partfootnames
  #colnames(part)=c("time","index","lat","lon","agl","grdht","foot","temp","temp0","swrad","zi","dens","pres","dmass")   #ndmass
  
  # get info form path srings
  ident=substring(pfbfns[mm],6,46)
  eventid=substring(pfbfns[mm],48,53) 
  
  print(ident)
  flog.info("Convolving for %s",ident, name='logger.b')
  info=id2info(ident) 
  
  #for different spatial resolutions between biosphere fluxes and background fluxes, two footprints are planned to be used, and to the end, we picked up the concentration values.
  #foot1=Trajecfoot(ident=ident,part=part, pathname="",foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,numpix.x=ncol2,numpix.y=nrow2,lon.ll=LLLon2,lat.ll=LLLat2,lon.res=ResLon2,lat.res=ResLat2)
  #foot1 40  66 240
  
  #foot=array(NA, dim=c(66,40,240))
  #for (i in 1:240)
  #    foot[,,i]=t(as.matrix(foot1[,,i]))  #113:152
  
  nc_close(footnc)
  
  if(length(foot)>100)
  {			
    inityr=as.numeric(substring(ident,1,4))
    initmo=as.numeric(substring(ident,6,7))
    initdy=as.numeric(substring(ident,9,10))
    inithr=as.numeric(substring(ident,12,13))
    initmi=as.numeric(substring(ident,15,16))
    inihgt=as.numeric(substring(ident,39,41))
    
    # get the time stamp for each foot step, going back time given by "endtime"
    xx=ISOdatetime(inityr,initmo,initdy,inithr,initmi,0,tz="UTC")
    yy=xx+(-foottimes*3600)
    cyy=as.character(yy)
    yrlist=substring(cyy,1,4)
    monlist=substring(cyy,6,7)
    daylist=substring(cyy,9,10)
    hrlist=substring(cyy,12,13)
    milist=substring(cyy,15,16)
    
    # get unique months and days
    daystring=paste(yrlist,monlist,daylist, sep="")
    udaystring=unique(daystring)
    yrmonstring=paste(yrlist,monlist, sep="")
    uyrmonstring=unique(yrmonstring)
    
    #current month
    sibyr=substring(uyrmonstring[1],1,4)
    sibmon=substring(uyrmonstring[1],5,6)	
    
    #----------------------------------------------------------------------------------
    # STEP 2: Read boundary conditions & use end points tracing to 
    # 1) get responding fluxes & convolve with footprints
    # 2) get the "actural" concentrations
    #----------------------------------------------------------------------------------
    # get endpoints
    endpts=footp$endpts
    endptna=footp$endptsnames 
    #endpts=ncvar_get(footnc,"endpts")
    #endptna=ncvar_get(footnc,"endptsnames")	
    colnames(endpts)=endptna   #=c("time","index","lat","lon","agl","grdht","temp","pres")  
    
    endptsdate=footp$endptsdate  #2010-08-04 20:18:00 UTC
    #endptsdate=ncvar_get(footnc,"endptsdate")
    latarr=endpts[,3]
    lonarr=endpts[,4]
    hgtarr=endpts[,5] + endpts[,6]  #agl+grdht
    
    # analyze the dates to recognize how many data to read
    dd=as.character(endptsdate)
    yrlist=substring(dd,1,4)
    monlist=substring(dd,6,7)
    daylist=substring(dd,9,10)
    hrlist=substring(dd,12,13)
    milist=substring(dd,15,16)
    
    # get unique months and days
    daystring=paste(yrlist,monlist,daylist, sep="")
    udaystring=unique(daystring)
    yrmonstring=paste(yrlist,monlist, sep="")
    uyrmonstring=unique(yrmonstring)
    
    #--------------------------------------------------------------
    # 2-1: Read fluxes and boundary data
    #--------------------------------------------------------------
    ndays=length(udaystring)
   
    if(boundary_flag=="CTNA")
        bouarr=array(0,dim=c(ndays,120,90,25,8))  #boundary : lon,lat,hgt,time
    if(boundary_flag=="CTE")
        bouarr=array(0,dim=c(ndays,360,180,25,8))
    
    #biospheric fluxes	
    
    nrow_flux1=181 #181
    ncol_flux1=360 #361 #288
    
    if(bioflux_flag == "SiB3")
    {
        gpp=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
        rec=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
        #pla=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
        #soi=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
    }

    #other fluxes #(360,180,8)	
    nrow_flux2=180
    ncol_flux2=360
    ocn=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))   
    fos=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    fir=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    
    if(bioflux_flag == "SiBCASA")
    {
        gpp=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
        rec=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    }

    if(bioflux_flag == "CT_OPT")
        bio=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    
    ntimes1=ndays*24  
    ntimes2=ndays*8  
    
    for(d in 1:ndays) 
    {
      datestr=udaystring[d]
      yr=substr(datestr,1,4)
      mn=substr(datestr,5,6)
      dy=substr(datestr,7,8)
    
      if(boundary_flag=="CTNA")
          bou=load.ncdf(paste(bounddir,"CT2013B.molefrac_glb3x2_",yr,"-",mn,"-",dy,".nc",sep=""))
          #co2(date, level, lat, lon) , ocn_flux_opt(date, lat, lon)
      if(boundary_flag=="CTE")
          bou=load.ncdf(paste(bounddir,"3d_molefractions_1x1_",yr,mn,dy,".nc",sep=""))
          
      bgf=load.ncdf(paste(bgfdir,"CT2013B.flux1x1.",yr,mn,dy,".nc",sep=""))

      if(bioflux_flag == "SiB3")
      {
          biof=load.ncdf(paste(sibdir,"SiB3.hourly.flux1x1.global.",yr,mn,dy,".nc",sep=""))  
          gpp[d,,,]=biof$gpp
          rec[d,,,]=biof$rtotal
          #pla[d,,,]=biof$ocs.gpp
          #soi[d,,,]=biof$ocs.soil
          remove(list=c("biof"))
      }

      if(bioflux_flag == "SiBCASA")
      {
          biof=load.ncdf(paste(sibdir,"SiBCASA.3hourly.flux1x1.global.",yr,mn,dy,".nc",sep=""))
          gpp[d,,,]=biof$gpp
          rec[d,,,]=biof$resp
          remove(list=c("biof"))
      }

      #flog.info("mean gpp=  %f", mean(gpp[d,,,]), name='logger.b')

 
      if(bioflux_flag == "CT_OPT")
          bio[d,,,]=bgf$bio.flux.opt

      if(boundary_flag=="CTNA")
      {
          bouarr[d,,,,]=bou$co2
          remove(list=c("bou"))
      }
      if(boundary_flag=="CTE")
      {
          bouarr[d,,,,]=bou$co2*1e6
          remove(list=c("bou"))
      }

      ocn[d,,,]=bgf$ocn.flux.opt        # here we can read less, so change the sizes of ocnarr
      fos[d,,,]=bgf$fossil.flux.imp
      fir[d,,,]=bgf$fire.flux.imp

      remove(list=c("bgf"))
    }
    
    #--------------------------------------------------------------
    # 2-2: prepare data for calculations
    #--------------------------------------------------------------
    dateall1=rep(ISOdatetime(0,0,0,0,0,0,tz="UTC"), ntimes1) 
    dateall2=rep(ISOdatetime(0,0,0,0,0,0,tz="UTC"), ntimes2)
        
    if(bioflux_flag == "SiB3")
    {
        tdim = ntimes1
        gppflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
        recflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
        #plaflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
        #soiflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    }

    if(bioflux_flag == "SiBCASA")
    {
        tdim = ntimes2
        gppflux=array(NA, dim=c(ncol2,nrow2,ntimes2))
        recflux=array(NA, dim=c(ncol2,nrow2,ntimes2))
    }
    
    if(bioflux_flag == "CT_OPT")
    {
        tdim = ntimes2
        bioflux=array(NA, dim=c(ncol2,nrow2,ntimes2))
    }
   
    ocnflux=array(NA, dim=c(ncol2,nrow2,ntimes2))  #ntimes1->ntimes2, revise on June 7,2015
    fosflux=array(NA, dim=c(ncol2,nrow2,ntimes2))
    firflux=array(NA, dim=c(ncol2,nrow2,ntimes2))

    #flog.info("tdim =  %f", tdim, name='logger.b')
    

    neeflux=array(NA, dim=c(ncol2,nrow2,tdim))
    neeflux1=array(NA, dim=c(ncol2,nrow2,tdim))
    neefluxarr=array(NA, dim=c(ncol2,nrow2,tdim,nsam))

    neeoutarr=array(NA, dim=c(nsam))
    fsimuarr=array(NA, dim=c(nsam))
   
   #----------------------------------------------------
    yr=rep(sibyr,ntimes1)
    mon=rep(sibmon,ntimes1)
    hrs=seq(1,ntimes1,1) 
    dy=ceiling(hrs/24)
    hr=(hrs-(dy-1)*24)*1-1
    time1=ISOdatetime(yr,mon,dy,hr,0,0,tz="UTC")  
    
    yr=rep(sibyr,ntimes2)
    mon=rep(sibmon,ntimes2)
    hrs=seq(1,ntimes2,1) 
    dy=ceiling(hrs/8)
    hr=(hrs-(dy-1)*8)*3-1
    time2=ISOdatetime(yr,mon,dy,hr,0,0,tz="UTC")  
    
    if(bioflux_flag == "SiB3")
    {
        for(hh in ntimes1:1)
        {
           dateall1[ntimes1-hh+1]=time1[hh]
           inxd=ceiling(hh/24)
           inxh=hh-(inxd-1)*24
           gppflux[,,ntimes1-hh+1]=gpp[inxd,52:117,113:152,inxh]  #transform matrix, 41:93; delete this on May 4,2015
           recflux[,,ntimes1-hh+1]=rec[inxd,52:117,113:152,inxh]  
           #plaflux[,,ntimes1-hh+1]=pla[inxd,52:117,113:152,inxh]
           #soiflux[,,ntimes1-hh+1]=soi[inxd,52:117,113:152,inxh]
        }
        # replace NA values with 0, as NA was used as zero values for datasets
        gppflux[,,]=replace(gppflux[,,],is.na(gppflux[,,]),0)
        recflux[,,]=replace(recflux[,,],is.na(recflux[,,]),0)
    }


    for(hh in ntimes2:1)
    {
      dateall2[ntimes2-hh+1]=time2[hh]
      inxd=ceiling(hh/8)
      inxh=hh-(inxd-1)*8
      ocnflux[,,ntimes2-hh+1]=ocn[inxd,52:117,113:152,inxh]  #transform matrix, delete this 
      fosflux[,,ntimes2-hh+1]=fos[inxd,52:117,113:152,inxh]  
      firflux[,,ntimes2-hh+1]=fir[inxd,52:117,113:152,inxh]   
     
      if(bioflux_flag == "SiBCASA")
      {
          gppflux[,,ntimes2-hh+1]=gpp[inxd,52:117,113:152,inxh]
          recflux[,,ntimes2-hh+1]=rec[inxd,52:117,113:152,inxh]
      }

      if(bioflux_flag == "CT_OPT")
          bioflux[,,ntimes2-hh+1]=bio[inxd,52:117,113:152,inxh]

    }
     # replace NA values with 0, as NA was used as zero values for datasets
     ocnflux[,,]=replace(ocnflux[,,],is.na(ocnflux[,,]),0)
     fosflux[,,]=replace(fosflux[,,],is.na(fosflux[,,]),0)
     firflux[,,]=replace(firflux[,,],is.na(firflux[,,]),0)
    
     if(bioflux_flag == "CT_OPT")
         bioflux[,,]=replace(bioflux[,,],is.na(bioflux[,,]),0)
   
     if(bioflux_flag == "SiBCASA")
     {
         gppflux[,,]=replace(gppflux[,,],is.na(gppflux[,,]),0)
         recflux[,,]=replace(recflux[,,],is.na(recflux[,,]),0)
     }
       
     if(bioflux_flag == "SiB3")
         neeflux[,,]=recflux[,,]-gppflux[,,]  # nee*1e6 for SiBCASA?

     if(bioflux_flag == "SiBCASA")
         neeflux[,,]=(recflux[,,]-gppflux[,,])*1e6

     if(bioflux_flag == "CT_OPT")
         neeflux[,,]= bioflux[,,]*1e6  # for scaling CT flux
   
    flog.info("bioflux_flag = %s", bioflux_flag, name='logger.b')
    flog.info("boundary_flag = %s", boundary_flag, name='logger.b')
    #flog.info("nee flux, max=%f, min=%f, mean=%f", max(neeflux),min(neeflux),mean(neeflux), name='logger.b')

    #--------------------------------------------------------------
    # 2-3: Convolving footprints with fluxes to get delta co2
    #--------------------------------------------------------------
    
    #**********************************
    # biosphere fluxes convolving
    #**********************************
    # sib3 flux, 1-hourly
    if(bioflux_flag == "SiB3")
    {
        datevalid=dateall1[1]+0.5*3600-(0:(ntimes1-1))*3600 	
        ixflux=(1:ntimes1)[(datevalid<=xx)][1:endtime]   #time backward
        
        #ocstmp=foot*plaflux[,,ixflux]
        #soiltmp=foot*soiflux[,,ixflux]
        gpptmp=foot*gppflux[,,ixflux]
        recotmp=foot*recflux[,,ixflux]
    
        #ocsout=sum(ocstmp,na.rm=TRUE)
        #soilout=sum(soiltmp,na.rm=TRUE)
        gppout=sum(gpptmp,na.rm=TRUE)
        recoout=sum(recotmp,na.rm=TRUE)

        remove(list=c("gppflux","recflux","gpptmp","recotmp","soiltmp","dateall1"))
        gc()

    }
    #dateall[ixflux], datevalid from dateall1, from time1, from SiB flux; xx from from footprint
    
    #****************************************************
    # scaling NEE, determine which fluxes need scaling
    #****************************************************
    
    # involve scaling factors to neeflux[,,ixflux], partly, for example as follow
    # xx = "2010-01-24 18:10:00 UTC" || b1 = "2010-01-19 12:00:00 UTC", b3 = "2010-01-29 12:00:00 UTC"
    
    xxleft=xx-ndays*24*3600
   
    #flog.info("xxleft = %s, b2=%s, b2-xxleft=%s", xxleft, b2, b2-xxleft,name='logger.b')
   

    if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
        flog.info("xxleft = %s, b2=%s, b2-xxleft=%s", xxleft, b2, b2-xxleft,name='logger.b')

    for (i in 1:nsam) 
    {
      if(xxleft<b2)  #flux with scaling factors located in first lag
      {    
        diff=abs(as.numeric(difftime(b2,xxleft,units='hours')))
        #flog.info("before diff = %f", diff, name='logger.b')

        if(bioflux_flag == "SiBCASA" || bioflux_flag == "CT_OPT")
            diff = diff/3.0
        diff=ceiling(diff)
        
        if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
           flog.info("after diff = %f", diff, name='logger.b')

        for (hh in 1:(diff)) 
            if(hh<=diff && hh>=1)
            {
              if(i==1 && ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")              
                 flog.info("neeflux[,,hh]=%f, scalefacarr1 = %f ",neeflux[,,hh], scalefacarr1[,,i], name='logger.b')
                
               neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr1[,,i])   #delete "t" transform
             }
        for (hh in (diff+1):tdim)
            if(hh<=tdim && hh>=diff+1)
            {
                if(i==1 && ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
                   flog.info("neeflux[,,hh]=%f, scalefacarr2 = %f",neeflux[,,hh], scalefacarr2[,,i],name='logger.b')
                neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr2[,,i])
            }
        neefluxarr[,,,i] = neeflux1[,,]  
      }      
      else  #all scaling within second lag
      {
        for (hh in 1:tdim) 
          neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr2[,,i])    
        
        neefluxarr[,,,i] = neeflux1[,,]   
      }
    }
   
    if(bioflux_flag == "SiBCASA" || bioflux_flag == "CT_OPT") 
    {
       datevalid=dateall2[1]+1.5*3600-(0:(ntimes2-1))*3*3600
       ixflux=(1:ntimes2)[(datevalid<=xx)][1:endtime]
    }

    if(bioflux_flag == "SiBCASA")
    {
       gpptmp=foot*gppflux[,,ixflux]*1e6
       recotmp=foot*recflux[,,ixflux]*1e6

       gppout=sum(gpptmp,na.rm=TRUE)
       recoout=sum(recotmp,na.rm=TRUE)

       remove(list=c("gpptmp","recotmp"))
       gc()
     }

     # delta co2 on NEE for all ensemble members
     for(i in 1:nsam)
     {
        neetmp=foot*neefluxarr[,,ixflux,i]   
        if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
            flog.info("neetmp = %f",neetmp, name='logger.b')

        neeout=sum(neetmp,na.rm=TRUE)
        neeoutarr[i]=neeout
      
        if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
            flog.info("neeoutarr[i] = %f",neeoutarr[i], name='logger.b')
    }

    remove(list=c("neetmp","neeout"))
    gc()
    
    #**********************************
    # component fluxes convolving
    #**********************************
    # bg fluxes, 3-hourly
    datevalid=dateall2[1]+1.5*3600-(0:(ntimes2-1))*3*3600	
    ixflux=(1:ntimes2)[(datevalid<=xx)][1:endtime]
    
    #need large amount of memories
    ocntmp=foot*ocnflux[,,ixflux]*1e6
    fostmp=foot*fosflux[,,ixflux]*1e6
    firtmp=foot*firflux[,,ixflux]*1e6
    
    ocnout=sum(ocntmp,na.rm=TRUE)
    fosout=sum(fostmp,na.rm=TRUE) 
    firout=sum(firtmp,na.rm=TRUE)

    if(bioflux_flag == "CT_OPT")
    {
        biotmp=foot*bioflux[,,ixflux]*1e6
        bioout=sum(biotmp,na.rm=TRUE)
    }
    
    remove(list=c("ocnflux","fosflux","firflux","ocntmp","fostmp","firtmp","foot","dateall2"))
    gc()
    
    #--------------------------------------------------------------
    # 2-4: Calculate boundary  
    #--------------------------------------------------------------
    #case 1: data from CT2013B
    if(boundary_flag == "CTNA" || boundary_flag == "CTE")
    {
      # boundary domain : lat 22.5~61.5, lon -128.5~-63.5
      latmin=-90   #22.5
      latmax=90    #61.5
      lonmin=-180  #-128.5
      lonmax=180   #-63.5
      npts=dim(endpts)[1]   
    
      # calulate concentrations
      pbou=0
      nval=0
      dmin=min(as.integer(substr(udaystring,7,8)))
      for(p in 1:npts) 
      { 
        dy=0
        goal=as.integer(daylist[p])
        for(uu in 1:ndays) 
        {
          dstr=udaystring[uu]
          ddy=substr(dstr,7,8)
        
          if(ddy<=goal)
            dy=dy+1
        }
     
      if(boundary_flag == "CTNA")
      {
         i=ceiling((lonarr[p]-lonmin)/3.0)
         j=ceiling((latarr[p]-latmin)/2.0)
      }
      if(boundary_flag == "CTE")
      {
         i=ceiling((lonarr[p]-lonmin)/1.0)
         j=ceiling((latarr[p]-latmin)/1.0)
      }

      # height matching
      k0=hgtarr[p]    #go upward stair stategy, one by one
      k=1
      for(l in 1:25) 
      {
        if(k0 > levelhgt[l] && k0 <= levelhgt[25])
          k=l+1    
        if(k0 > levelhgt[25])
          k=25
      }
      
      # 3-hourly matching /agt/alt
      hpos=ceiling((as.integer(hrlist[p])+0.000001+as.integer(milist[p])/60.0 )/3.0)
      
      tt=bouarr[dy,i,j,k,hpos]
      
      if(length(tt)==1)   #why sometimes we get a unnormal array?
      {
        pbou=pbou+tt  #sum  
        nval=nval+1
      }
      
     }
     fbou=pbou/nval
    }
    
    #case 2: data from Arlyn's emperical curtain
    if(boundary_flag == "EMP")
    {
       #bouf=paste(bounddir,year,"/",mon,"/","CO2.v201209_v1.init.hysplit.",year,mon,".txt",sep="")
       #print(bouf) 
       #val<-as.matrix(read.table(bouf,header=TRUE))
       for(p in 1:dim(val)[1])
           if(val[p,2]==eventid)
           {
              fbou=as.numeric(val[p,3])
              break
           }
    }
    
    #################################################
    # final results and output to files
    #################################################
    
    fnameout1=paste(outdir,"/","samples_simulated.",tb2,"00","_",tb3,"00","-",procs,".txt",sep="")
    fn=fnameout1
    outlist=NULL #c(inityr, initmo, initdy, inithr, initmi, inihgt, fbou)
    for(i in 1:nsam)
    { 
      deltaco2=neeoutarr[i]+ocnout+fosout+firout 
      fsimuarr[i]=(fbou+deltaco2)*1e-6        
      
      if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
          flog.info("fsimuarr[i] = %f",fsimuarr[i], name='logger.b')

      outlist<-c(outlist, fsimuarr[i])
    }

    #add component varibles:gpp,reco,ocn,fos,fir,boundary
    #outlist<-c(outlist,gppout,recoout,ocnout,fosout,firout,fbou,bioout)
    if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
        outlist<-c(outlist,gppout*1e-6,recoout*1e-6,ocnout*1e-6,fosout*1e-6,firout*1e-6,fbou*1e-6)
    if(bioflux_flag == "CT_OPT")
        outlist<-c(outlist,bioout*1e-6,ocnout*1e-6,fosout*1e-6,firout*1e-6,fbou*1e-6)


    out1<-c(ident,eventid,round(outlist,10))  
    if(newfile) 
      write.table(t(out1), file=fnameout1, append=F, col.names=F, row.names=F, quote=F)
    if(!newfile) 
      write.table(t(out1), file=fnameout1, append=T, col.names=F, row.names=F, quote=F)
    newfile=F
    
  }# end if foot NA
  
}#end loop mm

#-----------------------------------------------------------------------
# ouput rsults as *.nc files
#-----------------------------------------------------------------------
fin=paste(fn,sep="")  #already include eventid 
fout=paste(substring(fn,1,nchar(fn)-4),".nc",sep="")
data=as.matrix(read.table(fin,header=FALSE))

nobs=dim(data)[1]

if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
{
    nmem=dim(data)[2]-8     #9
    vals=data[,2:(nsam+8)]  #151+8 colums
}

if(bioflux_flag == "CT_OPT")
{
    nmem=dim(data)[2]-7     #9
    vals=data[,2:(nsam+7)]  #151+8 colums
}
# use a flag to specify for different flux simulation (SiB3 or CT_OPT)
write.results.netcdf(bioflux_flag,vals,nobs,nmem,fout)   #write simulated results once
#write.results.netcdf(vals,nobs,nmem,fout)

#create stilt.ok file
file.create(paste("stilt_",procs,".ok",sep=""))

#---------------------------------------------------------------------------------------
# calculate 10-days mean SiB3 biospheric fluxes and other background fluxes
#---------------------------------------------------------------------------------------
if (procs==0)
{

flog.info("Calculate mean background fluxes for each cycle 10 days", name='logger.b')

newfile=T #output results into a new file
ndays=10

#biospheric fluxes  
nrow=180
ncol=360

if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
{
    gpp=array(NA, dim=c(ncol,nrow,ndays))
    res=array(NA, dim=c(ncol,nrow,ndays))
}

if(bioflux_flag == "CT_OPT")
   bio=array(NA, dim=c(ncol,nrow,ndays))

ocn=array(NA, dim=c(ncol,nrow,ndays))
fos=array(NA, dim=c(ncol,nrow,ndays))
fir=array(NA, dim=c(ncol,nrow,ndays))

start=b2
end=b2+ndays*24*3600
str=paste("for fluxes from",start,"to",end,sep=" ")
flog.info(str, name='logger.b')

datelist = xx + (0:(ndays-1))*24*3600
#xx=xx+ndays*24*3600

for(d in 1:ndays) 
{
  datestr=datelist[d]
  yr=substr(datestr,1,4)
  mn=substr(datestr,6,7)
  dy=substr(datestr,9,10)
  
  if(bioflux_flag == "SiB3")
  {
      biof=load.ncdf(paste(sibdir,"SiB3.hourly.flux1x1.global.",yr,mn,dy,".nc",sep=""))  
      biofgpp=replace(biof$gpp, is.na(biof$gpp),0)                    # is.nan
      biofres=replace(biof$rtotal, is.na(biof$rtotal),0)
      gpp[,,d]=rowMeans(biofgpp[,1:180,], na.rm = FALSE, dims = 2)    # notice true or false
      res[,,d]=rowMeans(biofres[,1:180,], na.rm = FALSE, dims = 2)
  }

 if(bioflux_flag == "SiBCASA")
 {
      biof=load.ncdf(paste(sibdir,"SiBCASA.3hourly.flux1x1.global.",yr,mn,dy,".nc",sep=""))
      biofgpp=replace(biof$gpp, is.na(biof$gpp),0)                    # is.nan
      biofres=replace(biof$resp, is.na(biof$resp),0)
      gpp[,,d]=rowMeans(biofgpp[,1:180,], na.rm = FALSE, dims = 2)    # notice true or false
      res[,,d]=rowMeans(biofres[,1:180,], na.rm = FALSE, dims = 2)
  }

  bgf=load.ncdf(paste(bgfdir,"CT2013B.flux1x1.",yr,mn,dy,".nc",sep=""))
 
  bgfocn=replace(bgf$ocn.flux.opt, is.na(bgf$ocn.flux.opt),0)
  bgffos=replace(bgf$fossil.flux.imp, is.na(bgf$fossil.flux.imp),0)
  bgffir=replace(bgf$fire.flux.imp, is.na(bgf$fire.flux.imp),0)

  ocn[,,d]=rowMeans(bgfocn, na.rm = FALSE, dims = 2)
  fos[,,d]=rowMeans(bgffos, na.rm = FALSE, dims = 2)
  fir[,,d]=rowMeans(bgffir, na.rm = FALSE, dims = 2)

  if(bioflux_flag == "CT_OPT")
  {
      bgfbio=replace(bgf$bio.flux.opt, is.na(bgf$bio.flux.opt),0)
      bio[,,d]=rowMeans(bgfbio, na.rm = FALSE, dims = 2)
  }
 # gpp[,,d]=rowMeans(biof$gpp[,1:180,], na.rm = TRUE, dims = 2)    # notice true or false
 # res[,,d]=rowMeans(biof$rtotal[,1:180,], na.rm = TRUE, dims = 2)
 # ocn[,,d]=rowMeans(bgf$ocn.flux.opt, na.rm = TRUE, dims = 2)
 # fos[,,d]=rowMeans(bgf$fossil.flux.imp, na.rm = TRUE, dims = 2)
 # fir[,,d]=rowMeans(bgf$fire.flux.imp, na.rm = TRUE, dims = 2)
}

# calculate the mean values for these fluxes,mean by row/col
if(bioflux_flag == "SiB3")
{
   gppmean_hour=rowMeans(gpp, na.rm = FALSE, dims = 2)*1e-6
   resmean_hour=rowMeans(res, na.rm = FALSE, dims = 2)*1e-6
}

if(bioflux_flag == "SiBCASA")
{
   gppmean_hour=rowMeans(gpp, na.rm = FALSE, dims = 2)
   resmean_hour=rowMeans(res, na.rm = FALSE, dims = 2)
}

if(bioflux_flag == "CT_OPT")   #save the same NEE values for GPP and Re
{
   gppmean_hour=rowMeans(bio, na.rm = FALSE, dims = 2)
   resmean_hour=gppmean_hour
}

ocnmean_hour=rowMeans(ocn, na.rm = FALSE, dims = 2)
fosmean_hour=rowMeans(fos, na.rm = FALSE, dims = 2)
firmean_hour=rowMeans(fir, na.rm = FALSE, dims = 2)

# gppmean_hour=replace(gppmean_hour,is.nan(gppmean_hour),0)
# resmean_hour=replace(resmean_hour,is.nan(resmean_hour),0)
# ocnmean_hour=replace(ocnmean_hour,is.nan(ocnmean_hour),0)
# fosmean_hour=replace(fosmean_hour,is.nan(fosmean_hour),0)
# firmean_hour=replace(firmean_hour,is.nan(firmean_hour),0)


xvals <- -179.5:179.5
yvals <- -89.5:89.5

xdim <- ncdim_def( 'Lon', 'degree', xvals )
ydim <- ncdim_def( 'Lat', 'degree', yvals )

mv <- 0 # missing value
var_gpp <- ncvar_def( name="flux_gpp_prior_mean", units="mol/m2/sec", dim=list(xdim,ydim), longname="Gross Primary Productivity",missval=mv )
var_res <- ncvar_def( name="flux_res_prior_mean", units="mol/m2/sec", dim=list(xdim,ydim), longname="Total Ecosystem Respiration",missval=mv )
var_ocn <- ncvar_def( name="flux_ocean_prior_mean", units="mol/m2/sec", dim=list(xdim,ydim), longname="Ocean CO2 assimilation",missval=mv )
var_fos <- ncvar_def( name="flux_ff_prior_mean", units="mol/m2/sec", dim=list(xdim,ydim), longname="Fossil fuel CO2 emission",missval=mv )
var_fir <- ncvar_def( name="flux_fires_prior_mean", units="mol/m2/sec", dim=list(xdim,ydim), longname="Fires CO2 emission",missval=mv )


#if (procs==0)
#{
    output_fname=paste(outdir,"/","flux1x1_",tb2,"00","_",tb3,"00",".nc",sep="")

    ncid_new <- nc_create( output_fname, list(var_gpp,var_res,var_ocn,var_fos,var_fir))
    ncvar_put( ncid_new,var_gpp, gppmean_hour, start=c(1,1), count=c(ncol,nrow))
    ncvar_put( ncid_new,var_res, resmean_hour, start=c(1,1), count=c(ncol,nrow))
    ncvar_put( ncid_new,var_ocn, ocnmean_hour, start=c(1,1), count=c(ncol,nrow))
    ncvar_put( ncid_new,var_fos, fosmean_hour, start=c(1,1), count=c(ncol,nrow))
    ncvar_put( ncid_new,var_fir, firmean_hour, start=c(1,1), count=c(ncol,nrow))
    ncatt_put( ncid_new, 0, "Institution", "Centre for Isotope research, University of Groningen")
    ncatt_put( ncid_new, 0, "Contact", "wei.he@rug.nl")
    ncatt_put( ncid_new, 0, "Source", "CTDAS-STILT 1.0 background fluxes, generated from SiB3, SiBCASA and CT2013B products")      

    date = format(Sys.time(), "%b %d, %Y")
    user = Sys.getenv("LOGNAME")
    history = paste("created on",date,"by",user,sep=" ")
    ncatt_put( ncid_new, 0, "History", history)

    nc_close( ncid_new )

}
# Time marker
TEnd=proc.time()
Tot<-TEnd-TBegin
print(Tot)
