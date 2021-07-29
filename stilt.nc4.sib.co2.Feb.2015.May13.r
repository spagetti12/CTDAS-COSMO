# CarbonTracker Data Assimilation Shell (CTDAS) Copyright (C) 2017 Wouter Peters. 
# Users are recommended to contact the developers (wouter.peters@wur.nl) to receive
# updates of the code. See also: http://www.carbontracker.eu. 
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation, 
# version 3. This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this 
# program. If not, see <http://www.gnu.org/licenses/>. 

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

source("/Users/wei/co2simu/Rcode/output.ncdf4.r")

# #inputs and outputs
# sibdir="/Users/weihe/Documents/Concensimulation/sib3_fluxes/day1×1/"
# bgfdir="/Users/weihe/Documents/Concensimulation/bgfluxes/3-hourly/"
# footpdir="/Users/weihe/Documents/Concensimulation/stilt_footprints/"  #read arlyn's calculated data
# bounddir="/Users/weihe/Documents/Concensimulation/co2_boundary/co2_total/"
# outdir="/Users/weihe/Documents/Concensimulation/resultsforctdas/"
# 
# #scaling factors for 10-days time scale
# samdir="/Users/weihe/Documents/example_ctdas/input/" 

stiltfile="stilt.rc"
curdir="/Storage/CO2/wei/ctdas-stilt-proto/exec/da/stilt/"
fn=paste(curdir, stiltfile,sep="")
conf=as.matrix(read.table(fn,header=FALSE))

sibdir  = conf[1,3]
bgfdir  = conf[2,3]
footpdir= conf[3,3]
bounddir= conf[4,3]
outdir  = conf[5,3]   #output
samdir  = conf[6,3]   #input
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
#--------------------------------------------------------------------------------------------------------
#library(futile.logger)
#flog.appender(appender.file("/Users/weihe/Documents/Concensimulation/resultstest/stilt.log"), name='logger.b')

#flog.info("This writes to a logging %s", "file", name='logger.b')
#flog.warn("This statement has higher severity", name='logger.b')
#flog.fatal("This one is really scary", name='logger.b')
#--------------------------------------------------------------------------------------------------------
# site level: read scaling factors for all sites once
# member level: read boundary for all ensemble members once

#nsam= 150  
#ndayscycle=10   # one cycle of footprints
#year=2010
#month=1
#day=19

nsam=as.numeric(conf[7,3])    #get from stilt.txt file  =>replace with a dictionary style, not id order; also check parameters and report status
ndayscycle=as.numeric(conf[8,3])
startdate=conf[9,3]
year=as.numeric(substring(startdate,1,4))
month=as.numeric(substring(startdate,6,7))
day=as.numeric(substring(startdate,9,10))
#if(ndayscycle!=10)
 # flog.warn("Sorry, only support the case nndayscycle equals 10", name='logger.b') 

b2 = ISOdatetime(year,month,day,0,0,0,tz="UTC")
b1 = b2 - ndayscycle*24*3600      # for footprint
b3 = b1 + 2*ndayscycle*24*3600    # for flux

tb1 = paste(substring(b1,1,4),substring(b1,6,7),substring(b1,9,10),sep="")
tb2 = paste(substring(b2,1,4),substring(b2,6,7),substring(b2,9,10),sep="")
tb3 = paste(substring(b3,1,4),substring(b3,6,7),substring(b3,9,10),sep="")

scalefacarr1=array(NA, dim=c(ncol2,nrow2,nsam))   
scalefacarr2=array(NA, dim=c(ncol2,nrow2,nsam))   
# make CTDAS to generate domain for North America, read data partly?

#flog.info("Reading scaling factor files", name='logger.b')
for(i in 0:(nsam-1))    #parameters.000.2010010100_2010011100.nc
{
  if (i<10)
    ii = paste("00",i,sep="")
  if (i<100 && i>=10)
    ii = paste("0",i,sep="")
  if (i>=100)
    ii = i

  ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb1,"00","_",tb2,"00",".nc",sep="")) 
  scalefac <- ncvar_get(ncf,"parametermap",start=c(52,113),count=c(ncol2,nrow2))   #temporary 42  #real52:117,113:152,start=c(52,113),count=c(66,40)
  scalefacarr1[,,i+1]=scalefac
  
  ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb2,"00","_",tb3,"00",".nc",sep="")) 
  scalefac <- ncvar_get(ncf,"parametermap",start=c(52,113),count=c(ncol2,nrow2))  
  scalefacarr2[,,i+1]=scalefac  
}

#---------------------------------------------------------------------------------
# Centering at the state vector especially the ready-optimized cycle, look for data
# from corresponding period for optimzation
#---------------------------------------------------------------------------------
# according to b1, b2, decide which months
m1=substring(b2,6,7)
m2=substring(b3,6,7)
umon=unique(c(m1,m2))

#flog.info("Determing ready-use footprint files", name='logger.b')
fns=NULL
for(mm in 1:length(umon))  
{
  pfbpath=paste(footpdir,year,"/",umon[mm],"/",sep="")  #"stilt2010x01x01x18x59x33.4057Nx081.8334Wx00285.nc"
  tmp=list.files(pfbpath,pattern="nc", all=T)
  fns=c(fns,tmp)
}

# read needed event ids from sample_coordinates.***.nc files(e.g, sample_coordinates_2010011100_2010012100.nc)
ncf <- nc_open(paste(samdir,"sample_coordinates_",tb2,"00","_",tb3,"00",".nc",sep="")) 
eventidarr <- ncvar_get(ncf,"obs_num")

# filter before loop
pfbfns=NULL
fn=NULL
for(mm in 1:length(fns))  
{  
  #select: Uniform time to make it start from the same point
  #yr=substring(fns[mm],6,9)    #"2010-01-09 09:00:00 UTC"
  #mo=substring(fns[mm],11,12)
  #dy=substring(fns[mm],14,15)
  #bISO = ISOdate(yr,mo,dy)
  #diff1=as.numeric(bISO - b2)
  #diff2=as.numeric(bISO - b3)
  
  #if (diff1 >= 0 && diff2 <= 0)  # in the time span
  #  pfbfns = c(pfbfns,fns[mm])
  
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

#log=paste("For state vectors from ",b1,"to",b2,",",length(pfbfns),"observations have been found")
#flog.info(log, name='logger.b')

# for fluxes, the current 10 days using scaling factors
#flog.info("Start convolving for all footprint files", name='logger.b')
newfile=T #output results into a new file

for(mm in 1:length(pfbfns))  
{	  
  # Read pfb file
  #mm=1
  mo=substring(pfbfns[mm],11,12)
  fn=paste(footpdir,year,"/",mo,"/",pfbfns[mm],sep="")
  footp=load.ncdf(fn)  
  
  footnc <- nc_open(fn)
  foot <- ncvar_get(footnc,"foot1",start=c(41,12,1),count=c(ncol2,nrow2,-1))
  
  # Read footprint directly
  #ft=footp$foot1   #120*70*240  
  #foot=array(NA, dim=c(40,66,240))
  #for (i in 1:240)
  #  foot[,,i]=t(as.matrix(ft[41:106,12:51,i]))  #113:152
  
  # get info form path srings
  ident=substring(pfbfns[mm],6,46)
  eventid=substring(pfbfns[mm],48,53) 
  
  print(ident)
 # flog.info("Convolving for %s",ident, name='logger.b')
  info=id2info(ident)  #interpret info from file names 
  
  #--------------------------------------------------------
  # STEP 1: Creates footprint for individual particle runs
  #--------------------------------------------------------
  #foot1=Trajecfoot(ident=ident,part=part, pathname="",foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,
  #                 numpix.x=ncol1,numpix.y=nrow1,lon.ll=LLLon1,lat.ll=LLLat1,lon.res=ResLon1,lat.res=ResLat1)
  
  
  # for different spatial resolutions between biosphere fluxes and background fluxes, two footprints are planned to be used, and to the end, we picked up the concentration values.
  #foot2=Trajecfoot(ident=ident,part=part, pathname="",foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,
  #                 numpix.x=ncol2,numpix.y=nrow2,lon.ll=LLLon2,lat.ll=LLLat2,lon.res=ResLon2,lat.res=ResLat2)
  
  
  if(length(foot)>100)
  {			
    inityr=as.numeric(substring(ident,1,4))
    initmo=as.numeric(substring(ident,6,7))
    initdy=as.numeric(substring(ident,9,10))
    inithr=as.numeric(substring(ident,12,13))
    initmi=as.numeric(substring(ident,15,16))   #by He
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
    
    #kk=length(uyrmonstring)   #needs July (former month)
    
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
    bouarr=array(0,dim=c(ndays,120,90,25,8))  #boundary : lon,lat,hgt,time
    
    #biospheric fluxes	
    nrow_flux1=181
    ncol_flux1=361 #288
    gpp=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
    rec=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
    pla=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
    soi=array(NA, dim=c(ndays,ncol_flux1,nrow_flux1,24))
    
    #other fluxes #(360,180,8)	
    nrow_flux2=180
    ncol_flux2=360
    ocn=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))   
    fos=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    fir=array(NA, dim=c(ndays,ncol_flux2,nrow_flux2,8))
    
    ntimes1=ndays*24  
    ntimes2=ndays*8  
    
    for(d in 1:ndays) 
    {
      datestr=udaystring[d]
      yr=substr(datestr,1,4)
      mn=substr(datestr,5,6)
      dy=substr(datestr,7,8)
      #bou=load.ncdf(paste(bounddir,"CT2013B.molefrac_nam1x1_",yr,"-",mn,"-",dy,".nc",sep=""))
      bou=load.ncdf(paste(bounddir,"CT2013B.molefrac_glb3x2_",yr,"-",mn,"-",dy,".nc",sep=""))
      # co2(date, level, lat, lon)  /  ocn_flux_opt(date, lat, lon)
      
      biof=load.ncdf(paste(sibdir,"SiB3.hourly.flux1x1.global.",yr,mn,dy,".nc",sep=""))  
      bgf=load.ncdf(paste(bgfdir,"CT2013B.flux1x1.",yr,mn,dy,".nc",sep=""))
      
      bouarr[d,,,,]=bou$co2             # boundary CO2
      
      gpp[d,,,]=biof$gpp
      rec[d,,,]=biof$rtotal
      pla[d,,,]=biof$ocs.gpp
      soi[d,,,]=biof$ocs.soil
      
      ocn[d,,,]=bgf$ocn.flux.opt        # here we can read less, so change the sizes of ocnarr
      fos[d,,,]=bgf$fossil.flux.imp
      fir[d,,,]=bgf$fire.flux.imp
      
      remove(list=c("bou","biof","bgf"))
    }
    
    #--------------------------------------------------------------
    # 2-2: prepare data for calculations
    #--------------------------------------------------------------
    dateall1=rep(ISOdatetime(0,0,0,0,0,0,tz="UTC"), ntimes1) 
    dateall2=rep(ISOdatetime(0,0,0,0,0,0,tz="UTC"), ntimes2)

    gppflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    recflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    plaflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    soiflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    neeflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    neeflux1=array(NA, dim=c(ncol2,nrow2,ntimes1))
    #----------------------------------------------------
    neefluxarr=array(NA, dim=c(ncol2,nrow2,ntimes1,nsam))
    neeoutarr=array(NA, dim=c(nsam))
    fsimuarr=array(NA, dim=c(nsam))
    #----------------------------------------------------
    
    ocnflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    fosflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    firflux=array(NA, dim=c(ncol2,nrow2,ntimes1))
    
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
    
    for(hh in ntimes1:1)
    {
      dateall1[ntimes1-hh+1]=time1[hh]
      inxd=ceiling(hh/24)
      inxh=hh-(inxd-1)*24
      gppflux[,,ntimes1-hh+1]=gpp[inxd,52:117,113:152,inxh]  #transform matrix, 41:93; delete this on May 4,2015
      recflux[,,ntimes1-hh+1]=rec[inxd,52:117,113:152,inxh]  
      plaflux[,,ntimes1-hh+1]=pla[inxd,52:117,113:152,inxh]
      soiflux[,,ntimes1-hh+1]=soi[inxd,52:117,113:152,inxh]
    }
    neeflux[,,]=recflux[,,]-gppflux[,,]
    
    for(hh in ntimes2:1)
    {
      dateall2[ntimes2-hh+1]=time2[hh]
      inxd=ceiling(hh/8)
      inxh=hh-(inxd-1)*8
      ocnflux[,,ntimes2-hh+1]=ocn[inxd,52:117,113:152,inxh]  #transform matrix, delete this 
      fosflux[,,ntimes2-hh+1]=fos[inxd,52:117,113:152,inxh]  
      firflux[,,ntimes2-hh+1]=fir[inxd,52:117,113:152,inxh]   
    }	      
    
    #--------------------------------------------------------------
    # 2-3: Convolving footprints with fluxes to get delta co2
    #--------------------------------------------------------------
    # sib3 flux, 1-hourly
    datevalid=dateall1[1]+0.5*3600-(0:(ntimes1-1))*3600   #need to change, ntimes=1488	
    ixflux=(1:ntimes1)[(datevalid<=xx)][1:endtime]        #time backward
    #dateall[ixflux]
    
    #datevalid from dateall1, from time1, from SiB flux; xx from from footprint
    
    #need large amount of memories
    ocstmp=foot*plaflux[,,ixflux]
    soiltmp=foot*soiflux[,,ixflux]
    gpptmp=foot*gppflux[,,ixflux]     
    recotmp=foot*recflux[,,ixflux]   #dim(neetmp): 40  53 240
    
    #################################################
    # determine which fluxes need scaling factors
    #################################################
    
    # involve scaling factors to neeflux[,,ixflux], partly, for example as follow
    # xx = "2010-01-24 18:10:00 UTC" || b1 = "2010-01-19 12:00:00 UTC", b3 = "2010-01-29 12:00:00 UTC"

    xxleft=xx-ndays*24*3600
 
    for (i in 1:nsam) 
    {
      if(xxleft<b2)  #flux with scaling factors
      {    
        diff=as.numeric(difftime(b2,xxleft,units='hours'))  # needs scaling for first half
        diff=ceiling(diff)
        for (hh in 1:(diff)) 
          neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr1[,,i])   #delte "t" transform
 
        for (hh in (diff+1):ntimes1) 
          neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr2[,,i])    
        
        neefluxarr[,,,i] = neeflux1[,,]  
      }      
      else 
      {
        for (hh in 1:ntimes1) 
          neeflux1[,,hh] = neeflux[,,hh]*(scalefacarr2[,,i])    
        
        neefluxarr[,,,i] = neeflux1[,,]   
      }
    }

    # REMEMBER: the spatial resolutios of the scaling factors and fluxes are different, whiched needs to be uniformed afterwards 
    
    # delta co2 on NEE for all ensemble members
    for(i in 1:nsam)
    {
      neetmp=foot*neefluxarr[,,ixflux,i]   
      neeout=sum(neetmp,na.rm=TRUE)
      neeoutarr[i]=neeout
    }
    #------------------------------------------------
    
    ocsout=sum(ocstmp,na.rm=TRUE)
    soilout=sum(soiltmp,na.rm=TRUE) 
    gppout=sum(gpptmp,na.rm=TRUE)
    recoout=sum(recotmp,na.rm=TRUE)
    
    #remove(list=c("gppflux","recflux","soiflux","plaflux","gpptmp","recotmp","soiltmp","ocstmp","dateall1"))
    gc()
    
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
    
    #remove(list=c("ocnflux","fosflux","firflux","ocntmp","fostmp","firtmp","foot","dateall2"))
    gc()
    
    #--------------------------------------------------------------
    # 2-4: Calculate boundary  
    #--------------------------------------------------------------
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
      #if(lonarr[p]>=lonmin && lonarr[p]<=lonmax && latarr[p]>=latmin && latarr[p]<=latmax)
      #{
      #dy=as.integer(daylist[p])-dmin+1   #bug, who know if the dates are continuous or not
      dy=0
      goal=as.integer(daylist[p])
      for(uu in 1:ndays) 
      {
        dstr=udaystring[uu]
        ddy=substr(dstr,7,8)
        
        if(ddy<=goal)
          dy=dy+1
      }
      
      #i=ceiling((lonarr[p]-lonmin)/1.0+0.5)  #check if the position is right
      i=ceiling((lonarr[p]-lonmin)/3.0)
      #nval=nval+1  #notice!!! this (i,j) is for NA
      #notice: we maybe should use the data besides North America domain, not only nam1x1, but glb3x2...
      
      #j=ceiling((latarr[p]-latmin)/1.0+0.5)
      j=ceiling((latarr[p]-latmin)/2.0)
      
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
      
      #}
    }
    #print(nval)
    fbou=pbou/nval
    
    #deltaco2=(recoout-gppout)+ocnout+fosout+firout #contributions from the biosphere and background fluxes
    #fsimu=fbou+deltaco2 
    
    #################################################
    # final results and output to files
    #################################################
    
    fnameout1=paste(outdir,"/",tb2,"/","samples_simulated.",tb2,"00","_",tb3,"00",".txt",sep="")
    fn=fnameout1
    outlist=NULL #c(inityr, initmo, initdy, inithr, initmi, inihgt, fbou)
    for(i in 1:nsam)
    { 
      deltaco2=neeoutarr[i]+ocnout+fosout+firout 
      fsimuarr[i]=fbou+deltaco2        
      
      outlist<-c(outlist, fsimuarr[i])
    }
    
    out1<-c(ident,eventid,round(outlist,4))  
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
#source("/Users/weihe/Documents/Concensimulation/output.ncdf4.r")
#fn= "/Users/weihe/Documents/Concensimulation/resultstest/stilt_sib_ct_co2_2010-01-19_2010-01-29_simulated.txt"
fin=paste(fn,sep="")  #already include eventid 
fout=paste(substring(fn,1,nchar(fn)-4),".nc",sep="")
data=as.matrix(read.table(fin,header=FALSE))

nobs=dim(data)[1]
nmem=dim(data)[2]-2
vals=data[,2:(nsam+2)]  #151 colums
write.results.netcdf(vals,nobs,nmem,fout)   #write simulated results once
#-----------------------------------------------------------------------

#create stilt.ok file
file.create("stilt.ok")

# Time marker
TEnd=proc.time()
Tot<-TEnd-TBegin
print(Tot)
