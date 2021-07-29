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

source("/Users/wei/co2simu/Rcode/output.ncdf4.v2.3.r")

args <- commandArgs(trailingOnly = TRUE)
procs= as.integer(args[1])


#inputs and outputs
stiltfile = paste("stilt_",procs,".rc", sep="")

curdir="/Storage/CO2/wei/test_additive_fluxfac/ctdas-stilt-base-covL1000new-moreairsites-airlocR/exec/da/stilt/"
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
boundary_flag= "CTNA"     #"CTNA", "CTE", "EMP"

# optimization option (flux only, or flux + boundary)
opt_para_flag = "FLUX"    # "FLUX_BOU","FLUX"

#----------------------------------------------------------------------------------------------
#Convolve footprints with sib hourly fluxes, for observations from both Aircraft and towers
#----------------------------------------------------------------------------------------------
endtime=240                   #hourly for 10 days
foottimes=seq(0,endtime,1)   #vector of times (backtimes) in hours between which footprint is computed
zbot=0                         #lower vertical bound for influence projection, in meters agl
ztop=0                         #upper vertical bound for influence projection, in meters agl
#if ztop set to zero, *surface* influence will be calculated

#set up an equivalent domain among footprints,fluxes, but not CO2 boundary 
ncol2=120  #66                       #number of pixels in x directions in grid 
nrow2=70   #40                       #number of pixels in y directions in grid 
LLLon2=-169.5  #(-129)                  #lower left corner of grid 
LLLat2=10.5    #22                      #lower left corner of grid 
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
flog.info("======Launch WRF-STILT-based forward simulations======", name='logger.b')

flog.info("The biosphere flux is from %s", bioflux_flag, name='logger.b')
flog.info("The boundary is from %s", boundary_flag, name='logger.b')


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
nsf=4    # 4 parameters for BC adjustment
scalefacarr_bc1=array(NA, dim=c(nsf,nsam))
scalefacarr_bc2=array(NA, dim=c(nsf,nsam))
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
  {
    scalefacarr1[,,i+1] = 0
    scalefacarr_bc1[,i+1] = 0
  }
  if (b1>=b0)  #b1>b0
  {
    ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb1,"00","_",tb2,"00",".nc",sep="")) 
    scalefac <- ncvar_get(ncf,"parametermap",start=c(11,101),count=c(ncol2,nrow2))   #real52:117,113:152,start=c(52,113),count=c(66,40)
    scalefacarr1[,,i+1] = scalefac * 1e6 # to umol/m2/s
    scalefac <- ncvar_get(ncf,"parametervalues_bc")
    
    scalefacarr_bc1[,i+1] = scalefac
    print(scalefacarr_bc1[,i+1])

    nc_close(ncf)
  }
  
  ncf <- nc_open(paste(samdir,"parameters.",ii,".",tb2,"00","_",tb3,"00",".nc",sep="")) 
  scalefac <- ncvar_get(ncf,"parametermap",start=c(11,101),count=c(ncol2,nrow2))  
  scalefacarr2[,,i+1] = scalefac * 1e6 # to umol/m2/s
  scalefac <- ncvar_get(ncf,"parametervalues_bc")

  scalefacarr_bc2[,i+1] = scalefac
  print(scalefacarr_bc2[,i+1])

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
obs_id <- ncvar_get(ncf,"obs_id")

# filter before loop
pfbfns=NULL
fn=NULL
eventids=NULL
obsids=NULL
for(mm in 1:length(fns))  
{  
  
  eid=as.numeric(substring(fns[mm],48,53))
  for(i in 1:length(eventidarr))  
  { 
    if(eid==eventidarr[i])
    {
      pfbfns = c(pfbfns,fns[mm])
      obsids = c(obsids,obs_id[i])
      eventids = c(eventids,eventidarr[i])
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

#read EMP boundary monthly into memory

if(boundary_flag == "EMP")
{
    bouf=paste(bounddir,y1,"/",m1,"/","CO2.v201209_v1.init.stilt.",y1,m1,".txt",sep="")
    #print(bouf)
    val<-as.matrix(read.table(bouf,header=TRUE))
}

# Read precalculated BC footprints
#bcff=paste(footpdir,y1,"/",m1,"/","FP_BC_3km_time_spent_2.1.4bcs.after.norm.txt",sep="")
#bcfval<-as.matrix(read.table(bcff, header=TRUE))
#identarr<-bcfval[,1]
#bcfarr1<-bcfval[,3]
#bcfarr2<-bcfval[,4]
#bcfarr3<-bcfval[,5]
#bcfarr4<-bcfval[,6]
#bcfarr1=replace(bcfarr1,is.na(bcfarr1),0)
#bcfarr2=replace(bcfarr2,is.na(bcfarr2),0)
#bcfarr3=replace(bcfarr3,is.na(bcfarr3),0)
#bcfarr4=replace(bcfarr4,is.na(bcfarr4),0)

for(mm in 1:length(pfbfns))  
{
  #--------------------------------------------
  # STEP 1: Read footprints
  #--------------------------------------------
# flog.info('obsids: %s',obsids[mm], name='logger.b') 
# flog.info('eventids: %s',eventids[mm], name='logger.b')  #June 23 
  yr=substring(pfbfns[mm],6,9)
  mo=substring(pfbfns[mm],11,12)

  recefilebase=paste(footpdir,"receptor_info/receptor_info.",yr,"-",sep="")
  recefile=paste(recefilebase,mo,".txt",sep="")
  flog.info('receptor filename: %s',recefile, name='logger.b')
  recedata=as.matrix(read.table(recefile,header=TRUE))
  eidarr=recedata[,16] 
  eventid=substring(pfbfns[mm],48,53)
  typeflag=NULL
  len=length(eidarr)
  for(k in 1:len)
  {
      if(eidarr[k]==eventid)
      {
          typeflag = recedata[k,3]
          flog.info('type: %s', typeflag, name='logger.b')
          flog.info('eventid: %s', eid, name='logger.b')
      }
   }

  fn=paste(footpdir,yr,"/",mo,"/",pfbfns[mm],sep="")
  footp=load.ncdf(fn)  
  
  footnc=nc_open(fn)
  foot=ncvar_get(footnc,"foot1",start=c(1,1,1),count=c(ncol2,nrow2,-1))

  #print(mean(foot))
  
  # get particles
  #part=footp$partfoot 
  #partna=footp$partfootnames
  #colnames(part)=c("time","index","lat","lon","agl","grdht","foot","temp","temp0","swrad","zi","dens","pres","dmass")   #ndmass
  
  # get info form path srings
  ident=substring(pfbfns[mm],6,46)
  eventid=substring(pfbfns[mm],48,53) 
  
  print(ident)
  flog.info("Convolving for %s",ident, name='logger.b')
  #info=id2info(ident) 
  #print(info)

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
    #print(xx)
    yy=xx+(-foottimes*3600)   #reversed order
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
   
    udaystring=rev(udaystring)
    print("XXXudaystring 1st from ident name")
    print(udaystring)
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
    #hgtarr=endpts[,5] + endpts[,6]  #agl+grdht
    magl=endpts[,5]                 # meters above ground level
    masl=endpts[,5] + endpts[,6]    # meters above sea level 
    hgtarr=masl

    # analyze the dates to recognize how many data to read
    #dd=as.character(endptsdate)
    #yrlist=substring(dd,1,4)
    #monlist=substring(dd,6,7)
    #daylist=substring(dd,9,10)
    #hrlist=substring(dd,12,13)
    
    # get unique months and days
    #daystring=paste(yrlist,monlist,daylist, sep="")
    #udaystring=unique(daystring)
    #yrmonstring=paste(yrlist,monlist, sep="")
    #uyrmonstring=unique(yrmonstring)


    #print("daylist=") #it's wrong to use daylist independently, but to use daystring (yyyymmdd)
    #print(daylist)

    #--------------------------------------------------------------
    # 2-1: Read fluxes and boundary data
    #--------------------------------------------------------------
    ndays=length(udaystring)

    if(boundary_flag=="CTNA")
        bouarr=array(0,dim=c(ndays,120,90,25,8))  #boundary : lon,lat,hgt,time
    
    if(boundary_flag=="CTE")
        bouarr=array(0,dim=c(ndays,360,180,25,8))
    
    #biospheric fluxes	
    
    if(bioflux_flag == "SiB3")
    {
        nrow=181 #181
        ncol=360 #361 #288

        gpp=array(NA, dim=c(ndays,ncol,nrow,24))
        rec=array(NA, dim=c(ndays,ncol,nrow,24))
    }

    #other fluxes #(360,180,8)	
    nrow=180
    ncol=360
    ocn=array(NA, dim=c(ndays,ncol,nrow,8))   
    fos=array(NA, dim=c(ndays,ncol,nrow,8))
    fir=array(NA, dim=c(ndays,ncol,nrow,8))

    if(bioflux_flag == "SiBCASA")
    {
        gpp=array(NA, dim=c(ndays,ncol,nrow,8))
        rec=array(NA, dim=c(ndays,ncol,nrow,8))
    }

    if(bioflux_flag == "CT_OPT")
        bio=array(NA, dim=c(ndays,ncol,nrow,8))
    
    ntimes=ndays*24  

    for(d in 1:ndays) 
    {
      datestr=udaystring[d]
      yr=substr(datestr,1,4)
      mn=substr(datestr,5,6)
      dy=substr(datestr,7,8)

      #print("#####################")
      #print(datestr)

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
          #flog.info('d: %s, yr: %s, mn: %s, dy: %s',d,yr,mn,dy,name='logger.b')
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
    dateall=rep(ISOdatetime(0,0,0,0,0,0,tz="UTC"), ntimes)
        
    if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
    {
        gppflux=array(NA, dim=c(ncol2,nrow2,ntimes))
        recflux=array(NA, dim=c(ncol2,nrow2,ntimes))
    }
    
    if(bioflux_flag == "CT_OPT")
        bioflux=array(NA, dim=c(ncol2,nrow2,ntimes))
   
    ocnflux = array(NA, dim=c(ncol2,nrow2,ntimes))
    fosflux = array(NA, dim=c(ncol2,nrow2,ntimes))
    firflux = array(NA, dim=c(ncol2,nrow2,ntimes))

    neeflux   = array(NA, dim=c(ncol2,nrow2,ntimes))
    neeflux1  = array(NA, dim=c(ncol2,nrow2,ntimes))
    neefluxarr= array(NA, dim=c(ncol2,nrow2,ntimes,nsam))
     
    bc1  = array(NA, dim=c(nsf,ntimes))
    bcarr= array(NA, dim=c(nsf,ntimes,nsam))

    neeoutarr= array(NA, dim=c(nsam))
    bcoutarrtmp = array(NA, dim=c(nsf,nsam))
    bcoutarr = array(NA, dim=c(nsam))
    deltaco2arr = array(NA, dim=c(nsam))
    fbouarr  = array(NA, dim=c(nsam))
    fsimuarr = array(NA, dim=c(nsam))
   
   #----------------------------------------------------
    yr=rep(sibyr,ntimes)
    mon=rep(sibmon,ntimes)
    hrs=seq(1,ntimes,1) 
    dy=ceiling(hrs/24)
    hr=(hrs-(dy-1)*24)*1-1
    time=ISOdatetime(yr,mon,dy,hr,0,0,tz="UTC")  
    
    for(hh in ntimes:1)
    {
        dateall[ntimes-hh+1]=time[hh]
        
        if(bioflux_flag == "SiB3")
        {      
           inxd=ceiling(hh/24)
           inxh=hh-(inxd-1)*24
           gppflux[,,ntimes-hh+1]=gpp[inxd,11:130,101:170,inxh]
           recflux[,,ntimes-hh+1]=rec[inxd,11:130,101:170,inxh]  
        }
            
        inxd=ceiling(hh/24)
        inxh=ceiling( (hh-(inxd-1)*24)/3 )   #every three hours the same
       
        #print(hh)
        #print(inxh)
        #print(dim(ocn))
        ocnflux[,,ntimes-hh+1]=ocn[inxd,11:130,101:170,inxh]
        fosflux[,,ntimes-hh+1]=fos[inxd,11:130,101:170,inxh]
        firflux[,,ntimes-hh+1]=fir[inxd,11:130,101:170,inxh]
        
        if(bioflux_flag == "SiBCASA")
        {
           gppflux[,,ntimes-hh+1]=gpp[inxd,11:130,101:170,inxh]
           recflux[,,ntimes-hh+1]=rec[inxd,11:130,101:170,inxh]
           #flog.info('gppflux: %s, gpp: %s, %s, %s, %s',gppflux[30,30,ntimes-hh+1],gpp[inxd,82,143,inxh],ntimes-hh+1,inxd,inxh,name='logger.b')
        }

        if(bioflux_flag == "CT_OPT")
           bioflux[,,ntimes-hh+1]=bio[inxd,11:130,101:170,inxh]
    }
    
    # replace NA values with 0, as NA was used as zero values for datasets
    gppflux[,,]=replace(gppflux[,,],is.na(gppflux[,,]),0)
    recflux[,,]=replace(recflux[,,],is.na(recflux[,,]),0)
        
    
    ocnflux[,,]=replace(ocnflux[,,],is.na(ocnflux[,,]),0)
    fosflux[,,]=replace(fosflux[,,],is.na(fosflux[,,]),0)
    firflux[,,]=replace(firflux[,,],is.na(firflux[,,]),0)
    
     if(bioflux_flag == "CT_OPT")
         bioflux[,,]=replace(bioflux[,,],is.na(bioflux[,,]),0)
      

     if(bioflux_flag == "SiB3")
         neeflux[,,]=recflux[,,]-gppflux[,,]  # nee*1e6 for SiBCASA?

     if(bioflux_flag == "SiBCASA")
         neeflux[,,]=(recflux[,,]-gppflux[,,])*1e6  # to umol/m2/s

     if(bioflux_flag == "CT_OPT")
         neeflux[,,]= bioflux[,,]*1e6  # for scaling CT flux

    #flog.info("nee flux, max=%f, min=%f, mean=%f", max(neeflux),min(neeflux),mean(neeflux), name='logger.b')

    #--------------------------------------------------------------
    # 2-3: Convolving footprints with fluxes to get delta co2
    #--------------------------------------------------------------
    
    #**********************************
    # biosphere fluxes convolving
    #**********************************
    # hourly flux
    #print("XXXdateall")
    #print(dateall)  #reversed order
    
#    # note: ixflux should not be from 1:ntimes, should be more accurate in hours
#    dd=as.character(min(yy))
#    #yrlist=substring(dd,1,4)
#    #monlist=substring(dd,6,7)
#    #daylist=substring(dd,9,10)
#    hr=substring(dd,12,13)
#    mi=substring(dd,15,16)

#    n = floor(as.numeric(hr) + as.numeric(mi)/60)
#    n = 24-n
#    print("XXXn") 
#    print(n)
#    datevalid=dateall[1]+0.5*3600-(0:(ntimes-1))*3600
#    #ixflux=(1:ntimes)[(datevalid<=xx)][1:endtime]   #time backward
#    #ixflux=(ntimes:1)[(datevalid<=xx)][(endtime+n-1):n]
#    #ixflux=(ntimes:1)[(datevalid<=xx)][(endtime+n-1):(n)]
#    ixflux<-c(n:endtime+n-1)
##    ixflux=(1:ntimes)[(datevalid<=xx)][n:(endtime+n-1)]


     #dd=as.character(min(yy))
     dd=as.character(xx)
     #yrlist=substring(dd,1,4)
     #monlist=substring(dd,6,7)
     #daylist=substring(dd,9,10)
     hr=substring(dd,12,13)
     mi=substring(dd,15,16)

     #flog.info('dd: %s, length dd: %s,hr: %s, mi: %s',dd,length(dd),hr,mi,name='logger.b')

     n = 24-floor(as.numeric(hr) + as.numeric(mi)/60)
     if (nchar(dd,type='width') ==10)
     {
         n=24
         #flog.info('No hr and mi, set n to 24',name='logger.b')
     }

     #print("XXXn")
     #print(n)
     ##datevalid=dateall[1]+0.5*3600-(0:(ntimes-1))*3600
     #ixflux=(1:ntimes)[(datevalid<=xx)][1:endtime]   #time backward
     #ixflux=(1:ntimes)[(datevalid<=xx)][1:(endtime+n)]
    # flog.info('n: %s, endtime: %s',n,endtime,name='logger.b')
     ixflux <-array(n:(endtime+n-1),dim=c(240))

    #flog.info('ixflux endtime: %s,%s,%s',endtime,n,as.numeric(hr)+(as.numeric(hr)/60),name='logger.b')

    #for (ttt in 1:240)  
    #  flog.info('ixflux: %s,%s',ttt,ixflux[ttt],name='logger.b')



 
    if(bioflux_flag == "SiB3")   #footprint is backward
    {
        gpptmp=foot*gppflux[,,ixflux]
        recotmp=foot*recflux[,,ixflux]
    
        gppout=sum(gpptmp,na.rm=TRUE)
        recoout=sum(recotmp,na.rm=TRUE)

        remove(list=c("gpptmp","recotmp"))
        gc()
    }
    #dateall[ixflux], datevalid from dateall1, from time1, from SiB flux; xx from from footprint
    
    #****************************************************
    # scaling NEE, determine which fluxes need scaling
    #****************************************************
    
    # involve scaling factors to neeflux[,,ixflux], partly, for example as follow
    # xx = "2010-01-24 18:10:00 UTC" || b1 = "2010-01-19 12:00:00 UTC", b3 = "2010-01-29 12:00:00 UTC"
    
    #xxleft=xx-ndays*24*3600
    xx=ISOdatetime(inityr,initmo,initdy,inithr,initmi,0,tz="UTC")
    yy=xx+(-foottimes*3600)
    xxleft=min(yy)

    #flog.info("xxleft = %s, b2=%s, b2-xxleft=%s", xxleft, b2, b2-xxleft,name='logger.b')
   

    #if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
    #    flog.info("xxleft = %s, b2=%s, b2-xxleft=%s", xxleft, b2, b2-xxleft,name='logger.b')

    for (i in 1:nsam) 
    {
      if(xxleft<b2)  #flux with scaling factors located in first lag
      {    
        diff=abs(as.numeric(difftime(b2,xxleft,units='hours')))+24

        #if(bioflux_flag == "SiBCASA" || bioflux_flag == "CT_OPT")
        #    diff = diff/3.0
        diff=ceiling(diff)
        align=24*ceiling((ntimes-diff)/24)#-(ntimes-diff) 
        #if (i==1)
        #   flog.info('diff: %s, align: %s',diff,align,name='logger.b')
        #if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
        #   flog.info("after diff = %f", diff, name='logger.b')

        for (hh in 1:(align)) 
            if(hh<=align && hh>=1)
            {
              #if(i==1)              
              #   flog.info("hh: %s, neeflux[,,hh]=%f, scalefacarr2 = %f ",hh, neeflux[30,30,hh], scalefacarr2[30,30,i], name='logger.b')
                
                
               neeflux1[,,hh] = neeflux[,,hh]+(scalefacarr2[,,i])   #delete "t" transform
               
               dims=dim(scalefacarr_bc2)
               #flog.info("dims of scalefacarr_bc2: %s, %s", dims[1], dims[2], name='logger.b')
               
               bc1[,hh] = scalefacarr_bc2[,i]
             }

        for (hh in ( (align+1):ntimes ) )
            if(hh<=ntimes && hh>=align+1)
            {
                #if(i==1)
                #   flog.info("hh: %s, neeflux[,,hh]=%f, scalefacarr1 = %f",hh, neeflux[30,30,hh], scalefacarr1[30,30,i],name='logger.b')
                neeflux1[,,hh] = neeflux[,,hh]+(scalefacarr1[,,i])
                
                dims=dim(scalefacarr_bc1)
                #flog.info("dims of scalefacarr_bc1: %s, %s", dims[1], dims[2], name='logger.b')
                
                bc1[,hh] = scalefacarr_bc1[,i]
            }

        neefluxarr[,,,i] = neeflux1[,,]  
        dims=dim(bc1)
        #flog.info("dims of bc1: %s, %s", dims[1], dims[2], name='logger.b')

        bcarr[,,i]=bc1[,]
      }      
      else  #all scaling within second lag
      {
        for (hh in 1:ntimes) 
          neeflux1[,,hh] = neeflux[,,hh]+(scalefacarr2[,,i])
          bc1[,hh] = scalefacarr_bc2[,i]
        
        neefluxarr[,,,i] = neeflux1[,,]
        
        dims=dim(bc1)
        #flog.info("dims of bc1: %s, %s", dims[1], dims[2], name='logger.b')

        bcarr[,,i] = bc1[,]
      }
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
        #neetmp2=foot*neefluxarr[,,6:245,i]

        #if(i==1)
            #flog.info("neetmp = %f, %f",sum(neetmp,na.rm=TRUE), sum(neetmp2,na.rm=TRUE),name='logger.b')

        neeout=sum(neetmp,na.rm=TRUE)
        neeoutarr[i]=neeout
       
        dims=dim(bcarr)
        #flog.info("dims of bcarr: %s, %s, %s", dims[1], dims[2], dims[3], name='logger.b')
        bcoutarrtmp[,i]=bcarr[,1,i]    #nsf,ntimes,nsam
      
        #if(ident=="2010x01x22x19x27x40.0500Nx105.0040Wx00300")
        #    flog.info("neeoutarr[i] = %f",neeoutarr[i], name='logger.b')
    }

    remove(list=c("neetmp","neeout"))
    gc()
    
    #**********************************
    # component fluxes convolving
    #**********************************
    # bg fluxes, 3-hourly
    #datevalid=dateall2[1]+1.5*3600-(0:(ntimes2-1))*3*3600	
    #ixflux=(1:ntimes2)[(datevalid<=xx)][1:endtime]
    
    #need large amount of memories
    ocntmp=foot*ocnflux[,,ixflux]*1e6
    fostmp=foot*fosflux[,,ixflux]*1e6
    firtmp=foot*firflux[,,ixflux]*1e6

    ocnout=sum(ocntmp,na.rm=TRUE)
    fosout=sum(fostmp,na.rm=TRUE) 
    firout=sum(firtmp,na.rm=TRUE)

    #print("XXX Line 640")
    #print(dim(foot))
    #print(dim(fosflux))  #tailed by line 520
    #print(ixflux)
    #print(fosout)

    if(bioflux_flag == "CT_OPT")
    {
        biotmp=foot*bioflux[,,ixflux]*1e6
        bioout=sum(biotmp,na.rm=TRUE)
    }
    
    foot2=foot 
    remove(list=c("ocnflux","fosflux","firflux","ocntmp","fostmp","firtmp","foot","dateall"))
    gc()
    
    #--------------------------------------------------------------
    # 2-4: Calculate boundary  
    #--------------------------------------------------------------
    #case 1: data from CT2013B
    #if(boundary_flag == "CTNA" || boundary_flag == "CTE")
    #{
      # boundary domain : lat 22.5~61.5, lon -128.5~-63.5
      latmin = -90   #22.5
      latmax =  90   #61.5
      lonmin = -180  #-128.5
      lonmax = 180   #-63.5
      npts=dim(endpts)[1]   
    
      # calculate concentrations
      pbou=0
      nval=0
      #dmin=min(as.integer(substr(aystring,7,8)))
      #dmax=max(as.integer(substr(udaystring,7,8)))
      #dmin=min(endptsdate)
      #dmax=max(endptsdate)
 
      xx=ISOdatetime(inityr,initmo,initdy,inithr,initmi,0,tz="UTC")
      yy=xx+(-foottimes*3600)

      #print("XXXudaystring for selected flux period")
      #print(udaystring)
      #print("XXXdmin particles")
      #print(dmin)
      #print("XXXdmax particles")
      #print(dmax)
      #print("XXXdobs time")
      #print(xx)

      dd1=as.character(endptsdate)
      yrlist=substring(dd1,1,4)
      monlist=substring(dd1,6,7)
      daylist=substring(dd1,9,10)
      hrlist=substring(dd1,12,13)
      milist=substring(dd1,15,16)
      
      counter_height=0
      counter_latlon=0
      counter_latlonheight=0
      t_weight=0
      counter_bc1=0
      counter_bc2=0
      counter_bc3=0
      counter_bc4=0

      for(p in 1:npts) 
      { 
        #dy=0
        #goal=as.integer(daylist[p])
        #dy=ceiling( (as.integer(difftime(xx,endptsdate[p],units='hours')) )/24 )+1 #if boundary reversed
        dy=ceiling( (as.integer(difftime(endptsdate[p],min(yy),units='hours')) )/24 )+1
        #for(uu in 1:ndays) 
        #{
        #  dstr=udaystring[uu]
        #  ddy=as.integer(substr(dstr,7,8))
         
         #print(ddy)
          #print(goal)

         # if(ddy<=goal)
         #    dy=dy+1
        #}

        #print("XXXdy=")
        #print(dy)

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
        k0=hgtarr[p]
        k=1
        for(l in 1:25) 
        {
          if(k0 > levelhgt[l] && k0 <= levelhgt[l+1])
          {
              k=l+1
              break
          }
          if(k0 > levelhgt[25])
             k=25
        }
     
        # notice that in obspack.py,aircraft observations which are lower than 3 km have been excluded from event ids. Here take into account of particles higher than 3 km
         # aircraft footprints
        if (typeflag == "aircraft" && magl[p] > 3000)
        {
            counter_height = counter_height + 1
            if(latarr[i]< 10.5 && lonarr[i]>= (-169.5) && lonarr[i] < (-50.5))
                 counter_bc1 = counter_bc1 + 1
            if(latarr[i]> 80.5 && lonarr[i]>= (-169.5) && lonarr[i] < (-50.5))
                 counter_bc2 = counter_bc2 + 1
            if(lonarr[i]< (-169.5) && latarr[i] >= 10.5 && latarr[i] < 80.5)
                 counter_bc3 = counter_bc3 + 1
            if(lonarr[i]> (-50.5) && latarr[i] >= 10.5 && latarr[i] < 80.5)
                 counter_bc4 = counter_bc4 + 1
        }
        # surface footprints                
        if (typeflag == "surface") 
        {
              counter_height = counter_height +1
              if(latarr[i]< 10.5 && lonarr[i]>= (-169.5) && lonarr[i] < (-50.5))
                   counter_bc1 = counter_bc1 + 1
              if(latarr[i]> 80.5 && lonarr[i]>= (-169.5) && lonarr[i] < (-50.5))
                   counter_bc2 = counter_bc2 + 1
              if(lonarr[i]< (-169.5) && latarr[i] >= 10.5 && latarr[i] < 80.5)
                   counter_bc3 = counter_bc3 + 1
              if(lonarr[i]> (-50.5) && latarr[i] >= 10.5 && latarr[i] < 80.5)
                   counter_bc4 = counter_bc4 + 1
        }


        #if (lonarr[p] < -125 || lonarr[p] > -55 || latarr[p] < 20 || latarr[p] > 75)
        # {
        #    counter_latlon = counter_latlon +1
        #    if (k0>3000)
        #    {
        #        counter_latlonheight = counter_latlonheight +1
                
        #        diff = as.numeric(difftime(xx,endptsdate[i],units='hours'))
        #        #flog.info("time difference = %s", diff, name='logger.b')
                
        #        t_weight = t_weight + (1-diff/240)/npts #240)/500

        #    }
        # }
        # 3-hourly matching /agt/alt
        
        if(boundary_flag == "CTNA" || boundary_flag == "CTE")
        {
           hpos=ceiling((as.integer(hrlist[p])+1e-6+as.integer(milist[p])/60.0 )/3.0)

           #flog.info('p3:, %s, pbou: %s, dy: %s, i: %s, j: %s, k: %s, hpos: %s',p,pbou,dy,i,j,k,hpos,name='logger.b')

           tt=bouarr[dy,i,j,k,hpos]  #notice bouarr not be reversed
      
           if(length(tt)==1)   #why sometimes we get a unnormal array?
           {
              pbou=pbou+tt  #sum  
              nval=nval+1
           }
         }#end if
      }#end for

      if(boundary_flag == "CTNA" || boundary_flag == "CTE")
      {
         #print("pbou")
         #print(pbou)
         #print(nval)
         fbou=pbou/nval
      }

    #FP_BC_c = counter_latlonheight/npts*t_weight

    #bcoutarr = FP_BC_c*bcoutarr

    flog.info('sum of counters : %s',(counter_bc1+counter_bc2+counter_bc3+counter_bc4),name='logger.b')
    flog.info('num of particles : %s',npts,name='logger.b')
   
    sum=counter_bc1+counter_bc2+counter_bc3+counter_bc4

    # change code here, to read 4 precalculated BC footprints from files
    if(sum>0)
    {
        FP_BC_c1 = counter_bc1/sum  #npts   #*t_weight
        FP_BC_c2 = counter_bc2/sum  #npts   #*t_weight
        FP_BC_c3 = counter_bc3/sum  #npts   #*t_weight
        FP_BC_c4 = counter_bc4/sum  #npts   #*t_weight
    }
    else
    {
        FP_BC_c1=0
        FP_BC_c2=0
        FP_BC_c3=0
        FP_BC_c4=0
    
    }

    dims=dim(bcoutarrtmp)
    #flog.info('dims of bcoutarrtmp : %s, %s',dims[1],dims[2],name='logger.b')
    
    bcoutarr = FP_BC_c1*bcoutarrtmp[1,]+FP_BC_c2*bcoutarrtmp[2,]+FP_BC_c3*bcoutarrtmp[3,]+FP_BC_c4*bcoutarrtmp[4,]
    len=length(bcoutarr)
    #print(bcoutarr)
    flog.info('FP_BC_c1,2,3,4:%f,%f,%f,%f',FP_BC_c1,FP_BC_c2,FP_BC_c3,FP_BC_c4,name='logger.b')

    flog.info('value of bcoutarr[1] : %s',bcoutarr[1], name='logger.b')
    flog.info('length of bcoutarr : %s',len,name='logger.b')

    #flog.info('FP_BC_c based on lat, lon and altitude: %s',FP_BC_c,name='logger.b')
    #flog.info('standard deviation BC*FP_BC: %s',sd(FP_BCARR_c),name='logger.b')


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
   

    #for (i in 1:nsam)
    #{
    #    if(xxleft<b2)
    #        fbouarr[i] = fbou + (scalefacarr_bc1[i])
    #    else
    #        fbouarr[i] = fbou + (scalefacarr_bc2[i])
    #}

    #################################################
    # final results and output to files
    #################################################
    
    fnameout1=paste(outdir,"/","samples_simulated.",tb2,"00","_",tb3,"00","-",procs,".txt",sep="")
    fn=fnameout1
    outlist=NULL #c(inityr, initmo, initdy, inithr, initmi, inihgt, fbou)
    for(i in 1:nsam)
    { 
      deltaco2arr[i]=neeoutarr[i]+ocnout+fosout+firout
      
      if(typeflag == "aircraft")
          fsimuarr[i]=(fbou+bcoutarr[i])*1e-6
      if(typeflag == "surface")
          fsimuarr[i]=(fbou+bcoutarr[i]+deltaco2arr[i])*1e-6        

      if(i==1)
         flog.info("fsimuarr[i] = %f, fbou = %f,  bcoutarr[i] = %f, neeoutarr[i] = %f, deltaco2arr[i] = %f",fsimuarr[i]*1e6,fbou,bcoutarr[i],neeoutarr[i],deltaco2arr[i],name='logger.b')

      outlist<-c(outlist, fsimuarr[i])
    }

    for(i in 1:nsam)
       outlist<-c(outlist, deltaco2arr[i]*1e-6)   #neeoutarr[i]

    for(i in 1:nsam)
       outlist<-c(outlist, bcoutarr[i]*1e-6)
    
    #print(deltaco2arr)

    # output bcoutarr to netcdf files
    #fname=paste(outdir,"/","BC_state_adjustments_",tb2,"00","_",tb3,"00",".nc",sep="")

    #xdim <- ncdim_def( 'Lon', 'degree', 1:ncol2 )
    #ensdim <- ncdim_def( 'ensemble', 'unitless', 1:100 )

    #mv <- 0 # missing value
    #var_bc <- ncvar_def( name="BC_adjustments", units="ppm", dim=list(ensdim), longname="Ensemble of BC adjustments",missval=mv )
    #var_flux <- ncvar_def( name="conc_vari_flux_adjustments", units="ppm", dim=list(ensdim),longname="Ensemble of concentration variations by flux adjustments", missval=mv )
    
    #if(file.exists(fname))
    #    ncid_new <- nc_open(fname, write=TRUE)
    #else
    #    ncid_new <- nc_create(fname, list(var_bc, var_flux))

    #ncvar_put( ncid_new, var_bc, bcoutarr, start=c(1), count=c(100))
    #ncvar_put( ncid_new,var_flux,deltaco2arr, start=c(1), count=c(100))
    #ncatt_put( ncid_new, 0, "Institution", "Centre for Isotope Research, University of Groningen")
    #ncatt_put( ncid_new, 0, "Contact", "wei.he@rug.nl")
    #ncatt_put( ncid_new, 0, "Source", "CTDAS-WRF-STILT v1.0")

    #date = format(Sys.time(), "%b %d, %Y")
    #user = Sys.getenv("LOGNAME")
    #history = paste("created on",date,"by",user,sep=" ")
    #ncatt_put( ncid_new, 0, "History", history)

    #nc_close( ncid_new )


    #add component varibles:gpp,reco,ocn,fos,fir,boundary
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
    nmem=(dim(data)[2]-8)/3     #9
    vals=data[,2:(3*nsam+8)]  #151+8 colums
}

if(bioflux_flag == "CT_OPT")
{
    nmem=(dim(data)[2]-7)/3     #9
    vals=data[,2:(3*nsam+7)]  #151+8 colums
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
    ncatt_put( ncid_new, 0, "Source", "CarbonTracker-WRF-STILT 1.0 fluxes, generated from SiB3, SiBCASA and CT2013B products")      

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
