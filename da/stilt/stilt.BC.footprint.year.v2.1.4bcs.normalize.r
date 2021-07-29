#------------------------------------------------------------------------------------------------
# Calculation for BC sensitivity for CarbonTracker-Lagrange
# including general case (lat/lon/alt/time) and normalized case, which is a try for solving  
# the problem of relative small values
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


# constants for level-height transfrom
levelhgt=c(34.5,111.9,256.9,490.4,826.4,1274.1,1839.0,2524.0,3329.9,4255.6,5298.5,6453.8,7715.4,9076.6,
           + 10533.3,12108.3,13874.2,15860.1,18093.2,20590.0,24247.3,29859.6,35695.0,42551.5,80000.0)
#data from http://www.esrl.noaa.gov/gmd/ccgg/carbontracker-ch4/documentation_tm5.html

library(ncdf4)
library(futile.logger)

logfile = "stilt_BC_footprint.log"
flog.appender(appender.file(logfile), name='logger.b')
flog.info("======Launch BC footprints calculations======", name='logger.b')

newfile=T 

datadir="/Storage/ctdas-wrfstilt/new_stilt_footprints/flask/2010/"
eidfilebase="/Storage/ctdas-wrfstilt/stilt_footprints/receptor_info/receptor_info.2010-"


for(n in 1:12)
{
    if(n<10)
        mon=paste("0",n,sep="")
    else
        mon=n
    
    pfbpath=paste(datadir,mon,"/",sep="")
    pfbfns=list.files(pfbpath,pattern="nc", all=T)


    for(mm in 1:length(pfbfns))  
    {
        #--------------------------------------------
        # STEP 1: Read footprints
        #--------------------------------------------
        fn=paste(pfbpath,pfbfns[mm],sep="")
        footp=load.ncdf(fn)  
  
        # get info form path srings
        ident=substring(pfbfns[mm],6,46)
        eventid=substring(pfbfns[mm],48,53) 
        height=as.numeric(substring(pfbfns[mm],42,46))

        print(ident)
        flog.info(" %s",ident, name='logger.b')
        #info=id2info(ident) 
        #print(info)
  
        #if(length(foot)>100)
        #{			
        inityr=as.numeric(substring(ident,1,4))
        initmo=as.numeric(substring(ident,6,7))
        initdy=as.numeric(substring(ident,9,10))
        inithr=as.numeric(substring(ident,12,13))
        initmi=as.numeric(substring(ident,15,16))
        inihgt=as.numeric(substring(ident,39,41))      #bug, height 38,41
        inihgt=as.numeric(substring(ident,37,41)) 
     
        # get the time stamp for each foot step, going back time given by "endtime"
        xx=ISOdatetime(inityr,initmo,initdy,inithr,initmi,0,tz="UTC")
    

        # get endpoints
        endpts=footp$endpts
        endptna=footp$endptsnames 
        #endpts=ncvar_get(footnc,"endpts")
        #endptna=ncvar_get(footnc,"endptsnames")	
        colnames(endpts)=endptna   #=c("time","index","lat","lon","agl","grdht","temp","pres")  
    
        endptsdate=footp$endptsdate  #2010-08-04 20:18:00 UTC
        #endptsdate=ncvar_get(footnc,"endptsdate")
    
        #diff = as.numeric(difftime(xx,endptsdate,units='hours'))
        #flog.info("time difference = %s", diff, name='logger.b')

        latarr=endpts[,3]
        lonarr=endpts[,4]
        hgtarr=endpts[,5] + endpts[,6]  #agl+grdht

        #The domain we want to consider is a box around North America, with boundaries at
        #125W, 55W, 20N, 75N, and 3km altitude
        num_in1=0
        num_in2=0
        num_in3=0
        num_out1=0
        num_out2=0
        num_out3=0
        t_weight1=0
        
        counter_bc1=0
        counter_bc2=0
        counter_bc3=0
        counter_bc4=0

        for (i in 1:length(lonarr))
        {
           
           if (hgtarr[i]>3000)
           {
                counter_height = counter_height +1
                if(latarr[i]<20 && lonarr[i]>= (-125) && lonarr[i] <= (-55))
                     counter_bc1 = counter_bc1 + 1
                if(latarr[i]>75 && lonarr[i]>= (-125) && lonarr[i] <= (-55))
                     counter_bc2 = counter_bc2 + 1
                if(lonarr[i]< (-125) && latarr[i] >= 20 && latarr[i] <= 75)
                     counter_bc3 = counter_bc3 + 1
                if(lonarr[i]> (-55) && latarr[i] >= 20 && latarr[i] <= 75)
                     counter_bc4 = counter_bc4 + 1
           }

           if(latarr[i]<20 || latarr[i]>75 || lonarr[i]< -125 || lonarr[i]> -55) # && hgtarr[i]<=3000 )
           {
              if(hgtarr[i]>3000)
              {
                #print(c(latarr[i],lonarr[i],hgtarr[i]))
                num_out3=num_out3+1
          
               diff = as.numeric(difftime(xx,endptsdate[i],units='hours'))
               flog.info("time difference = %s", diff, name='logger.b')
          
               # linear
               t_weight1 = t_weight1 +  (1-diff/240)/500
              }      

           }

           if(height<500)
               flag="T"
           else
               flag="A"
       }


    fnameout1=paste(pfbpath,"FP_BC_3km_time_spent_2.1.4bcs.prep01.norm.txt",sep="")
    flog.info("Output to: ",fnameout1, name='logger.b')
    outlist<-c(counter_bc1/500*t_weight1,counter_bc2/500*t_weight1,counter_bc3/500*t_weight1,counter_bc4/500*t_weight1)
    out1<-c(ident,eventid,round(outlist,6),flag)  
    if(newfile) 
      write.table(t(out1), file=fnameout1, append=F, col.names=F, row.names=F, quote=F)
    if(!newfile) 
      write.table(t(out1), file=fnameout1, append=T, col.names=F, row.names=F, quote=F)
    newfile=F
    
     #}# end if foot NA
  
  }#end loop mm
  
  # for every month, normalize 

  fnameout2=paste(pfbpath,"FP_BC_3km_time_spent_2.1.4bcs.prep02.norm.txt",sep="")
  flog.info("Output to: ",fnameout2, name='logger.b')

  footfile<-fnameout1
  print(footfile)
  footdata<-as.matrix(read.table(footfile,header=FALSE))
  ident=footdata[,1]
  evn1=footdata[,2]
  bcf_1=footdata[,3]
  bcf_2=footdata[,4]
  bcf_3=footdata[,5]
  bcf_4=footdata[,6]
  d1=dim(footdata)
  lin1=d1[1]
  #sitenew=NULL
  #airtowflag=NULL

  # according to evn, find site name, and using these recorded infomation to normalize
  
  
  eidfile<-paste(eidfilebase,mon,".txt",sep="")
  eiddata<-as.matrix(read.table(eidfile,header=TRUE))
  evn2=eiddata[,15]
  site=eiddata[,11]
  airtower=eiddata[,3]

  d2=dim(eiddata)
  lin2=d2[1]

  #out2=NULL

  newfile2=T
  for (i in 1:lin1)
  {
     id1 = evn1[i]
    
     sitenew=NULL
     airtowflag=NULL
     out2=NULL

     for (j in 1:lin2)
     {
       id2 = evn2[j]
       
       if(id1==id2)
       {
           sitenew=site[j]
           airtowflag=airtower[j]
           break
       }
     }

    out2=c(ident[i],evn1[i],bcf_1[i],bcf_2[i],bcf_3[i],bcf_4[i],sitenew,airtowflag)
    
    header=c("ident","eventid","bc1","bc2","bc3","bc4","site","type")

    if(newfile2)
          #write.table(t(header), file=fnameout2, append=F, col.names=F, row.names=F, quote=F)
          write.table(t(out2), file=fnameout2, append=F, col.names=header, row.names=F, quote=F)
    if(!newfile2)
          write.table(t(out2), file=fnameout2, append=T, col.names=F, row.names=F, quote=F)
    newfile2=F

  }
 
  footfile<-fnameout2
  
  footdata<-as.matrix(read.table(footfile,header=TRUE))
  ident=footdata[,1]
  evn1=footdata[,2]
  #bcf=footdata[,3]
  site=footdata[,7]
  type=footdata[,8]

  d3=dim(footdata)
  lin3=d3[1]

  fnameout3=paste(pfbpath,"FP_BC_3km_time_spent_2.1.4bcs.after.norm.txt",sep="")
  flog.info("Output to: ",fnameout3, name='logger.b')
 

for(kk in 1:4)
{

  amtsum=0
  baosum=0
  lefssum=0
  sctsum=0
  strsum=0
  wbisum=0
  wgcsum=0
  wktsum=0
  carsum=0
  etlsum=0
  lefasum=0
  nhasum=0
  tgcsum=0
  thdsum=0
  amtcount=0
  baocount=0
  lefscount=0
  sctcount=0
  strcount=0
  wbicount=0
  wgccount=0
  wktcount=0
  carcount=0
  etlcount=0
  lefacount=0
  nhacount=0
  tgccount=0
  thdcount=0

  for(kk in 1:4)
  {
     bcf=footdata[,2+kk]

   for (i in 1:lin3)
   {

     if(type[i]=="surface" && site[i]=="AMT")
     {
         amtsum=amtsum+as.numeric(bcf[i])
         amtcount=amtcount+1
     }

     if(type[i]=="surface" && site[i]=="BAO")
     {
         baosum=baosum+as.numeric(bcf[i])
         baocount=baocount+1
     }

     if(type[i]=="surface" && site[i]=="LEF")
     {
         lefssum=lefssum+as.numeric(bcf[i])
         lefscount=lefscount+1
     }

     if(type[i]=="surface" && site[i]=="SCT")
     {
         sctsum=sctsum+as.numeric(bcf[i])
         sctcount=sctcount+1
     }

     if(type[i]=="surface" && site[i]=="STR")
     {
         strsum=strsum+as.numeric(bcf[i])
         strcount=strcount+1
     }

     if(type[i]=="surface" && site[i]=="WBI")
     {
         wbisum=wbisum+as.numeric(bcf[i])
         wbicount=wbicount+1
     }

     if(type[i]=="surface" && site[i]=="WGC")
     {
         wgcsum=wgcsum+as.numeric(bcf[i])
         wgccount=wgccount+1
     }

     if(type[i]=="surface" && site[i]=="WKT")
     {
         wktsum=wktsum+as.numeric(bcf[i])
         wktcount=wktcount+1
     }

     if(type[i]=="aircraft" && site[i]=="CAR")
     {
         carsum=carsum+as.numeric(bcf[i])
         carcount=carcount+1
     }

     if(type[i]=="aircraft" && site[i]=="ETL")
     {
         etlsum=etlsum+as.numeric(bcf[i])
         etlcount=etlcount+1
     }
     
     if(type[i]=="aircraft" && site[i]=="LEF")
     {
         lefasum=lefasum+as.numeric(bcf[i])
         lefacount=lefacount+1
     }
     
     if(type[i]=="aircraft" && site[i]=="NHA")
     {
         nhasum=nhasum+as.numeric(bcf[i])
         nhacount=nhacount+1
     }
     
     if(type[i]=="aircraft" && site[i]=="TGC")
     {
         tgcsum=tgcsum+as.numeric(bcf[i])
         tgccount=tgccount+1
     }

     if(type[i]=="aircraft" && site[i]=="THD")
     {
         thdsum=thdsum+as.numeric(bcf[i])
         thdcount=thdcount+1
     }
     
  }
  
  newfile3=T
  for (i in 1:lin3)  # for the site AMT, lin3 should be replaced with amtcount
  {
      tmp=bcf[i]  #footdata[,2+kk]
      if(type[i]=="surface" && site[i]=="AMT")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(amtsum/amtcount)
      if(type[i]=="surface" && site[i]=="BAO")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(baosum/baocount)
      if(type[i]=="surface" && site[i]=="LEF")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(lefssum/lefscount)
      if(type[i]=="surface" && site[i]=="SCT")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(sctsum/sctcount)
      if(type[i]=="surface" && site[i]=="STR")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(strsum/strcount)
      if(type[i]=="surface" && site[i]=="WBI")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(wbisum/wbicount)
      if(type[i]=="surface" && site[i]=="WGC")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(wgcsum/wgccount)
      if(type[i]=="surface" && site[i]=="WKT")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(wktsum/wktcount)
      if(type[i]=="aircraft" && site[i]=="CAR")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(carsum/carcount)
      if(type[i]=="aircraft" && site[i]=="ETL")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(etlsum/etlcount)
      if(type[i]=="aircraft" && site[i]=="LEF")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(lefasum/lefacount)
      if(type[i]=="aircraft" && site[i]=="NHA")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(nhasum/nhacount)
      if(type[i]=="aircraft" && site[i]=="TGC")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(tgcsum/tgccount)
      if(type[i]=="aircraft" && site[i]=="THD")
            footdata[i,2+kk]=as.numeric(footdata[i,2+kk])/(thdsum/thdcount)   
  }

}

for (i in 1:lin3)
{
      if(site[i]=="AMT"||site[i]=="BAO"||site[i]=="LEF"||site[i]=="SCT"||site[i]=="STR"||site[i]=="WBI"||site[i]=="WGC"
          ||site[i]=="WKT"||site[i]=="CAR"||site[i]=="ETL"||site[i]=="NHA"||site[i]=="TGC"||site[i]=="THD")
      {
          outlist=c(as.numeric(footdata[i,3]),as.numeric(footdata[i,4]),as.numeric(footdata[i,5]),as.numeric(footdata[i,6]))
          out3=c(ident[i],evn1[i],tmp, round(outlist,6), site[i], type[i])
     
          header=c("ident","eventid","bc1","bc2","bc3","bc4","site","type")
          if(newfile3)
                write.table(t(out3), file=fnameout3, append=F, col.names=header, row.names=F, quote=F)
          if(!newfile3)
                write.table(t(out3), file=fnameout3, append=T, col.names=F, row.names=F, quote=F)
      newfile3=F
      }

  }

  


}

