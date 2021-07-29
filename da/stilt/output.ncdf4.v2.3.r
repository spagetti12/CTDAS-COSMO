#******************************************************************************************
# Function for writing a netcdf4 file
# created by W.He on April 15, 2015 
#****************************************************************************************** 

write.results.netcdf <- function(bioflux_flag,vals,nobs,nmem,output_fname) { 
    
  obsdim <- ncdim_def( 'obsnum', '', 1:nobs, unlim=FALSE, create_dimvar=FALSE ) #IMPORTANT TO USE create_dimvar
  memdim <- ncdim_def( 'nmembers', '', 1:nmem, unlim=FALSE, create_dimvar=FALSE)
  
  ccgg_evn=array(NA,dim=c(nobs))
  flask=array(NA,dim=c(nmem,nobs))
  deltacarbon=array(NA,dim=c(nmem,nobs))
  deltabc=array(NA,dim=c(nmem,nobs))

  if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA") 
  {
      gpp=array(NA,dim=c(nobs))
      reco=array(NA,dim=c(nobs))
  }

  if(bioflux_flag == "CT_OPT")
      bio=array(NA,dim=c(nobs))

  ocn=array(NA,dim=c(nobs))
  fos=array(NA,dim=c(nobs))
  fir=array(NA,dim=c(nobs))
  bou=array(NA,dim=c(nobs))
  
  #mv <- 0 # missing value #ccgg_evn
  var_ccgg_evn <- ncvar_def( name="obs_num", units=" ", dim=list(obsdim), longname="NOAA_CCGG_Event_Number", prec="integer" )
  var_flask <- ncvar_def( name="flask", units="mol tracer (mol air)^-1", dim=list(memdim,obsdim), longname="mole_fraction_of_trace_gas_in_air", prec="float" )
  var_deltacarbon <- ncvar_def( name="deltaco2", units="mol tracer (mol air)^-1", dim=list(memdim,obsdim),   longname="delta_co2_air", prec="float" )
  var_deltabc <- ncvar_def( name="deltabc", units="mol tracer (mol air)^-1", dim=list(memdim,obsdim), longname="delta_boundary_conditions", prec="float" )

 if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
 {
     var_gpp <- ncvar_def( name="gpp", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from gpp", prec="float" )
     var_reco <- ncvar_def( name="reco", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from ecosystem respiration", prec="float" )
 }

if(bioflux_flag == "CT_OPT")
     var_bio <- ncvar_def( name="bio", units="mol (mol air)^-1", dim=list(obsdim), longname="simulation using CT optimized bio flux as a reference", prec="float" )

var_ocn <- ncvar_def( name="ocn", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from ocean", prec="float" )
  var_fos <- ncvar_def( name="fossil", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from fossil", prec="float" )
  var_fir <- ncvar_def( name="fire", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from fire", prec="float" )
  var_bou <- ncvar_def( name="boundary", units="mol (mol air)^-1", dim=list(obsdim), longname="boundary concentration", prec="float" )


  for(i in 1:nobs)
  { 
    if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
    {
        ccgg_evn[i]=vals[i,1]
        gpp[i]=vals[i,3*nmem+2]
        reco[i]=vals[i,3*nmem+3]
        ocn[i]=vals[i,3*nmem+4]
        fos[i]=vals[i,3*nmem+5]
        fir[i]=vals[i,3*nmem+6]
        bou[i]=vals[i,3*nmem+7]
    }
    if(bioflux_flag == "CT_OPT")
    {
       ccgg_evn[i]=vals[i,1]
       bio[i]=vals[i,3*nmem+2]
       ocn[i]=vals[i,3*nmem+3]
       fos[i]=vals[i,3*nmem+4]
       fir[i]=vals[i,3*nmem+5]
       bou[i]=vals[i,3*nmem+6]
    }
  }
  
  for(j in 1:nmem)
    for(i in 1:nobs)
    {
      flask[j,i]=vals[i,j+1]
      deltacarbon[j,i]=vals[i,nmem+j+1]
      deltabc[j,i]=vals[i,2*nmem+j+1]
    }

print(deltabc)

if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
    ncid_new <- nc_create( output_fname, list(var_ccgg_evn,var_flask,var_deltacarbon,var_deltabc,var_gpp,var_reco,var_ocn,var_fos,var_fir,var_bou))

if(bioflux_flag == "CT_OPT")
     ncid_new <- nc_create( output_fname, list(var_ccgg_evn,var_flask,var_deltacarbon,var_deltabc,var_bio,var_ocn,var_fos,var_fir,var_bou))

  ncvar_put( ncid_new,var_ccgg_evn, ccgg_evn, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_flask, flask, start=c(1,1), count=c(nmem,nobs))
  ncvar_put( ncid_new,var_deltacarbon, deltacarbon, start=c(1,1), count=c(nmem,nobs))
  ncvar_put( ncid_new,var_deltabc, deltabc, start=c(1,1), count=c(nmem,nobs))

  if(bioflux_flag == "SiB3" || bioflux_flag == "SiBCASA")
  {
      ncvar_put( ncid_new,var_gpp, gpp, start=c(1), count=c(nobs))
      ncvar_put( ncid_new,var_reco, reco, start=c(1), count=c(nobs))
  }
  if(bioflux_flag == "CT_OPT")
      ncvar_put( ncid_new,var_bio, bio, start=c(1), count=c(nobs))

  ncvar_put( ncid_new,var_ocn, ocn, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_fos, fos, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_fir, fir, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_bou, bou, start=c(1), count=c(nobs))
  ncatt_put( ncid_new, 0, "Institution", "Centre for Isotope research, University of Groningen")
  ncatt_put( ncid_new, 0, "Contact", "wei.he@rug.nl")
  ncatt_put( ncid_new, 0, "Source", "CarbonTracker-WRF-STILT 1.0")      

  date = format(Sys.time(), "%b %d, %Y")
  user = Sys.getenv("LOGNAME")
  history = paste("created on",date,"by",user,sep=" ")
  ncatt_put( ncid_new, 0, "History", history)

  nc_close( ncid_new )    
}
