#******************************************************************************************
# Function for writing a netcdf4 file
# created by W.He on April 15, 2015 
#****************************************************************************************** 

write.results.netcdf <- function(vals,nobs,nmem,output_fname) { 
    
  obsdim <- ncdim_def( 'obsnum', '', 1:nobs, unlim=FALSE, create_dimvar=FALSE ) #IMPORTANT TO USE create_dimvar
  memdim <- ncdim_def( 'nmembers', '', 1:nmem, unlim=FALSE, create_dimvar=FALSE)
  
  #obs_num=array(NA,dim=c(nobs))
  ccgg_evn=array(NA,dim=c(nobs))
  flask=array(NA,dim=c(nmem,nobs))
  gpp=array(NA,dim=c(nobs))
  reco=array(NA,dim=c(nobs))
  ocn=array(NA,dim=c(nobs))
  fos=array(NA,dim=c(nobs))
  fir=array(NA,dim=c(nobs))
  bou=array(NA,dim=c(nobs))
  bio=array(NA,dim=c(nobs))

  #mv <- 0 # missing value #ccgg_evn
  #var_obs_num <- ncvar_def( name="obs_num", dim=list(obsdim), longname="Unique_Dataset_observation_index_number", units=" ", prec="integer" )
  var_ccgg_evn <- ncvar_def( name="obs_num", units=" ", dim=list(obsdim), longname="NOAA_CCGG_Event_Number", prec="integer" )
  var_flask <- ncvar_def( name="flask", units="mol tracer (mol air)^-1", dim=list(memdim,obsdim), longname="mole_fraction_of_trace_gas_in_air", prec="float" )
  
  var_gpp <- ncvar_def( name="gpp", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from gpp", prec="float" )
  var_reco <- ncvar_def( name="reco", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from ecosystem respiration", prec="float" )
  var_ocn <- ncvar_def( name="ocn", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from ocean", prec="float" )
  var_fos <- ncvar_def( name="fossil", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from fossil", prec="float" )
  var_fir <- ncvar_def( name="fire", units="mol (mol air)^-1", dim=list(obsdim), longname="contribution from fire", prec="float" )
  var_bou <- ncvar_def( name="boundary", units="mol (mol air)^-1", dim=list(obsdim), longname="boundary concentration", prec="float" )
  var_bio <- ncvar_def( name="bio", units="mol (mol air)^-1", dim=list(obsdim), longname="simulation using CT optimized bio flux as a reference", prec="float" )

  for(i in 1:nobs)
  { 
    #obs_num[i]=i
    ccgg_evn[i]=vals[i,1]
    gpp[i]=vals[i,nmem+2]
    reco[i]=vals[i,nmem+3]
    ocn[i]=vals[i,nmem+4]
    fos[i]=vals[i,nmem+5]
    fir[i]=vals[i,nmem+6]
    bou[i]=vals[i,nmem+7]
    bio[i]=vals[i,nmem+8]
  }
  
  for(j in 1:nmem)
    for(i in 1:nobs)
      flask[j,i]=vals[i,j+1]
  
  ncid_new <- nc_create( output_fname, list(var_ccgg_evn,var_flask,var_gpp,var_reco,var_ocn,var_fos,var_fir,var_bou,var_bio))
  ncvar_put( ncid_new,var_ccgg_evn, ccgg_evn, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_flask, flask, start=c(1,1), count=c(nmem,nobs))
  ncvar_put( ncid_new,var_gpp, gpp, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_reco, reco, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_ocn, ocn, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_fos, fos, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_fir, fir, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_bou, bou, start=c(1), count=c(nobs))
  ncvar_put( ncid_new,var_bou, bio, start=c(1), count=c(nobs))
  ncatt_put( ncid_new, 0, "Institution", "Centre for Isotope research, University of Groningen")
  ncatt_put( ncid_new, 0, "Contact", "wei.he@rug.nl")
  ncatt_put( ncid_new, 0, "Source", "CarbonTracker-STILT 1.0")      

  date = format(Sys.time(), "%b %d, %Y")
  user = Sys.getenv("LOGNAME")
  history = paste("created on",date,"by",user,sep=" ")
  ncatt_put( ncid_new, 0, "History", history)

  nc_close( ncid_new )    
}
