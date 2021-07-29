"""CarbonTracker Data Assimilation Shell (CTDAS) Copyright (C) 2017 Wouter Peters. 
Users are recommended to contact the developers (wouter.peters@wur.nl) to receive
updates of the code. See also: http://www.carbontracker.eu. 

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation, 
version 3. This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with this 
program. If not, see <http://www.gnu.org/licenses/>."""
#!/usr/bin/env python
import sys
sys.path.append('../../')
import os
import numpy as np
import string
import datetime as dt
import logging
import re
import da.tools.io4 as io
import netCDF4 as cdf
import matplotlib.pyplot as plt 
import matplotlib
#from maptools import *
from matplotlib.font_manager import FontProperties

fontsize = 10

def nice_lat(cls,format='html'):
    #
    # Convert latitude from decimal to cardinal
    #
    if cls > 0:
        h = 'N'
    else:
        h = 'S'
  
    dec, deg = np.math.modf(cls)

    #return string.strip('%2d %2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))
    if format == 'python':
        return string.strip('%3d$^\circ$%2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))
    if format == 'html':
        return string.strip('%3d&deg%2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))

def nice_lon(cls,format='html'):
    #
    # Convert longitude from decimal to cardinal
    #
    if cls > 0:
        h = 'E'
    else:
        h = 'W'
  
    dec, deg = np.math.modf(cls)

    #return string.strip('%3d %2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))
    if format == 'python':
        return string.strip('%3d$^\circ$%2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))
    if format == 'html':
        return string.strip('%3d&deg%2d\'%s' % (abs(deg), round(abs(60 * dec), 0), h))

def nice_alt(cls):
    #
    # Reformat elevation or altitude
    #
    #return string.strip('%10.1f masl' % round(cls, -1))
    return string.strip('%i masl' %cls)


def summarize_obs(analysisdir, printfmt='html'):
    """***************************************************************************************
    Call example:

    python summarize_obs.py 

    Option printfmt    : [tex,scr,html] print summary table in latex, terminal, or html format 

    Other options are all those needed to create a dacycle object

    OR:

    call directly from a python script as:

    q=summarize_obs(dacycle,printfmt='html')

    ***************************************************************************************"""

    sumdir = os.path.join(analysisdir, 'summary')
    if not os.path.exists(sumdir):
        logging.info("Creating new directory " + sumdir)
        os.makedirs(sumdir)

    mrdir = os.path.join(analysisdir, 'data_molefractions')
    if not os.path.exists(mrdir):
        logging.error("Input directory does not exist (%s), exiting... " % mrdir)
        return None

    mrfiles = os.listdir(mrdir)
    infiles = [os.path.join(mrdir, f) for f in mrfiles if f.endswith('.nc')]

    if printfmt == 'tex': 
        print('\\begin{tabular*}{\\textheight}{l l l l r r r r}')
        print('Code &  Name & Lat, Lon, Elev & Lab &  N (flagged) & $\\sqrt{R}$  &Inn \\XS &Bias\\\\')
        print('\hline\\\\ \n\multicolumn{8}{ c }{Semi-Continuous Surface Samples}\\\\[3pt] ')
        fmt =('%8s  & ' + ' %55s  & ' + '%20s &' + '%6s &' + ' %4d (%d)  & ' + ' %5.2f  & ' + ' %5.2f & ' + '%+5.2f  \\\\')
    elif printfmt == 'html':
        tablehead = \
              "<TR>\n <TH> Site code </TH> \
                   <TH> Sampling Type </TH> \
                   <TH> Lab. </TH> \
                   <TH> Country </TH> \
                   <TH> Lat, Lon, Elev. (m ASL) </TH> \
                   <TH> No. Obs. Available </TH>  \
                   <TH> No. Obs. Assimilated </TH>  \
                   <TH> &#8730;R (&mu;mol mol<sup>-1</sup>) </TH>  \
                   <TH> &#8730;HPH  (&mu;mol mol<sup>-1</sup>) </TH> \
                   <TH> H(x)-y (&mu;mol mol<sup>-1</sup>) </TH> \
                   <TH> H(x)-y (JJAS) (&mu;mol mol<sup>-1</sup>) </TH>  \
                   <TH> H(x)-y (NDJFMA) (&mu;mol mol<sup>-1</sup>) </TH>  \
                   <TH> Inn. &#935;<sup>2</sup></TH> \
                   <TH> Site code </TH>\n \
               </TR>\n"

        fmt = """<TR> \n \
                <TD><a href='javascript:LoadCO2Tseries("%s")'>%s </a></TD>\
                <TD>%s</TD>\
                <TD>%s</TD>\
                <TD>%40s</TD>\
                <TD>%s</TD>\
                <TD>%d</TD>\
                <TD>%d</TD>\
                <TD>%+5.2f</TD>\
                <TD>%+5.2f</TD>\
                <TD>%+5.2f&plusmn;%5.2f</TD>\
                <TD>%+5.2f&plusmn;%5.2f</TD>\
                <TD>%+5.2f&plusmn;%5.2f</TD>\
                <TD bgcolor=%s>%+5.2f</TD>\
                <TD>%s</TD>\n \
               </TR>\n"""
    elif printfmt == 'scr':
        print('Code   Site     NObs flagged  R  Inn X2')
        fmt = '%8s ' + ' %55s  %s %s' + ' %4d ' + ' %4d ' + ' %5.2f ' + ' %5.2f'

    table = []

    for infile in infiles:
            #if not 'mlo_surface-insitu' in infile: continue
            #if not 'poc' in infile: continue
            logging.debug( infile )
            f = io.CT_CDF(infile, 'read')
            date = f.get_variable('time')
            obs = f.get_variable('value') * 1e6
            mdm = f.get_variable('modeldatamismatch') * 1e6
            simulated_fc = f.get_variable('modelsamplesmean_forecast') * 1e6
            simulated = f.get_variable('modelsamplesmean') * 1e6
            simulated_std = f.get_variable('modelsamplesstandarddeviation_forecast') * 1e6
            hphtr = f.get_variable('totalmolefractionvariance_forecast') * 1e6 * 1e6
            flag = f.get_variable('flag_forecast')
            obs_avail = len(np.ma.compressed(mdm))

            pydates = np.array([dt.datetime(1970, 1, 1) + dt.timedelta(seconds=int(d)) for d in date])
 
            sampled = (np.ma.getmaskarray(simulated) == False)
 
            pydates = pydates.compress(sampled)
            simulated = simulated.compress(sampled)
            simulated_fc = simulated_fc.compress(sampled)
            obs = obs.compress(sampled)
            mdm = mdm.compress(sampled)
            hphtr = hphtr.compress(sampled)
            flag = flag.compress(sampled)
            
            if f.site_code.upper() == 'POC': 
                pydates = pydates.compress(simulated > -9000)
                simulated_fc = simulated_fc.compress(simulated > -9000)
                obs = obs.compress(simulated > -9000)
                mdm = mdm.compress(simulated > -9000)
                hphtr = hphtr.compress(simulated > -9000)
                flag = flag.compress(simulated > -9000)
                simulated = simulated.compress(simulated > -9000)
            
            if mdm.mean() > 900:
                pydates = pydates.compress(flag != 2)
                simulated_fc = simulated_fc.compress(flag != 2)
                simulated = simulated.compress(flag != 2)
                obs = obs.compress(flag != 2)
                mdm = mdm.compress(flag != 2)
                hphtr = hphtr.compress(flag != 2)
                obs_assim = 0
            else:
                pydates = pydates.compress(flag == 0)
                simulated_fc = simulated_fc.compress(flag == 0)
                simulated = simulated.compress(flag == 0)
                obs = obs.compress(flag == 0)
                mdm = mdm.compress(flag == 0)
                hphtr = hphtr.compress(flag == 0)
                obs_assim = len(np.ma.compressed(mdm))
                 
            summer = [i for i, d in enumerate(pydates) if d.month in [6, 7, 8, 9] ]
            winter = [i for i, d in enumerate(pydates) if d.month in [11, 12, 1, 2, 3, 4] ]
            
            diff = ((simulated - obs).mean())
            diffsummer = ((simulated - obs).take(summer).mean())
            diffwinter = ((simulated - obs).take(winter).mean())
            diffstd = ((simulated - obs).std())
            diffsummerstd = ((simulated - obs).take(summer).std())
            diffwinterstd = ((simulated - obs).take(winter).std())
            #chi_summer = ((simulated - obs)**2/mdm).take(summer).mean() 
            #chi_winter = ((simulated - obs)**2/mdm).take(winter).mean() 
            #n_summer = simulated.take(summer).shape[0]
            #n_winter = simulated.take(winter).shape[0]
            #print 'summer: %0.2f, %0.2f, %0.2f, %i'%(diffsummer,diffsummerstd,chi_summer,n_summer)
            #print 'winter: %0.2f, %0.2f, %0.2f, %i'%(diffwinter,diffwinterstd,chi_winter,n_winter)
            chi_sq = -99
            if mdm.mean() < 900:
                chi_sq = ((simulated_fc - obs)**2/hphtr).mean() 
                #chi_sq = ((simulated - obs)**2/mdm).mean() 
            if mdm.mean() > 900:
                chi_clr = '#EEEEEE'
            elif chi_sq > 1.2:
                chi_clr = '#ff0000'
            elif chi_sq < 0.5:
                chi_clr = '#ff7f00'
            else: chi_clr = '#00cc00'

            if abs(f.site_latitude) > 9999:
                location = 'Variable'
            else:location = nice_lat(f.site_latitude,'html') + ', ' + nice_lon(f.site_longitude,'html') + ', ' + nice_alt(f.site_elevation)
            if 'site_country' in f.ncattrs():
                country = f.site_country
            else: country = 'Multiple'    

            if printfmt == 'html':
                ss = (f.dataset_name[4:],
                    f.site_code.upper(),
                    f.dataset_project,
                    f.lab_1_abbr,
                    country,
                    location,
                    obs_avail,
                    obs_assim,
                    mdm.mean(),
                    np.sqrt((simulated_std ** 2).mean()),
                    diff, diffstd,
                    diffsummer, diffsummerstd,
                    diffwinter, diffwinterstd,
                    chi_clr, chi_sq,
                    f.site_code.upper())
            
            table.append(ss)
            f.close()

    if printfmt == 'tex':
        saveas = os.path.join(sumdir, 'site_table.tex')
        f = open(saveas, 'w')
    elif printfmt == 'html':
        saveas = os.path.join(sumdir, 'site_table.html')
        f = open(saveas, 'w')
        txt = "<meta http-equiv='content-type' content='text/html;charset=utf-8' />\n"
        f.write(txt)
        txt = "<table border=1 cellpadding=2 cellspacing=2 width='100%' bgcolor='#EEEEEE'>\n"
        f.write(txt)

    f.write(tablehead)

    for i, ss in enumerate(table):

        f.write(fmt % ss)
        if (i + 1) % 15 == 0:
            f.write(tablehead)

    if printfmt == 'tex': 
        f.write('\cline{2-8}\\\\')
        f.write('\hline \\\\')
        f.write('\end{tabular*}')
    else:
        txt = "\n</table>"
        f.write(txt)
    f.close()

    logging.info("File written with summary: %s" % saveas)

def make_map(analysisdir): #makes a map of amount of assimilated observations per site

    sumdir = os.path.join(analysisdir, 'summary')
    if not os.path.exists(sumdir):
        logging.info("Creating new directory " + sumdir)
        os.makedirs(sumdir)

    mrdir = os.path.join(analysisdir, 'data_molefractions')
    if not os.path.exists(mrdir):
        logging.error("Input directory does not exist (%s), exiting... " % mrdir)
        return None

    mrfiles = os.listdir(mrdir)
    infiles = [os.path.join(mrdir, f) for f in mrfiles if f.endswith('.nc')]
   
    lats=[]
    lons=[]
    labs=[]
    nobs=[]
    for files in infiles:
        f=cdf.Dataset(files)
        if f.variables['modeldatamismatch'][:].max() < 0.001:
            sim = f.variables['modelsamplesmean'][:]
            flag = f.variables['flag_forecast'][:]
            sim = sim.compress(flag != 2)
            sampled = (np.ma.getmaskarray(sim) == False)
            sim = sim.compress(sampled)
            lats.append(f.site_latitude)
            lons.append(f.site_longitude)
            labs.append(f.site_code)
            nobs.append(len(sim))
        f.close()
    
    lats = np.array(lats)
    lons = np.array(lons)
    labs = np.array(labs)
    nobs = np.array(nobs)

    saveas = os.path.join(sumdir, 'networkmap')
    logging.info("Making map: %s" % saveas)

    fig = plt.figure(1,figsize=(20,12))
    ax = fig.add_axes([0.05,0.1,0.9,0.8])
    #m,nx,ny = select_map('Global Cylinder')
    m,nx,ny = select_map('Europe Conformal')
    m.drawcountries()
    m.drawcoastlines()
    parallels = arange(-90.,91,30.)
    m.drawparallels(parallels,color='grey',linewidth=0.5,dashes=[1,0.001],labels=[1,0,0,1],fontsize=16)
    meridians = arange(-180.,181.,60.)
    m.drawmeridians(meridians,color='grey',linewidth=0.5,dashes=[1,0.001],labels=[1,0,0,1],fontsize=16)

    #for lon,lat,name,n in zip(lons,lats,names,nobs):
    count = 0
    for i in range(len(lats)):
        if nobs[i]   < 250: 
            n = 0
            c = 'blue'
        elif nobs[i] < 500:
            n = 1
            c = 'green'
        elif nobs[i] < 750: 
            n = 2
            c = 'orange'
        elif nobs[i] < 1000: 
            n = 3
            c = 'brown'
        else: 
            n = 4
            c = 'red'
        if lons[i] > -900:
            x,y = m(lons[i],lats[i])
            ax.plot(x,y,'o',color=c,markersize=12+1.5*n)#,markeredgecolor='k',markeredgewidth=2)
            #ax.annotate(labs[i],xy=m(lons[i],lats[i]),xycoords='data',fontweight='bold')
        else:    
            x,y = m(169,87-count)
            ax.plot(x,y,'o',color=c,markersize=12+1.5*n)
            ax.annotate(labs[i],xy=m(172,86-count),xycoords='data',fontweight='bold')
            count = count + 4

    fig.text(0.15,0.945,u'\u2022',fontsize=35,color='blue')
    fig.text(0.16,0.95,': N<250',fontsize=24,color='blue')
    fig.text(0.30,0.94,u'\u2022',fontsize=40,color='green')
    fig.text(0.31,0.95,': N<500',fontsize=24,color='green')
    fig.text(0.45,0.94,u'\u2022',fontsize=45,color='orange')
    fig.text(0.46,0.95,': N<750',fontsize=24,color='orange')
    fig.text(0.60,0.939,u'\u2022',fontsize=50,color='brown')
    fig.text(0.61,0.95,': N<1000',fontsize=24,color='brown')
    fig.text(0.75,0.938,u'\u2022',fontsize=55,color='red')
    fig.text(0.765,0.95,': N>1000',fontsize=24,color='red')
    ax.set_title('Assimilated observations',fontsize=24)

    font0= FontProperties(size=15,style='italic',weight='bold')
    txt='CarbonTracker Europe\n $\copyright$ Wageningen University'
    clr='green'
    fig.text(0.82,0.01,txt,ha='left',font_properties = font0, color=clr )
    saveas=os.path.join(sumdir,'networkmap.png')
    fig.savefig(saveas,dpi=200)
    saveas=os.path.join(sumdir,'networkmap.large.png')
    fig.savefig(saveas,dpi=300)
    close(fig)

def summarize_stats(dacycle):
    """
    Summarize the statistics of the observations for this cycle
    This includes X2 statistics, RMSD, and others for both forecast and
    final fluxes
    """
    

    sumdir = os.path.join(dacycle['dir.analysis'], 'summary')
    if not os.path.exists(sumdir):
        logging.info("Creating new directory " + sumdir)
        os.makedirs(sumdir)

    # get forecast data from optimizer.ddddd.nc

    startdate = dacycle['time.start'] 
    dacycle['time.sample.stamp'] = "%s" % (startdate.strftime("%Y%m%d"),)
    infile = os.path.join(dacycle['dir.output'], 'optimizer.%s.nc' % dacycle['time.sample.stamp'])

    if not os.path.exists(infile):
        logging.error("File not found: %s" % infile)
        raise IOError

    f = io.CT_CDF(infile, 'read')
    sites = f.get_variable('sitecode')
    y0 = f.get_variable('observed') * 1e6
    hx = f.get_variable('modelsamplesmean_prior') * 1e6
    dF = f.get_variable('modelsamplesdeviations_prior') * 1e6
    HPHTR = f.get_variable('totalmolefractionvariance').diagonal() * 1e6 * 1e6
    R = f.get_variable('modeldatamismatchvariance').diagonal() * 1e6 * 1e6
    flags = f.get_variable('flag')
    f.close()

    HPHT = dF.dot(np.transpose(dF)).diagonal() / (dF.shape[1] - 1.0)
    rejected = (flags == 2.0)

    sitecodes = [string.join(s.compressed(), '').strip() for s in sites]


    # calculate X2 per observation for this time step

    x2 = []
    for i, site in enumerate(sitecodes):
        
        x2.append((y0[i] - hx[i]) ** 2 / HPHTR[i])

    x2 = np.ma.masked_where(HPHTR == 0.0, x2)

    # calculate X2 per site
    saveas = os.path.join(sumdir, 'x2_table_%s.html' % dacycle['time.sample.stamp'])
    logging.info("Writing HTML tables for this cycle (%s)" % saveas)
    f = open(saveas, 'w')
    txt = "<meta http-equiv='content-type' content='text/html;charset=utf-8' />\n"
    f.write(txt)
    txt = "<table border=1 cellpadding=2 cellspacing=2 width='100%' bgcolor='#EEEEEE'>\n"
    f.write(txt)
    tablehead = \
          "<TR>\n <TH> Site code </TH> \
               <TH> N<sub>obs</sub> </TH>  \
               <TH> N<sub>rejected</sub> </TH>  \
               <TH> &#8730;R (&mu;mol mol<sup>-1</sup>) </TH>  \
               <TH> &#8730;HPH<sup>T</sup> (&mu;mol mol<sup>-1</sup>) </TH>  \
               <TH> H(x)-y (&mu;mol mol<sup>-1</sup>) </TH> \n \
               <TH> X2 </TH> \n \
           </TR>\n"

    fmt = """<TR> \n \
            <TD>%s</TD>\
            <TD>%d</TD>\
            <TD>%d</TD>\
            <TD>%+5.2f</TD>\
            <TD>%+5.2f</TD>\
            <TD>%+5.2f&plusmn;%5.2f</TD>\
            <TD>%5.2f</TD>\n \
           </TR>\n"""

    f.write(tablehead)

    set_sites = set(sitecodes)
    set_sites = np.sort(list(set_sites))

    for i, site in enumerate(set_sites):
        sel = [i for i, s in enumerate(sitecodes) if s == site]
        ss = (site, len(sel), rejected.take(sel).sum(), np.sqrt(R.take(sel)[0]), np.sqrt(HPHT.take(sel).mean()), (hx - y0).take(sel).mean(), (hx - y0).take(sel).std(), x2.take(sel).mean(),)
        #print site,sel,x2.take(sel)

        f.write(fmt % ss)
        if (i + 1) % 15 == 0:
            f.write(tablehead)

    txt = "\n</table>"
    f.write(txt)
    f.close()

    # Now summarize for each site across time steps

    if not dacycle['time.start'] >= dt.datetime(2008, 12, 29):
        return

    logging.info("Writing HTML tables for each site")
    for site in set_sites:
        saveas = os.path.join(sumdir, '%s_x2.html' % site)
        f = open(saveas, 'w')
        logging.debug(saveas)
        txt = "<meta http-equiv='content-type' content='text/html;charset=utf-8' />\n"
        f.write(txt)
        txt = "<table border=1 cellpadding=2 cellspacing=2 width='100%' bgcolor='#EEEEEE'>\n"
        f.write(txt)
        tablehead = \
              "<TR>\n <TH> From File </TH> \
                   <TH> Site </TH>  \
                   <TH> N<sub>obs</sub> </TH>  \
                   <TH> N<sub>rejected</sub> </TH>  \
                   <TH> &#8730;R (&mu;mol mol<sup>-1</sup>) </TH>  \
                   <TH> &#8730;HPH<sup>T</sup> (&mu;mol mol<sup>-1</sup>) </TH>  \
                   <TH> H(x)-y (&mu;mol mol<sup>-1</sup>) </TH> \n \
                   <TH> X2 </TH> \n \
               </TR>\n"
        f.write(tablehead)

        files = os.listdir(sumdir)
        x2_files = [fil for fil in files if fil.startswith('x2')]
        for htmlfile in x2_files:
            lines = grep(site, os.path.join(sumdir, htmlfile))
            for line in lines:
                f.write('<TR>\n')
                f.write('<TD>' + htmlfile + '</TD>')
                f.write(line + '\n')
                f.write('</TR>\n')

        txt = "\n</table>"
        f.write(txt)
        f.close()


def grep(pattern, fil):
    fileObj = open(fil, 'r')
    r = []
    for line in fileObj:
        if re.search(pattern, line):
            r.append(line)
    return r

# main body if called as script

if __name__ == '__main__':    # started as script

#    sys.path.append('../../')

    logging.root.setLevel(logging.DEBUG)
    analysisdir = "/store/empa/em05/parsenov/output"
    
    summarize_obs(analysisdir)
    #make_map(analysisdir)

    sys.exit(0)


