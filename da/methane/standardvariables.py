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
standard_variables = { 'bio_flux_prior' : {'name'        : 'bio_flux_prior',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, terrestrial biosphere, not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_opt' : {'name'          : 'bio_flux_opt',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, terrestrial biosphere , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_prior' : {'name'     : 'fossil_flux_prior',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, anthropogenic , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_opt' : {'name'       : 'fossil_flux_opt',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, anthropogenic , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                        'ocn_flux_imp' : {'name'         : 'ocn_flux_imp',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, open ocean , imposed ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fire_flux_imp' : {'name'         : 'fire_flux_imp',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, biomass burning , imposed ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'term_flux_imp' : {'name'         : 'term_flux_imp',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of methane, termites , imposed ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_methane_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_prior_cov' : {'name'    : 'bio_flux_prior_cov',\
                                         'units'         : '[mol m-2 s-1]^2' ,\
                                         'long_name'     : 'Covariance of surface flux of methane, terrestrial vegetation , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_prior_ensemble' : {'name'    : 'bio_flux_prior_ensemble',\
                                         'units'         : '[mol m-2 s-1]' ,\
                                         'long_name'     : 'Ensemble of surface flux of methane, terrestrial vegetation , not optimized ', \
                                         'comment'       : "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance", \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_opt_cov' : {'name'      : 'bio_flux_opt_cov',\
                                         'units'         : '[mol m-2 s-1]^2' ,\
                                         'long_name'     : 'Covariance of surface flux of methane, terrestrial vegetation , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_opt_ensemble' : {'name'    : 'bio_flux_opt_ensemble',\
                                         'units'         : '[mol m-2 s-1]' ,\
                                         'long_name'     : 'Ensemble of surface flux of methane, terrestrial vegetation , optimized ', \
                                         'comment'       : "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance", \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_prior_cov' : {'name'    : 'fossil_flux_prior_cov',\
                                         'units'         : '[mol m-2 s-1]^2' ,\
                                         'long_name'     : 'Covariance of surface flux of methane, anthropogenic , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_prior_ensemble' : {'name'    : 'fossil_flux_prior_ensemble',\
                                         'units'         : '[mol m-2 s-1]' ,\
                                         'long_name'     : 'Ensemble of surface flux of methane, anthropogenic , not optimized ', \
                                         'comment'       : "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance", \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_opt_cov' : {'name'      : 'fossil_flux_opt_cov',\
                                         'units'         : '[mol m-2 s-1]^2' ,\
                                         'long_name'     : 'Covariance of surface flux of methane, anthropogenic , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_opt_ensemble' : {'name'    : 'fossil_flux_opt_ensemble',\
                                         'units'         : '[mol m-2 s-1]' ,\
                                         'long_name'     : 'Ensemble of surface flux of methane, anthropogenic , optimized ', \
                                         'comment'       : "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance", \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'decimal_date' :  {'name'         : 'decimal_date',\
                                         'units'         : 'years' ,\
                                         'long_name'     : 'dates and times', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'date', \
                                         'dims'          : (), \
                                         'dtype'         : 'double', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'date' :         {'name'          : 'date',\
                                         'units'         : 'days since 2000-01-01 00:00:00 UTC' ,\
                                         'long_name'     : 'UTC dates and times', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'date', \
                                         'dims'          : (), \
                                         'dtype'         : 'double', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'idate' :        {'name'          : 'idate',\
                                         'units'         : 'yyyy MM dd hh mm ss ' ,\
                                         'long_name'     : 'integer components of date and time', \
                                         'standard_name' : 'calendar_components', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'dims'          : (), \
                                         'dtype'         : 'int', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'latitude' :     {'name'          : 'latitude',\
                                         'units'         : 'degrees_north ' ,\
                                         'long_name'     : 'latitude', \
                                         'standard_name' : 'latitude', \
                                         'comment'       : 'center of interval',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'longitude' :     {'name'         : 'longitude',\
                                         'units'         : 'degrees_east ' ,\
                                         'long_name'     : 'longitude', \
                                         'standard_name' : 'longitude', \
                                         'comment'       : 'center of interval',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'height' :        {'name'         : 'height',\
                                         'units'         : 'masl ' ,\
                                         'long_name'     : 'height_above_ground_level', \
                                         'standard_name' : 'height_above_ground_level', \
                                         'comment'       : 'value is meters above sea level',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'cell_area' :     {'name'         : 'cell_area',\
                                         'units'         : 'm2 ' ,\
                                         'long_name'     : 'horizontal ara of a gridcell', \
                                         'standard_name' : 'cell_area', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ch4' :           {'name'         : 'ch4',\
                                         'units'         : 'micromol mol-1 ' ,\
                                         'long_name'     : 'mole_fraction_of_methane_in_air', \
                                         'standard_name' : 'mole_fraction_of_methane_in_air', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'meanstate' :     {'name'         : 'statevectormean',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'mean_value_of_state_vector', \
                                         'standard_name' : 'mean_value_of_state_vector', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ensemblestate':  {'name'         : 'statevectorensemble',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'ensemble_value_of_state_vector', \
                                         'standard_name' : 'ensemble_value_of_state_vector', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'meanstate_prior' : {'name'       : 'statevectormean_prior',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'mean_value_of_state_vector_prior', \
                                         'standard_name' : 'mean_value_of_state_vector_prior', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ensemblestate_prior':  {'name'         : 'statevectorensemble_prior',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'ensemble_value_of_state_vector_prior', \
                                         'standard_name' : 'ensemble_value_of_state_vector_prior', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'meanstate_opt' : {'name'       : 'statevectormean_opt',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'mean_value_of_state_vector_optimized', \
                                         'standard_name' : 'mean_value_of_state_vector_opt', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ensemblestate_opt':  {'name'         : 'statevectorensemble_opt',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'ensemble_value_of_state_vector_optimized', \
                                         'standard_name' : 'ensemble_value_of_state_vector_opt', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'unknown' :      {'name'          : '',\
                                         'units'         : '' ,\
                                         'long_name'     : '', \
                                         'standard_name' : '', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                     }




