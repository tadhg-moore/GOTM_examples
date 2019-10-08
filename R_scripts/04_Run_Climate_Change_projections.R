#First set your working directory to the lake folder
setwd("...\\feeagh")

library(GOTMr)
library(gotmtools)
library(ggplot2)
library(ggpubr)
library(reshape)

#install.packages('lubridate')
#install.packages('Hmisc')
library(lubridate)
library(Hmisc)

##### ISIMIP Data
# Example data from the Intersectoral Impact Model Intercomparison Project (ISIMIP: https://www.isimip.org/)


# Meteorological data from Global Circulation Model's (GCM) has been generated under 
# different concentrations of CO2.
# Historical: Historical climate and CO2 concentration. 
# rcp26, #Future climate and CO2 concentration from RCP2.6 	- limit warming to 1.5°C
# rcp60 #Future climate and CO2 concentration from RCP6.0 	- warming up to 3°C




# Copy in config file
file.copy(from = 'isimip_files/gotm_ISIMIP.yaml', to = 'gotm.yaml', overwrite = T)


#Files
gotm_yaml <- 'gotm.yaml' #GOTM config file

dt <- 3600 # time step for integration [s]
depth <- 46.8 # water depth [m]
nlev <- round(depth/0.5)  # number of layers
t_prof_file <- 'init_tprof.dat'  # path to file with series of profiles
met_file <- 'met_file.dat' # path to file with time series
u10_factor <- 1 # scale factor to be applied to values read from file
v10_factor <- 1 # scale factor to be applied to values read from file
at_factor <- 1 # scale factor to be applied to values read from file
ap_factor <- 1 # scale factor to be applied to values read from file
rh_factor <- 1 # scale factor to be applied to values read from file
hum_type <- 1 # humidity metric [1=relative humidity (%), 2=wet-bulb temperature, 3=dew point temperature, 4=specific humidity (kg/kg)]
cc_factor <- 1 # scale factor to be applied to values read from file
sw_factor <- 1 # scale factor to be applied to values read from file
precip_factor <- 1 # scale factor to be applied to values read from file
ice_model <- 0 # model [0=none, 1=Lebedev (1938), 2=MyLake, 3=Winton]
k_min <- 3.6e-6 # minimum turbulent kinetic energy [m^2/s^2]
zeta_method <- 0 ## method [0=constant, 1=from tidal constituents, 2=from file]
zeta_file <- 'lake_level_2010-2016.dat' # path to file with time series
zeta_offset <- -10.723 # offset to be added to values read from file
A <- 0.55 # non-visible fraction of shortwave radiation
g1 <- 0.17 # e-folding depth of non-visible shortwave radiation
g2 <- 3.18 # e-folding depth of visible shortwave radiation
use_fabm <- 'false' # enable FABM 
restart <- 'false' # initialize simulation with state stored in restart.nc 

###### Input values into 'gotm.yaml'
input_yaml(file = gotm_yaml, label = 'time', key = 'dt', value = dt)
input_yaml(file = gotm_yaml, label = 'location', key = 'depth', value = depth)
input_yaml(file = gotm_yaml, label = 'grid', key = 'nlev', value = nlev)
input_yaml(file = gotm_yaml, label = 'temperature', key = 'file', value = t_prof_file)
input_yaml(file = gotm_yaml, label = 'u10', key = 'scale_factor', value = u10_factor)
input_yaml(file = gotm_yaml, label = 'v10', key = 'scale_factor', value = v10_factor)
input_yaml(file = gotm_yaml, label = 'airt', key = 'scale_factor', value = at_factor)
input_yaml(file = gotm_yaml, label = 'airp', key = 'scale_factor', value = ap_factor)
input_yaml(file = gotm_yaml, label = 'hum', key = 'scale_factor', value = rh_factor)
input_yaml(file = gotm_yaml, label = 'hum', key = 'type', value = hum_type)
input_yaml(file = gotm_yaml, label = 'cloud', key = 'scale_factor', value = cc_factor)
input_yaml(file = gotm_yaml, label = 'swr', key = 'scale_factor', value = sw_factor)
input_yaml(file = gotm_yaml, label = 'precip', key = 'scale_factor', value = precip_factor)
input_yaml(file = gotm_yaml, label = 'u10', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'v10', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'airt', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'airp', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'hum', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'cloud', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'swr', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'precip', key = 'file', value = met_file)
input_yaml(file = gotm_yaml, label = 'ice', key = 'model', value = ice_model)
input_yaml(file = gotm_yaml, label = 'turb_param', key = 'k_min', value = k_min)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'method', value = zeta_method)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'file', value = zeta_file)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'offset', value = zeta_offset)
input_yaml(file = gotm_yaml, label = 'A', key = 'constant_value', value = A)
input_yaml(file = gotm_yaml, label = 'g1', key = 'constant_value', value = g1)
input_yaml(file = gotm_yaml, label = 'g2', key = 'constant_value', value = g2)



#### Historical Scenario
# Copy in met file
file.copy(from = 'isimip_files/GFDL_historical_met_file.dat', to = 'met_file.dat', overwrite = T)

scan_timeseries(file = met_file, header = T)
start <- '1965-01-01 00:00:00' # start date and time [yyyy-mm-dd HH:MM:SS;]
stop <- '2005-12-31 00:00:00' # stop date and time [yyyy-mm-dd HH:MM:SS;]

## Open 'init_tprof.dat' and change the date to match start

input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)

# Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 

out <- 'output.nc'

# Extract surface water temperature
wtemp_hist <- get_vari(ncdf = out, 'temp')[,1:2] # Extract just surface temperature
wtemp_hist <- wtemp_hist[(wtemp_hist[,1] >= '1970-01-01'),] # Remove first 5 years 
colnames(wtemp_hist)[2] <- 'surface'
wtemp_hist$yday <- yday(wtemp_hist[,1])

p1 <- ggplot(wtemp_hist, aes(yday, surface))+
  stat_summary(geom = "line", fun.y = mean) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3)+  theme_bw()
p1


#### RCP 2.6 Scenario
# Copy in met file
file.copy(from = 'isimip_files/GFDL_rcp26_met_file.dat', to = 'met_file.dat', overwrite = T)

scan_timeseries(file = met_file, header = T)
start <- "2060-01-01 00:00:00" # start date and time [yyyy-mm-dd HH:MM:SS;]
stop <- "2099-12-31 00:00:00" # stop date and time [yyyy-mm-dd HH:MM:SS;]

## Open 'init_tprof.dat' and change the date to match start

input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)

# Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 

# Extract surface water temperature
wtemp_rcp26 <- get_vari(ncdf = out, 'temp')[,1:2] # Extract just surface temperature
wtemp_rcp26 <- wtemp_rcp26[(wtemp_rcp26[,1] >= '2065-01-01'),] # Remove first 5 years 
colnames(wtemp_rcp26)[2] <- 'surface'
wtemp_rcp26$yday <- yday(wtemp_rcp26[,1])

p2 <- ggplot(wtemp_hist, aes(yday, surface))+
  stat_summary(geom = "line", fun.y = mean, aes(colour = 'Hist')) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3,
               aes(fill = 'Hist'))+
  stat_summary(data = wtemp_rcp26, geom = "line", fun.y = mean, aes(colour = 'RCP2.6')) +
  stat_summary(data = wtemp_rcp26, geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3, aes(fill = 'RCP2.6'))+
  theme_bw()
p2


#### RCP 6.0 Scenario
# Copy in met file
file.copy(from = 'isimip_files/GFDL_rcp60_met_file.dat', to = 'met_file.dat', overwrite = T)

scan_timeseries(file = met_file, header = T)
start <- "2060-01-01 00:00:00" # start date and time [yyyy-mm-dd HH:MM:SS;]
stop <- "2099-12-31 00:00:00" # stop date and time [yyyy-mm-dd HH:MM:SS;]

## Open 'init_tprof.dat' and change the date to match start

input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)

# Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 

# Extract surface water temperature
wtemp_rcp60 <- get_vari(ncdf = out, 'temp')[,1:2] # Extract just surface temperature
wtemp_rcp60 <- wtemp_rcp60[(wtemp_rcp60[,1] >= '2065-01-01'),] # Remove first 5 years 
colnames(wtemp_rcp60)[2] <- 'surface'
wtemp_rcp60$yday <- yday(wtemp_rcp60[,1])

p3 <- ggplot(wtemp_hist, aes(yday, surface))+
  stat_summary(geom = "line", fun.y = mean, aes(colour = 'Hist')) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3,
               aes(fill = 'Hist'))+
  stat_summary(data = wtemp_rcp26, geom = "line", fun.y = mean, aes(colour = 'RCP2.6')) +
  stat_summary(data = wtemp_rcp26, geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3, aes(fill = 'RCP2.6'))+
  stat_summary(data = wtemp_rcp60, geom = "line", fun.y = mean, aes(colour = 'RCP6.0')) +
  stat_summary(data = wtemp_rcp60, geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3, aes(fill = 'RCP6.0'))+
  ylab('Temperature (°C)')+
  xlab('')+
  ggtitle('Climatic Annual Surface Temperature')+
  guides(colour = F, fill = guide_legend(title = 'Scenario'))+
  theme_bw()
p3
