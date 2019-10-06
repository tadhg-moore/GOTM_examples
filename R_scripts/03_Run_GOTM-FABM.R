#First set your working directory to the lake folder
setwd("...\\feeagh") 

# Load libraries
library(GOTMr)
library(gotmtools)
library(ggplot2)
library(ggpubr)
library(reshape)

#Files
gotm_yaml <- 'gotm.yaml' # GOTM config file
obs.file <- 'wtemp_dly_2010_2015.dat' # Temperature profile file
met_file <- 'met_2010_2015.dat' # path to file with time series
fabm_yaml <- 'fabm.yaml' # FABM config file - you can view in RStudio


start <- '2015-01-01 00:00:00' # start date and time [yyyy-mm-dd HH:MM:SS;]
stop <- '2015-12-31 00:00:00' # stop date and time [yyyy-mm-dd HH:MM:SS;]

###### Input values into 'gotm.yaml'
input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)

#Input calibrated parameters
input_yaml(file = gotm_yaml, label = 'heat', key = 'scale_factor', value = 1.1)
input_yaml(file = gotm_yaml, label = 'swr', key = 'scale_factor', value = 1.2)
input_yaml(file = gotm_yaml, label = 'turb_param', key = 'k_min', value = 5.87e-6)
input_yaml(file = gotm_yaml, label = 'g1', key = 'constant_value', value = 0.17)
input_yaml(file = gotm_yaml, label = 'g2', key = 'constant_value', value = 3)

# Switch on FABM
input_yaml(file = gotm_yaml, label = 'fabm', key = 'use', value = 'true')
input_yaml(file = gotm_yaml, label = 'fabm', key = 'repair_state', value = 'true')

#### Create inital temp file
init_prof(obs_file = obs.file, date = start, tprof_file = 'init_tprof.dat')

# Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 

#### View NetCDF file
out <- 'output.nc' # Output from GOTM

#### Plot water temperature
p1 <- plot_wtemp(out, size =2, zlab = '°C')+
  ggtitle('water temperature')

#### Load in observed data
obs <- load_obs(obs.file = obs.file)
obs <- obs[(obs$date >= start & obs$date <= stop),] # Subset to same time as modelled
mod_wtemp <- get_vari(ncdf = out, var = 'temp') # Extract potential temp
z <- get_vari(out, 'z') # Extract depths corresponding to temp
mod <- setmodDepths(mod.val = mod_wtemp, mod.dep = z, obs = obs) # Match time and depths to obs data

#### Model diagnostic plot
g1 <- diag_plots(mod = mod, obs = obs, colourblind = F)
g1

#View list of variables
list_vars(out, long = T)
list_vars(out, long = F)

#### Extract and plot total Chl-a
p2 <- plot_vari(out, "total_chlorophyll_calculator_result", size = 2, zlab = "mg chl a/m3")+
  ggtitle("total_chlorophyll_calculator result")
p2

#### Extract and plot Oxygen
p3 <- plot_vari(out, "selma_o2" , size = 2, zlab = "mmol O2/m3")+
  ggtitle("selma oxygen")
p3

#### Plot all variables
g1 <- ggarrange(p1,p2,p3, nrow =3, align = 'v')
g1

