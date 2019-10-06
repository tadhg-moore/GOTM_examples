#First set your working directory to the lake folder
setwd("...\\feeagh")

library(GOTMr)
library(gotmtools)
library(ggplot2)
library(ggpubr)
library(reshape)

#Files
gotm_yaml <- 'gotm.yaml' #GOTM config file
obs.file <- 'wtemp_dly_2010_2015.dat' #Temperature profile file

start <- '2013-01-01 00:00:00' # start date and time [yyyy-mm-dd HH:MM:SS;]
stop <- '2015-12-31 00:00:00' # stop date and time [yyyy-mm-dd HH:MM:SS;]
dt <- 3600 # time step for integration [s]
depth <- 46.8 # water depth [m]
nlev <- round(depth/0.5)  # number of layers
t_prof_file <- 'init_tprof.dat'  # path to file with series of profiles
met_file <- 'met_2010_2015.dat' # path to file with time series
u10_factor <- 1 # scale factor to be applied to values read from file
v10_factor <- 1 # scale factor to be applied to values read from file
at_factor <- 1 # scale factor to be applied to values read from file
ap_factor <- 1 # scale factor to be applied to values read from file
rh_factor <- 1 # scale factor to be applied to values read from file
hum_type <- 3 # humidity metric [1=relative humidity (%), 2=wet-bulb temperature, 3=dew point temperature, 4=specific humidity (kg/kg)]
cc_factor <- 1 # scale factor to be applied to values read from file
sw_factor <- 1 # scale factor to be applied to values read from file
precip_factor <- 1 # scale factor to be applied to values read from file
ice_model <- 0 # model [0=none, 1=Lebedev (1938), 2=MyLake, 3=Winton]
k_min <- 3.6e-6 # minimum turbulent kinetic energy [m^2/s^2]
zeta_method <- 2 # method [0=constant, 1=from tidal constituents, 2=from file]
zeta_file <- 'lake_level_2010-2016.dat' # path to file with time series
zeta_offset <- -10.723 # offset to be added to values read from file
A <- 0.55 # non-visible fraction of shortwave radiation
g1 <- 0.5 # e-folding depth of non-visible shortwave radiation
g2 <- 1.6 # e-folding depth of visible shortwave radiation
use_fabm <- 'false' # enable FABM 
restart <- 'false' # initialize simulation with state stored in restart.nc 


#### View input files
met <- load_input(input_file = met_file, header = T)
windows() # Opens new window for viewing plots
par(mfrow = c(3,3)) # Sets it to 3x3
# Loop through the columns and plot each column
for(i in 2:ncol(met)){
  plot(met$X.date, met[,i], type = 'l', main = colnames(met)[i], xlab = 'Date')
}

plot_hypso(hypsograph_file = 'hypsograph.dat')

#Lake Level
lev <- load_input('lake_level_2010-2016.dat', header = T)
par(mfrow = c(1,1))
plot(lev[,1], lev$level, type = 'l', main = 'Lake Level', ylab = '(m)')


#### Create inital temp file
init_prof(obs_file = obs.file, date = start, tprof_file = 'init_tprof.dat')

###### Input values into 'gotm.yaml'
input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)
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
input_yaml(file = gotm_yaml, label = 'ice', key = 'model', value = ice_model)
input_yaml(file = gotm_yaml, label = 'turb_param', key = 'k_min', value = k_min)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'method', value = zeta_method)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'file', value = zeta_file)
input_yaml(file = gotm_yaml, label = 'zeta', key = 'offset', value = zeta_offset)
input_yaml(file = gotm_yaml, label = 'A', key = 'constant_value', value = A)
input_yaml(file = gotm_yaml, label = 'g1', key = 'constant_value', value = g1)
input_yaml(file = gotm_yaml, label = 'g2', key = 'constant_value', value = g2)

# Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 

# Now we are going to load in data fron the NetCDF file into R
out <- 'output.nc' # output file

list_vars(ncdf = out, long = T) # View long names of variables in netCDF
list_vars(ncdf = out, long = F) #View short names of variables in netCDF, corresponds to positions of long names

lnams <- list_vars(ncdf = out, long = T) # Store long names
index <- which(lnams == "potential temperature") # Store index of where temp is
snam <- list_vars(ncdf = out, long = F)[index] # Extracts short name
snam

mod_wtemp <- get_vari(ncdf = out, var = snam) # Extract potential temp
head(mod_wtemp) # View first 6 rows

z <- get_vari(out, 'z') # Extract depths corresponding to temp
head(z)

p1 <- plot_wtemp(out, size = 2) # Plot water temperature, adjust 'size' to fill empty space
p1 <- p1 + xlab('') + ylab('Depth (m)') + ggtitle('Feeagh - Modelled Heatmap') # Add axis labels and title
p1 # Prints plot to screen


### But is our model any good?
# We need to compare this with actual observed data
obs <- load_obs(obs.file = obs.file)
head(obs, 14)
long_heatmap(obs)
long_lineplot(obs)

# Sub set to same period as modelled
obs <- obs[(obs$date >= start & obs$date <= stop),]
p2 <- long_lineplot(obs, main = 'Feeagh - Observed')

ggarrange(p1, p2, nrow = 2) #View plots next to each other

# To do a direct comparison we will need to extract the same depths as observed
mod2 <- setmodDepths(mod.val = mod_wtemp, mod.dep = z, obs = obs)
head(mod2, 14)
ylims <- range(mod2[,3], obs[,3])

p3 <- long_lineplot(mod2, main = ' Feeagh - Modelled')
p3 <- p3 + coord_cartesian(ylim = ylims)
p2 <- p2 + coord_cartesian(ylim = ylims)
ggarrange(p3, p2, nrow = 2) #View plots next to each other

g1 <- diag_plots(mod = mod2, obs = obs, colourblind = F)
g1
sum_stat(mod2, obs, depth = T)

# Uh-oh, it looks like our model is not replicating observed conditions very well!
# Can you think of any potential sources of error?
# How can we adjust our model to account of for these sources of error?

# Calibration
db_file <- 'feeagh_calib_resuts.db'
config_file <- 'calib_config.xml'

pars <- get_param(dbFile = db_file, acpyXML = config_file)
summary(pars)
mlt <- melt(pars, id.vars = c('id', 'run','time', 'lnlikelihood', 'correlation', 'rmse'))
mlt$plike <- mlt$lnlikelihood - min(mlt$lnlikelihood)

ggplot(mlt, aes(value, rmse, colour = run))+
  geom_point()+
  facet_wrap(~variable, scales = 'free')+
  coord_cartesian(ylim = c(0,3))+
  theme_bw()
ggplot(mlt, aes(value, plike, colour = run))+
  geom_point()+
  facet_wrap(~variable, scales = 'free')+
  scale_y_log10()+
  coord_cartesian(ylim = c(max(mlt$plike),max(mlt$plike)-10000))+
  theme_bw()

best_par <- pars[which.max(pars$lnlikelihood),]

input_yaml(file = gotm_yaml, label = 'heat', key = 'scale_factor', value = best_par$`surface/fluxes/heat/scale_factor`)
input_yaml(file = gotm_yaml, label = 'swr', key = 'scale_factor', value = 1.2)
input_yaml(file = gotm_yaml, label = 'turb_param', key = 'k_min', value = best_par$`turbulence/turb_param/k_min`)
input_yaml(file = gotm_yaml, label = 'g1', key = 'constant_value', value = best_par$`light_extinction/g1/constant_value`)
input_yaml(file = gotm_yaml, label = 'g2', key = 'constant_value', value = best_par$`light_extinction/g2/constant_value`)

run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat') 
plot_wtemp(out, size = 2)

# Extract temp and format again to match observed data
mod_wtemp <- get_vari(ncdf = out, var = 'temp') # Extract potential temp
z <- get_vari(out, 'z') # Extract depths corresponding to temp
mod3 <- setmodDepths(mod.val = mod_wtemp, mod.dep = z, obs = obs)

g2 <- diag_plots(mod = mod3, obs = obs, colourblind = F)
g2 #Now we have a much better fit
sum_stat(mod3, obs, depth = T)


#Examine an event
obs_sub <- obs[(obs[,1] >= '2015-07-01' & obs[,1] <= '2015-08-01'),]
mod_sub <- mod3[(mod3[,1] >= '2015-07-01' & mod3[,1] <= '2015-08-01'),]
ylims <- range(mod_sub[,3], obs_sub[,3])


p5 <- long_lineplot(obs_sub, main = 'Feeagh - Observed')+ coord_cartesian(ylim = ylims)
p6 <- long_lineplot(mod_sub, main = ' Feeagh - Modelled')+ coord_cartesian(ylim = ylims)
ggarrange(p5,p6)

g3 <- diag_plots(mod = mod_sub, obs = obs_sub, colourblind = F)
g3 #Now we have a much better fit
sum_stat(mod_sub, obs_sub, depth = T)
