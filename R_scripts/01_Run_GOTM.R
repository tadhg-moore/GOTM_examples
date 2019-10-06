#First set your working directory to the lake folder
setwd("...\\feeagh")

library(GOTMr)
library(gotmtools)

#View 'gotm.yaml' file  #Hold 'Ctrl' and left-click 'gotm.yaml'
gotm_yaml <- 'gotm.yaml'

#Check GOTM version
gotm_version()
# Should get the same input as below...
# ------------------------------------------------------------------------
#   GOTM:    5.4.0 (unknown branch)
# YAML:    0.1.0 (unknown branch)
# flexout: 0.1.0 (unknown branch)
# STIM:     (unknown branch)
# FABM:    unknown (unknown branch)
# NetCDF: 3.6.1-beta1 of Mar  7 2018 17:17:53 $
#   ------------------------------------------------------------------------
#   Compiler: Intel 19.0.0.20181018


# If you want to work with a template version of GOTM you can find one
template_files()
# ".../GOTMr/extdata" Go here on your computer and you will find template files if you want to set up your own simulation.


#Before running GOTM we will step through the 'gotm.yaml' file...
#################################################################
#################################################################



#Run GOTM
run_gotm(yaml_file = gotm_yaml)
# system('run_gotm.bat')  #Use this command if the GOTMr package
# isn't installed

##Check output in the console to make sure that GOTM has ran successfully

# ...
# ------------------------------------------------------------------------
#   GOTM finished on 2019/09/27 at 14:23:12
# ------------------------------------------------------------------------
#   CPU time:                       2.152814      seconds
# Simulated time/CPU time:        13404595.6534740     
# ------------------------------------------------------------------------
#   GOTM:    5.4.0 (unknown branch)
# YAML:    0.1.0 (unknown branch)
# flexout: 0.1.0 (unknown branch)
# STIM:     (unknown branch)
# FABM:    unknown (unknown branch)
# NetCDF: 3.6.1-beta1 of Mar  7 2018 17:17:53 $
#   ------------------------------------------------------------------------
#   Compiler: Intel 19.0.0.20181018

#Run run_gotm
system2("run_gotm.bat")


##If you get the above output, SUCCESS!
##You can now run GOTM in R

##If you want to save the text from the GOTM run
gotm_out <- run_gotm(yaml_file = gotm_yaml, verbose = T)
# gotm_out <- system2("run_gotm.bat", stdout = T)
gotm_out
write.table(gotm_out, file='got_out.txt', row.names = F, quote = F, col.names = F)

#################################################################
############# View GOTM Output ##################################

# There is a program 'PyNcView' which has been developed as a graphical
# viewer for NetCDF files #https://sourceforge.net/projects/pyncview/

# In your directory there is now an 'output.nc' file. Open this with PyNcView
# to view your NetCDF file.



#################################################################
############# Editing parameters ################################

# Before we start editing the the 'gotm.yaml' we are going to make a
# master copy
file.copy(from = 'gotm.yaml', to = 'gotm_master.yaml', overwrite = FALSE)
# [1] TRUE  

# 1. Go into 'gotm.yaml', increase the surface/meteo/swr/scale_factor (L110)
# from 1 to 1.5
# 2. SAVE the file, then 
# 3. Run GOTM
# 4. view NetCDF file and save plot of 'potential temperature'

# # Repeat steps 1-4 by increasing scale factors  for airt, u10, v10, 
# surface/fluxes/heat/scale_factor (L56) and save potential temperature 
# in the 'plots' folder. Hypothesize what you expect will happen with
# each change.
# 
# If you have time also try varying the light extinction parameters
# g1[0-2] and g2[0.5-3] and see what effect that has.

run_gotm(yaml_file = gotm_yaml)
# system2("run_gotm.bat")



#################################################################
############# Editing parameters within R ####################


start <- '2010-04-01 00:00:00'
stop <- '2010-10-01 00:00:00'


input_yaml(file = gotm_yaml, label = 'time', key = 'start', value = start)
input_yaml(file = gotm_yaml, label = 'time', key = 'stop', value = stop)
input_yaml(file = gotm_yaml, label = 'location', key = 'name', value = 'feeagh')
input_yaml(file = gotm_yaml, label = 'roughness', key = 'charnock', value = 'false')

# Create initial temperature file for simulation
init_prof(obs_file = 'wtemp_dly_2010_2015.dat', date = start, tprof_file = 'init_tprof.dat')

run_gotm(yaml_file = gotm_yaml)
# system2("run_gotm.bat")

#View netCDF file to see the model output. 
