# Tinderbox_drought

Repository to share code used for analyses of the 2017-19 Tinderbox Drought in southeast Australia

- impact based drought indicator/ML_Database_All_AWRA_MOf_and_3MPrecip - used to train machine learning model

- drought_metrics/sp*_calc_gridded* - used to calculate SPI-3 and SPEI-3

- drought_metrics/plot_drought_focusReg_meanSPI.ipynb & drought_metrics/plot_drought_focusReg_meanSPEI.ipynb - used to create spatial plot of mean SPI-3 and SPEI-3

- water_cycle_analysis/*.ipynb - used to plot monthly/seasonal anomalies in water cycle variables & save the anomaly numbers in csv files

- veg_analysis/*.R - code to analyse satelliet based vegetation datasets

- land_atmosphere_feedbacks/soil_moisture_drought_impact_on_summer_temperature.py - used to plot difference map of soil drought impact on summer temperature

- synoptic_analysis/Fig8a_code.ncl - used to plot the differences in anomalous proportions of heavy rainfall

- drought_break_probability/
  - identify_gridlevel_drought_events.ipynb & drought_events_calc_smDeficit.ipynb: identify soil moisture droughts and calulate deficits
  - fit_logiReg_gridded_varyThresh_model4_sm_drought_parallel_tb.py : fit logistic regression models to estimate probability of exceedance of deficits
  - plot_drought_break_probability.ipynb : plot modelled drought breaking probability results

- climate_change_analysis/Fig_4A_and_9_Devanand_etal.ipynb: code used to perform bootstrapping of observations, and analyses of climate change signal from CMIP6 data 
