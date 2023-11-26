# Tinderbox_drought

Repository to share code used for analyses of the 2017-19 Tinderbox Drought in southeast Australia

Please add the a line in this readme file about the code you commit along the lines of 'file_name.R - used to calculate SPI3'

1. impact based drought indicator/ML_Database_All_AWRA_MOf_and_3MPrecip - used to train machine learning model

2. drought_metrics/sp*_calc_gridded* - used to calculate SPI-3 and SPEI-3

3. drought_metrics/plot_drought_focusReg_meanSPI.ipynb & drought_metrics/plot_drought_focusReg_meanSPEI.ipynb - used to create spatial plot of mean SPI-3 and SPEI-3

4. water_cycle_analysis/precipitation_and_ET_monthly_plots.ipynb - used to plot monthly anomalies in precipitation and ET (first panel of Figure 2) & save the anomaly numbers in csv files

5. land_atmosphere_feedbacks/soil_moisture_drought_impact_on_summer_temperature.py - used to plot difference map of soil drought impact on summer temperature

6. synoptic_analysis/Fig8a_code.ncl - used to plot the differences in anomalous proportions of heavy rainfall shown in Fig 8A. 

7. drought_break_probability/
  - identify_gridlevel_drought_events.ipynb & drought_events_calc_smDeficit.ipynb: identify soil moisture droughts and calulate deficits
  - fit_logiReg_gridded_varyThresh_model4_sm_drought_parallel_tb.py : fit logistic regression models to estimate probability of exceedance of deficits
  - plot_drought_break_probability.ipynb : plot modelled drought breaking probability results 
