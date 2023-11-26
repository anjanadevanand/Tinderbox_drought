import xarray as xr
import numpy as np
import pandas as pd
import os
from statsmodels.formula.api import glm
import statsmodels.formula.api as smf
import statsmodels.api as sm
import itertools
import datetime
from dask.distributed import Client, LocalCluster
from dask.distributed import Client,Scheduler
from dask_jobqueue import SLURMCluster
import time
import glob
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from joblib import Parallel, delayed

# #########################
# LOADED ALL ARRAYS FOR THE SUBSET CASE - probably can't do this with the full data

def create_filepath_oneTime(ds, prefix='filename', root_path="."):
    """
    Generate a filepath when given an xarray dataset
    """
    time_str = ds.time.dt.strftime("%Y-%m-%d").data
    filepath = f'{root_path}/{prefix}_{time_str}.nc'
    return filepath

# define a function to fit the GLM model
def fit_logistReg_3Pred_oneThres_handleNan(y, x1, x2, x3, predictors, thres, formula, x1_new, x2_new, x3_new):
    '''Function to fit a logistic regression model to estimate exceedence probability
    using a thres argument. If the thres argument is nan, the grid is not in drought.
    '''
    n_predictors = len(predictors) + 1
    GLM_params = np.empty(n_predictors)
    GLM_pvalues = np.empty(n_predictors)
    GLM_probability = np.empty(1)
    GLM_aic = np.empty(1)

    if np.isnan(thres):
        GLM_probability[:] = np.nan
        GLM_params[:] = np.nan
        GLM_pvalues[:] = np.nan
        GLM_aic[:] = np.nan
    else:
        y_binary = np.where(y >= thres, 1, 0)
        y_binary = y_binary[~np.isnan(y)]
        if (sum(y_binary) < 4):                      # too few data points for estimation
            GLM_probability[:] = 0
            GLM_params[:] = 0
            GLM_pvalues[:] = np.nan
            GLM_aic[:] = np.nan
        else:                                        # logistic regression fit
            # create a dataframe of reponse and predictors
            x_dict = {predictors[0]:x1[~np.isnan(y)], predictors[1]:x2[~np.isnan(y)], predictors[2]:x3[~np.isnan(y)]}
            x = pd.DataFrame(x_dict)
            x['response'] = y_binary

            x_new_dict = {predictors[0]:x1_new, predictors[1]:x2_new, predictors[2]:x3_new}
            x_new = pd.DataFrame(x_new_dict, index = [0])

            model = glm(formula, x, family=sm.families.Binomial())
            try:
                model_GLM = model.fit()
                GLM_probability[:] = model_GLM.predict(x_new)
                GLM_params[:] = model_GLM.params.values
                GLM_pvalues[:] = model_GLM.pvalues.values
                GLM_aic[:] = model_GLM.aic
            except PerfectSeparationError:          # this error occurs at longer timescales with fewer data points
                GLM_probability[:] = np.nan
                GLM_params[:] = np.nan
                GLM_pvalues[:] = np.nan
                GLM_aic[:] = np.nan
    return GLM_params, GLM_pvalues, GLM_probability, GLM_aic

def identify_start_year(glm_dir, start_lat, end_lat, start_lon, end_lon, start_year_default):
    start_year = start_year_default
    # first check the yearly files in the top level directory
    fnames = glm_dir + 'GLM_results_lat' + str(start_lat) + '*_lon' + str(start_lon) + '*.nc'
    yearly_files = sorted(glob.glob(fnames))
    if (len(yearly_files)>0):
        last_file = yearly_files[-1]
        # split the last file name to get the corresponding year
        start_year = int(((last_file.split('/')[-1]).split('_')[-1]).split('.')[0])
        if start_year == 2020:
            start_year = None

    # check incomplete daily files in the corresponding sub-folder
    sub_dir = 'lat' + str(start_lat) + '_' + str(end_lat) + '_lon' + str(start_lon) + '_' + str(end_lon)
    daily_files = sorted(glob.glob(glm_dir + sub_dir + '/by_day/GLM_results_*.nc'))
    if (len(daily_files)>0):
        last_file = daily_files[-1]
        # split the last file name to get the corresponding date
        start_date = ((last_file.split('/')[-1]).split('_')[-1]).split('.')[0]
        start_year_temp = int(start_date.split('-')[0])
        #if there are daily files that haven't been concated to yearly files yet (comment if condition to identify start year from daily files only)
        if start_year_temp >= start_year:
            # it is better to start from the year before because parallel jobs may not have completed that year completely
            start_year = start_year_temp - 1
        ############### To run only one missing year; also comment the above if loop
        # start_year = start_year_temp
        ###############
        start_mon = int(start_date.split('-')[1])
        start_day = int(start_date.split('-')[2])
        if start_year_temp == 2020 & start_mon == 12 & start_day == 31:
            start_year = None
    return start_year

if __name__ == '__main__':
    
    subset_time = False #True 
    
    # spatial subset
    start_lat = float(os.environ['start_lat'])
    end_lat = float(os.environ['end_lat'])
    
    start_lon = float(os.environ['start_lon'])
    end_lon = float(os.environ['end_lon'])
    
    lat_slice = slice(start_lat, end_lat)
    lon_slice = slice(start_lon, end_lon)
    sub_dir = 'lat' + str(start_lat) + '_' + str(end_lat) + '_lon' + str(start_lon) + '_' + str(end_lon)
    
    # lat_list = [-33, -34, -35, -36]
    # lon_list = [141, 142, 143, 144]
    
    main_dir = '/g/data/w97/ad9701/drought_2017to2020/drought_breakProb/awra/'

    varname = 'sm' #'P'   # the name of the directory
    vname = 'sm_diff'   #'precip'  # the name of the files and variable inside the files
    fname = vname + '_*_*_*.nc'

    iW = int(os.environ['iWeek'])
    print(iW)
    start_yr = int(os.environ['start_yr'])
    end_yr = int(os.environ['end_yr'])
    drght_dir = os.environ['glm_dir'] #specify the full path here 'GLM_results_model3_subset/deficits_basedOn_futureDay/'

    start_yr = identify_start_year(glm_dir = drght_dir + varname + '_week' + str(iW) + '/', start_lat = start_lat, end_lat = end_lat, start_lon = start_lon, end_lon = end_lon, start_year_default = start_yr)
    print('start_yr=')
    print(start_yr)
    
    ########### to run for only a single year - identified using the 'identify_start_year' function 
    #end_yr = start_yr
    ###########
    print('end_yr=')
    print(end_yr)
    if start_yr is not None:
        # select thresholds
        # load the threshold data file & select the drought period of interest
        sm_deficit_files = os.environ['sm_deficit_files'] #'sm_droughts/deficits_basedOn_futureDay/sm_events_[1-2]*.nc'
        events_file = os.environ['sm_events_files'] #'sm_droughts/events_[1-2]*.nc'

        ds_thresh = xr.open_mfdataset(main_dir + sm_deficit_files)
        ds_events = xr.open_mfdataset(main_dir + events_file)   # binary data indicating whether the grid is in drought
        # print(ds_events)

        #### NOTE: training dataset does not include the drought period - using the 1980 to 2016 baseline
        overall_time_slice = slice('1980-01-01', '2016-12-31')
        drght_time_slice = slice(str(start_yr)+'-01-01', str(end_yr)+'-12-31')
        #drght_dir = os.environ['glm_dir'] #specify the full path here 'GLM_results_model3_subset/deficits_basedOn_futureDay/'

        # select the thresholds for the time periods of the drought
        thresName = 'sm_deficit'
        # da_thresh = ds_thresh[thresName].sel(time = drght_time_slice, timescale = iW).persist()
        da_thresh = ds_thresh[thresName].sel(time = drght_time_slice, lat = lat_slice, lon = lon_slice, timescale = iW).persist()

        ############################################
        # GET THE SST PREDICTORS
        ############################################

        # get the sst data
        sst_dir = '/g/data/w97/ad9701/p_prob_analysis/sst_data/'
        pNames = ['soi', 'sami', 'dmi', 'nino34_anom', 'nino4_anom']
        pFiles = ['soi_monthly.nc', 'newsam.1957.2021.nc', 'dmi.had.long.data.nc', 'nino34.long.anom.data.nc', 'nino4.long.anom.data.nc']
        for p in np.arange(len(pNames)):
            ds_temp = xr.open_dataset(sst_dir+pFiles[p])
            if (p>0):
                ds_p[pNames[p]]=ds_temp[pNames[p]]
            else:
                ds_p = ds_temp
            del ds_temp

        # select the predictors to include in the model
        predSel = ['soi', 'dmi', 'sm']
        formula = 'response ~ soi + dmi + sm'
        parameter = ['Intercept']
        parameter.extend(predSel)

        # select the sst predictors corresponding to the dates of the thresholds data
        thresh_time_bymon = np.array(pd.to_datetime(da_thresh.time).to_period('M').to_timestamp().floor('D'))
        da_p1_current = ds_p['soi'].sel(time = thresh_time_bymon).persist()
        da_p2_current = ds_p['dmi'].sel(time = thresh_time_bymon).persist()
        # da_p3_current = ds_events['sm_drought'].sel(time = da_thresh.time, lat = lat_list, lon = lon_list).persist()

        ############################################
        # START A LOCAL CLUSTER
        ############################################
        # from dask.distributed import Client, LocalCluster
        # cluster = LocalCluster()
        # client = Client(cluster)
        # client

        ############################################
        # PERFORM CALCULATIONS FOR THE MAIN SET
        ############################################

        # get data
        data_dir = '/g/data/w97/ad9701/p_prob_analysis/processed_data/awra/' + varname + '_week' + str(iW) + '/'
        print(data_dir)
        ds_var_temp = xr.open_mfdataset(data_dir + fname) #, chunks = {'lat':400, 'lon':400})
        ds = ds_var_temp.sel(time = overall_time_slice)

        # da_var_temp = ds[vname].reindex(lat=ds.lat[::-1]).chunk(chunks = {'lat':40,'lon':40,'time':-1}).rename({'time':'hist_time'})
        # da_var_temp = ds[vname].reindex(lat=ds.lat[::-1]).sel(lat = lat_list, lon = lon_list).chunk(chunks = {'lat':40,'lon':40,'time':-1}).rename({'time':'hist_time'})
        da_var_temp = ds[vname].sel(lat = lat_slice, lon = lon_slice).chunk(chunks = {'time':-1}).rename({'time':'hist_time'})

        # if subset_time:
        #     time_sel=slice('1981-01-01', '2020-05-31')
        #     ds = ds.sel(time = time_sel)
        #     da_var_temp = da_var_temp.sel(hist_time = time_sel)

        # da_drought = ds_events['sm_drought'].rename({'time':'hist_time'}).sel(hist_time = da_var_temp.hist_time).chunk(chunks = {'hist_time':-1})
        da_drought = ds_events['sm_drought'].rename({'time':'hist_time'}).sel(lat = lat_slice, lon = lon_slice, hist_time = da_var_temp.hist_time).chunk(chunks = {'hist_time':-1})
        da_var_drought = da_var_temp.where(da_drought == 1).persist().groupby('hist_time.season')     # the no-drought data points are set to nan

        # select predictors for the same time points as the P-E or P-E-Q data at multi-weekly timescale
        da_time_bymon = np.array(pd.to_datetime(ds.time).to_period('M').to_timestamp().floor('D'))
        ds_p_sel = ds_p.sel(time = da_time_bymon)
        ds_p1_sel_gb = ds_p_sel['soi'].rename({'time':'hist_time'}).persist().groupby('hist_time.season')
        ds_p2_sel_gb = ds_p_sel['dmi'].rename({'time':'hist_time'}).persist().groupby('hist_time.season')
        # ds_p3_sel_temp = ds_events['sm_drought'].sel(time = ds.time, lat = lat_list, lon = lon_list).rename({'time':'hist_time'})
        # ds_p3_sel_gb = ds_p3_sel_temp.chunk(chunks = {'lat':40,'lon':40,'hist_time':-1}).persist().groupby('hist_time.season')

        ############################################
        # Add SM as a predictor
        ############################################

        # intial sm values corresponding to P-E-Q accumulations - did not use this, just read in all the required values from the AWRA full dataset. (One day difference problem)
        # data_dir_sm = main_dir + 'sm_week'+str(iW)+'/' + sub_dir + '/'
        # fname_sm_init = 'sm_init_*_*_*.nc'
        # ds_sm_init = xr.open_mfdataset(data_dir_sm + fname_sm_init, chunks = {'lat':400, 'lon':400})
        # da_sm_init_temp = ds_sm_init['sm'].reindex(lat=ds_sm_init.lat[::-1]).sel(lat = lat_list, lon = lon_list).chunk(chunks = {'lat':40,'lon':40,'time':-1}).rename({'time':'hist_time'})
        # da_sm_init = da_sm_init_temp.persist().groupby('hist_time.season')

        # get all sm values from the awra data
        awra_dir = '/g/data/fj8/BoM/AWRA/DATA/SCHEDULED-V6/processed/values/day/'
        file_names = 'sm_[1-2]*.nc'
        ds_temp = xr.open_mfdataset(awra_dir + file_names, chunks = {'lat':400,'lon':400})
        # converting the datatypes of SM to match P
        lat_new = np.float32(ds_temp['latitude'])
        lon_new = np.float32(ds_temp['longitude'])
        da_sm = ds_temp['sm'].rename({'latitude':'lat','longitude':'lon'}).assign_coords(lat = lat_new, lon = lon_new)

        # values at times corresponding to da_var to estimate parameters of the glm model
        # ds_p3_sel_temp = da_sm.reindex(lat=da_sm.lat[::-1]).sel(time = ds.time).rename({'time':'hist_time'})
        # ds_p3_sel_temp = da_sm.reindex(lat=da_sm.lat[::-1]).sel(time = ds.time, lat = lat_list, lon = lon_list).rename({'time':'hist_time'})
        ds_p3_sel_temp = da_sm.sel(time = ds.time, lat = lat_slice, lon = lon_slice).rename({'time':'hist_time'})
        ds_p3_sel_gb = ds_p3_sel_temp.chunk(chunks = {'hist_time':-1}).transpose('lat','lon','hist_time').persist().groupby('hist_time.season')

        # current values at times corresponding to da_thresh for prediction
        da_p3_current = da_sm.sel(time = da_thresh.time, lat = lat_slice, lon = lon_slice).persist()
        # da_p3_current = da_sm.sel(time = da_thresh.time).persist()

        full_dir_path = drght_dir + '/' + varname + '_week' + str(iW) + '/' + sub_dir + '/by_day/'
        if not os.path.exists(full_dir_path):
            os.makedirs(full_dir_path) 

        start_day = 0
        # check for existing half completed files - this doesn't work properly for parallel jobs
        # glm_files = glob.glob(full_dir_path + 'GLM_results_*' + str(start_yr) + '-*.nc')
        # nfiles = len(glm_files)
        # if nfiles < 367:
        #    if nfiles > 4:
        #        start_day = nfiles - 4
        #    else:
        #        start_day = 0
        dask_gufunc_kwargs = {'output_sizes':{"glm_parameter": len(parameter)}} #, 'allow_rechunk': True} #, 'time':1}}
        input_core_dims=[["hist_time"], ["hist_time"], ["hist_time"], ["hist_time"], ["predictors"], [], [], [], [], []]
        output_core_dims=[["glm_parameter"], ["glm_parameter"], [], []]

        def do_workflow_current_time(i_time, da_thresh=da_thresh, da_var_drought=da_var_drought, ds_p1_sel_gb=ds_p1_sel_gb, 
                                     ds_p2_sel_gb=ds_p2_sel_gb, ds_p3_sel_gb=ds_p3_sel_gb, predSel=predSel, formula=formula, 
                                     da_p1_current=da_p1_current, da_p2_current=da_p2_current, da_p3_current=da_p3_current, 
                                     input_core_dims = input_core_dims, output_core_dims = output_core_dims,
                                     dask_gufunc_kwargs=dask_gufunc_kwargs, parameter=parameter):

            seas = da_thresh['time.season'].values[i_time]
            da_logistReg = xr.apply_ufunc(
                fit_logistReg_3Pred_oneThres_handleNan,        # first the function, this function returns a tuple (GLM params, GLM pvalues, GLM modelled probabilities)
                da_var_drought[seas],                          # function arg
                ds_p1_sel_gb[seas].values,
                ds_p2_sel_gb[seas].values,
                ds_p3_sel_gb[seas].values,
                predSel,                                     #      "
                da_thresh.isel(time = i_time),                                  #      "
                formula,                                     #      "
                [da_p1_current.isel(time = i_time).values],                    #      "
                [da_p2_current.isel(time = i_time).values], #      "
                da_p3_current.isel(time = i_time),
                input_core_dims=input_core_dims, #["sample_time"], ["sample_time"]],   # list with one entry per arg, these are the dimensions not to be broadcast
                output_core_dims=output_core_dims,                                # dimensions of the output
                vectorize=True,                                                                                 # broadcast over non-core dimensions of the input object?
                dask="parallelized",                                                                                                               # enable dask?
                dask_gufunc_kwargs=dask_gufunc_kwargs,                     
                output_dtypes=[float, float, float, float]
            )

            # assign co-ordinates add metadata
            new_coords_dict = {'glm_parameter':parameter} #, 'current_time':[da_thresh['current_time'][i_time]]}    
            ds_all = da_logistReg[2].rename('glm_probability').to_dataset()
            ds_all['glm_params'] = da_logistReg[0].rename('glm_params').assign_coords(new_coords_dict)
            ds_all['glm_pvalues'] = da_logistReg[1].rename('glm_pvalues').assign_coords(new_coords_dict)
            ds_all['glm_aic'] = da_logistReg[3].rename('glm_aic')
            ds_all[predSel[0]] = da_p1_current.isel(time = i_time)
            ds_all[predSel[1]] = da_p2_current.isel(time = i_time)

            out_file = create_filepath_oneTime(ds_all, prefix = 'GLM_results_' + '_'.join(predSel), root_path = full_dir_path)
            #print(out_file)
            ds_all.to_netcdf(out_file)
            return i_time

        Parallel(n_jobs=28)(delayed(do_workflow_current_time)(i_time) for i_time in range(start_day, len(da_thresh.time)))
