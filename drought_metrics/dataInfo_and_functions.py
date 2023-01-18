mswep_info = dict(full_path = ['/g/data/fj4/SatellitePrecip/MSWEP_V280/Past/Monthly/', '/g/data/fj4/SatellitePrecip/MSWEP_V280/NRT/Monthly/'], 
                  file_name = '*.nc',
                  varname = 'precipitation',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = False,
                  land_mask = False,
                   mon_file = 'precipitation_mswep_monthly_1979_2021.nc')
chirps_info = dict(full_path='/g/data/w97/ad9701/CHIRPS-2.0/global_monthly/netcdf/',
                    file_name = 'chirps-v2.0.monthly.nc',
                    varname = 'precip',
                    lat_slice = slice(-44, -10),
                    lon_slice = slice(112, 154),
                    units = 'mm month-1',
                    rename_latlon = True,
                    land_mask = True,
                     mon_file = 'chirps_monthly_1981_2021.nc')
agcd_info = dict(full_path = '/g/data/zv2/agcd/v1/precip/total/r005/01month/',
                file_name = 'agcd_v1_precip_total_r005_monthly_*.nc',
                varname = 'precip',
                lat_slice = slice(-44, -10),
                lon_slice = slice(112, 154),
                units = 'mm month-1',
                rename_latlon = False,
                land_mask = False,
                mon_file = 'agcd_monthly_1900_2020.nc')
# Using PET to calculate SPEI, so variable is set to PET for now
gleam_info = dict(full_path = '/g/data/ua8/GLEAM_v3-5/v3-5a/monthly/', 
                  file_name = 'Ep_1980-2020_GLEAM_v3.5a_MO.nc',
                  varname = 'Ep',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = False,
                  mon_file = 'PminusPET_gleam_monthly_1980_2020.nc')
awra_info = dict(full_path = '/g/data/fj8/BoM/AWRA/DATA/SCHEDULED-V6/', 
                  file_name = 'e0_*.nc',
                  varname = 'e0',
                  lat_slice = slice(-10, -44),
                  lon_slice = slice(112, 154),
                  units = 'mm month-1',
                  rename_latlon = True,
                  mon_file = 'PminusPET_awra_monthly_1911_2020.nc')
esacci_info = dict(full_path = ['/g/data/w97/ad9701/ESACCI/3_unzip/' + str(i) + '/' for i in range(1979, 2021)],
                   file_name = 'ESACCI-SOILMOISTURE*.nc',
                   varname = 'sm',
                   lat_slice = slice(-10, -44),
                   lon_slice = slice(112, 154),
                   units = 'm3 m-3',
                   rename_latlon = False,
                   mon_file = 'esacci_monthly_1979_2020.nc')

# all variables of interest
# ro: runoff (mm)
# sm: soil moisture (1m depth)
# e: Evaporation

era5land_info = dict(full_path = '/g/data/w97/ad9701/p_prob_analysis/era5-land_daily/e/',
                   file_name = 'e_era5-land_oper_sfc_*.nc',
                   varname = 'e',
                   lat_slice = slice(-10, -44),
                   lon_slice = slice(112, 154),
                   units = 'mm',
                   rename_latlon = True,
                   mon_file = '')


alldata_dict = dict(mswep = mswep_info, chirps = chirps_info, agcd = agcd_info, gleam = gleam_info, awra = awra_info, esacci = esacci_info, era5land = era5land_info)

# GLEAM variables
# E: Actual evaporation
# Eb: Bare soil evaporation
# Ei: Interception loss
# Ep: Potential evaporation
# Es: snow sublimation
# Et: transpiration
# Ew: Open-water evaporation
# S: Evaporative stress
# SMroot: rootzone soil mositure
# SMsurf: surface soil moisture


import xarray as xr
import numpy as np
import sys
import glob

def save_monthly_data_byyear(alldata_dict, data_name, out_dir, calc_from_daily = False, calc_fun = 'sum'):
    '''
    Used in cases where there are too may files that need to be processed directory where each directory is a year.
    '''
    if type(data_name) == list:
        data_list = data_name
    elif type(data_name) == str:
        data_list = [data_name]
    else:
        sys.exit("data_name should be a list or a string")
        
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
        
    for d in data_list:
        for path in alldata_dict[d]['full_path']:
            print('working on directory ' + str(path))
            da_AU = get_da_simple(path,
                                  alldata_dict[d]['file_name'],
                                  alldata_dict[d]['varname'],
                                  alldata_dict[d]['rename_latlon'],
                                  alldata_dict[d]['lat_slice'],
                                  alldata_dict[d]['lon_slice'])                     
            if calc_from_daily:
                if calc_fun == 'sum':
                    da_AU_mon = da_AU.resample(time="M").sum()
                elif calc_fun == 'mean':
                    da_AU_mon = da_AU.resample(time="M").mean()
                out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU_mon['time.year'].min().values) + '.nc'
                da_AU_mon.to_netcdf(out_dir + out_file)
            else:
                out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU['time.year'].min().values) + '.nc'
                da_AU.to_netcdf(out_dir + out_file)
    return None

def get_da_simple(full_path, file_name, varname, rename_latlon, lat_slice, lon_slice):
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
    star_loc = file_name.find('*')
    if star_loc == -1:
        ds = xr.open_dataset(full_path + file_name)
    else:
        ds = xr.open_mfdataset(full_path + file_name)

    if rename_latlon:
        da_AU = ds[varname].rename(latlon_rename).sel(lat = lat_slice, lon = lon_slice)
    else:
        da_AU = ds[varname].sel(lat = lat_slice, lon = lon_slice)
    return da_AU

def save_monthly_data(alldata_dict, data_name, out_dir, calc_from_daily = False, calc_fun = 'sum'):
    if type(data_name) == list:
        data_list = data_name
    elif type(data_name) == str:
        data_list = [data_name]
    else:
        sys.exit("data_name should be a list or a string")
        
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
        
    for d in data_list:
        da_AU = get_da(alldata_dict, d)
    if calc_from_daily:
        if calc_fun == 'sum':
            da_AU_mon = da_AU.resample(time="M").sum()
        elif calc_fun == 'mean':
            da_AU_mon = da_AU.resample(time="M").mean()
        print('Calculated monthly data')
        out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU_mon['time.year'].min().values) + '_' + str(da_AU_mon['time.year'].max().values) + '.nc'
        da_AU_mon.to_netcdf(out_dir + out_file)
    else:
        out_file = alldata_dict[d]['varname'] + '_' + d + '_monthly_' + str(da_AU['time.year'].min().values) + '_' + str(da_AU['time.year'].max().values) + '.nc'
        da_AU.to_netcdf(out_dir + out_file)
    return None

def get_da(alldata_dict, d):
    latlon_rename = {'latitude': 'lat', 'longitude': 'lon'}
    fnames = None
    if type(alldata_dict[d]['full_path']) == list:    # data has to be read from multiple locations
        fnames = []
        for p in alldata_dict[d]['full_path']:
            fnames.extend(glob.glob(p + alldata_dict[d]['file_name']))
        # print(fnames)
                
    if fnames is None:
        star_loc = alldata_dict[d]['file_name'].find('*')
        if star_loc == -1:
            ds = xr.open_dataset(alldata_dict[d]['full_path'] + alldata_dict[d]['file_name'])
        else:
            ds = xr.open_mfdataset(alldata_dict[d]['full_path'] + alldata_dict[d]['file_name'])
            # print('completed reading daily files')
    else:
        ds = xr.open_mfdataset(fnames)
    if alldata_dict[d]['rename_latlon']:
        da_AU = ds[alldata_dict[d]['varname']].rename(latlon_rename).sel(lat = alldata_dict[d]['lat_slice'], lon = alldata_dict[d]['lon_slice'])
    else:
        da_AU = ds[alldata_dict[d]['varname']].sel(lat = alldata_dict[d]['lat_slice'], lon = alldata_dict[d]['lon_slice'])
    # print('returning daily data to main func')
    return da_AU

# regridding to common resolution for averaging
import xesmf as xe

def regrid_all_from_list(da_list, lat = np.arange(-10.125, -44.125, -0.25), lon = np.arange(112.125, 154.125, 0.25)):
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lat),
            "lon": (["lon"], lon),
        }
    )
    shape_tuple = (len(lat), len(lon))    
    
    da_list_regrid = []  
    for da in da_list:
        da = da.chunk({'lat':-1, 'lon':-1})
        if (da.values.shape[1:] == shape_tuple) | (da.values.shape == shape_tuple):     # assumed that same shape means data is already at the desired resolution
            da_list_regrid.append(da)
        else:
            regridder = xe.Regridder(da, ds_out, 'bilinear')
            da_reg = regridder(da)
            da_list_regrid.append(da_reg)        
    return da_list_regrid


from cartopy.util import add_cyclic_point
import matplotlib
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import BoundaryNorm
import geopandas as gpd

# packages create shapfiles from contours in the data
from shapely import geometry
from shapely.geometry import Point, Polygon
from matplotlib import cm
import fiona
import os,json
from descartes.patch import PolygonPatch
from fiona.crs import from_epsg


def draw_pcolormesh_Npanels(ds_list, cmap='viridis', levels=None, ncols = 5, nrows = 2, subplot_title=None, main_title=None, out_dir=None, out_figname=None, fig_format= 'png', fig_dpi = 300,
                            add_cbar = True, cbar_extend='both', cbar_label = '', panel_width=3.66, panel_height=4, xticks = np.arange(115,155,10), yticks = np.arange(-40,-10,5), shapefile = None, shapefile_color='black', sh_lwd=2, ds_hatch_list=None, hatches=["..", ".."], hatch_zorder=1, longer_cbar=False, thicker_cbar=False, cbar_ticks=None, cbar_ticklabels=None,
                           projection=ccrs.PlateCarree(), yName='lat', xName='lon', xlim=None, ylim=None):
    '''
    Generic function to create spatial plot figure containing any number of panels, specified using ncols and nrows (defaut set to 5 cols x 2 rows)
    The size of the figure is set using number of panels and individual panel_width & panel height
    Uses matplotlib.pyplot.pcolormesh & adds a common horizontal colorbar at the bottom if add_cbar is True
    '''
    
    if len(ds_list) < (ncols*nrows):
        unwanted_axes = list(range(len(ds_list), ncols*nrows))
    else:
        unwanted_axes = None
        
    ds = ds_list[0]
    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw={'projection': projection},
                            figsize=(panel_width*ncols,panel_height*nrows)) #width, height
    if nrows*ncols > 1:
        axs = axs.flatten()
    else:
        axs = [axs]


    
    if levels is None:
        levels = np.linspace(ds.min(), ds.max(), 5)
    
    norm = None
    if type(cmap) == str:
        cmap_func = plt.colormaps[cmap]
        ncolors = cmap_func.N
        norm = BoundaryNorm(levels, ncolors=ncolors, extend=cbar_extend)
    elif type(cmap) == list:  # a list of colors
        cmap_func = cmap
    else:
        cmap_func = cmap
        ncolors = cmap_func.N
        norm = BoundaryNorm(levels, ncolors=ncolors, extend=cbar_extend)

    for i in np.arange(len(ds_list)):
        #print(i)
        # FAR TOO SMOOTH AND NOT REPRESENTATIVE FOR SPARSE DATA
        # cs=axs[i].contourf(ds_list[i]['lon'],ds_list[i]['lat'],ds_list[i],levels,
        #                       transform = ccrs.PlateCarree(),
        #                       cmap=cmap, extend='both')   #cmap options: coolwarm,
        cs=axs[i].pcolormesh(ds_list[i][xName],ds_list[i][yName],ds_list[i], #,levels=levels,
                      transform = ccrs.PlateCarree(),
                      cmap=cmap, norm=norm, rasterized=True) #, extend='both')   #cmap options: coolwarm,
        if ds_hatch_list is not None:
            axs[i].contourf(ds_hatch_list[i][xName],ds_hatch_list[i][yName],ds_hatch_list[i],levels=1, hatches=hatches,zorder=hatch_zorder, colors='none')   #cmap options: coolwarm,

        # Draw the coastines for each subplot
        axs[i].coastlines()
        axs[i].add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
        if subplot_title is not None:
            axs[i].set_title(subplot_title[i], pad = 2)
        
        if xlim is not None:
            axs[i].set_xlim(xlim)
        
        if ylim is not None:
            axs[i].set_ylim(ylim)      
        
        # COMMENTED FOR TESTING
        #if projection == ccrs.PlateCarree():
        # gl = axs[i].gridlines(crs=projection, draw_labels=False,
        #           linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

# NEW WAY: remove later if not required
#         axs[i].gridlines(xlocs=xticks, ylocs=yticks, draw_labels=False, linewidth=0.5, color='k', alpha=0.5, linestyle='--')
#         axs[i].tick_params(axis='both',labelsize=12,direction='out',right=False,top=False)

        #axs[i].set_extent([-40, -32, 142, 150], crs=ccrs.PlateCarree())
    
        gl = axs[i].gridlines(xlocs=xticks, ylocs=yticks, draw_labels=False, dms=False, x_inline=False, y_inline=False, linewidth=0.5, alpha=0.5, linestyle='--')
        # gl.xlabels_top = False
        # gl.ylabels_right = False

        
        # COMMENTED FOR TESTING
        # gl.xlines = False
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        # gl.xlabel_style = {'rotation': 90}
        # gl.xlabel_style = {'size': 15, 'color': 'gray'}
        # gl.ylabel_style = {'orientation': 'horizontal'}

        # if projection == ccrs.PlateCarree():
        if nrows == 1:
            gl.bottom_labels = True
        else:
            if i in range(ncols*(nrows-1), nrows*ncols):
                # Longitude labels
                gl.bottom_labels = True
                # axs[i].set_xticks(xticks, crs=ccrs.PlateCarree())
                # lon_formatter = cticker.LongitudeFormatter()
                # axs[i].xaxis.set_major_formatter(lon_formatter)

        if i in list(range(0, ncols*nrows, ncols)):
            # Latitude labels
            gl.left_labels = True
            # axs[i].set_yticks(yticks, crs=ccrs.PlateCarree())
            # lat_formatter = cticker.LatitudeFormatter()
            # axs[i].yaxis.set_major_formatter(lat_formatter)

        # axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
        #           linewidth=2, color='gray', alpha=0.5) #, linestyle='--')
        
        #ax.set_extent((30, 55, 20, 45), crs=ccrs.PlateCarree())

        # add a shapefile if specified
        if shapefile is not None:
            shdf = gpd.read_file(shapefile)
            #shdf = shdf.to_crs("EPSG:4326")
            #projection = ccrs.PlateCarree() #central_longitude=0) 
            if projection == ccrs.PlateCarree():
                axs[i].add_geometries(shdf.geometry,
                          projection,
                          facecolor='none',
                          edgecolor=shapefile_color, linewidth=sh_lwd, zorder=5)
            else:
                axs[i].add_geometries(shdf.to_crs(projection).geometry,
                      projection,
                      facecolor='none',
                      edgecolor=shapefile_color, linewidth=sh_lwd, zorder=5)
                
                    
        
            # Delete the unwanted axes
    if unwanted_axes is not None:
        for i in unwanted_axes:
            fig.delaxes(axs[i])

    # # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        if (ncols < 3) | (longer_cbar): #add a colorbar covering a larger proportion of the figure
            if thicker_cbar:
                cbar_ax = fig.add_axes([0.1, 0.15, 0.8, 0.05])
            else:
                cbar_ax = fig.add_axes([0.1, 0.15, 0.8, 0.03])
        else: #shorter colorbar
            if thicker_cbar:
                cbar_ax = fig.add_axes([0.1, 0.15, 0.8, 0.05])
            else:
                cbar_ax = fig.add_axes([0.3, 0.15, 0.4, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)
        if cbar_ticks is not None:
            cbar.set_ticks(cbar_ticks)
            if cbar_ticklabels is not None:
                cbar.set_ticklabels(cbar_ticklabels)

    if main_title is not None:
        plt.suptitle(main_title)
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname + '.' + fig_format, format = fig_format, dpi = fig_dpi, rasterized = True, bbox_inches='tight')
            #plt.savefig(out_dir + out_figname + '.svg', format = 'svg')
        else:
            plt.savefig(out_dir + 'figure.' + fig_format, format = fig_format, dpi = fig_dpi, rasterized = True)
            #plt.savefig(out_dir + 'figure.svg', format = 'svg')
            # NOTE: The svg plots take a lot of time. submit as jobs on Gadi if required
    return fig

def draw_spatial_plot_3panels(ds_list, cmap, levels, subplot_title, main_title, out_dir, out_figname, add_cbar = True, cbar_label = ''):
    ds = ds_list[0]
    fig, axs = plt.subplots(nrows=1,ncols=3,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(24,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    for i in np.arange(len(ds_list)):
        cs=axs[i].contourf(ds_list[i]['lon'],ds_list[i]['lat'],ds_list[i],levels,
                              transform = ccrs.PlateCarree(),
                              cmap=cmap) #,extend='both')   #cmap options: coolwarm,

        # Draw the coastines for each subplot
        axs[i].coastlines()
        axs[i].add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')

        axs[i].set_title(subplot_title[i])

        # Longitude labels
        axs[i].set_xticks(xticks, crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        axs[i].xaxis.set_major_formatter(lon_formatter)
        axs[i].set_xlim(xlim)

        # Latitude labels
        axs[i].set_yticks(yticks, crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        axs[i].yaxis.set_major_formatter(lat_formatter)
        axs[i].set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.85, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.3, 0.15, 0.4, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    plt.suptitle(main_title)
    plt.savefig(out_dir + out_figname)

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

def draw_spatial_plot(ds, cmap, levels, main_title=None, out_dir=None, out_figname=None, add_cbar = True, cbar_label = ''):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat
    
    # xticks = np.arange(113,156,1)  #lon
    # yticks = np.arange(-40,-10,1)   #lat

    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,

    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    if main_title is not None:
        axs.set_title(main_title)

    gl = axs.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right=True
    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
 
########################################################
#     # Longitude labels
#     axs.set_xticks(xticks, crs=ccrs.PlateCarree())
#     lon_formatter = cticker.LongitudeFormatter()
#     axs.xaxis.set_major_formatter(lon_formatter)
#     axs.set_xlim(xlim)

#     # Latitude labels
#     axs.set_yticks(yticks, crs=ccrs.PlateCarree())
#     lat_formatter = cticker.LatitudeFormatter()
#     axs.yaxis.set_major_formatter(lat_formatter)
#     axs.set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)
        
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')


def draw_spatial_plot_addcontours(ds, cmap, levels, main_title, out_dir, out_figname, 
                                 ds_contour_list, contour_level, contour_labels, contour_colors, add_cbar = True, cbar_label = ''):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,
    
    contour_elements_list = []
    for i in range(len(ds_contour_list)):
        cs_contour=axs.contour(ds_contour_list[i]['lon'],ds_contour_list[i]['lat'],ds_contour_list[i],contour_level,
                    transform = ccrs.PlateCarree(),
                    colors=contour_colors[i]) #,extend='both')   #cmap options: coolwarm,
        h1,_ = cs_contour.legend_elements()
        contour_elements_list.append(h1[0])
    
    axs.legend(contour_elements_list, contour_labels, loc='lower left')
    # plt.clabel(cs_contour, inline=1, fontsize=10)
    # cs_contour.collections[0].set_label(contour_label)
    # plt.legend(loc='lower left')
    
    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    axs.set_title(main_title)

    # Longitude labels
    axs.set_xticks(xticks, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs.xaxis.set_major_formatter(lon_formatter)
    axs.set_xlim(xlim)

    # Latitude labels
    axs.set_yticks(yticks, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs.yaxis.set_major_formatter(lat_formatter)
    axs.set_ylim(ylim)

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.3, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    plt.savefig(out_dir + out_figname)
    
def draw_spatial_plot_12panels(ds_list, cmap, levels, subplot_title=None, main_title=None, out_dir=None, out_figname=None, 
                               add_cbar = True, cbar_label = ''): #, unwanted_axes=None):
    npanels = len(ds_list)
    if npanels < 12:
        unwanted_axes = list(range(npanels, 12))
    else:
        unwanted_axes = None
        
    ds = ds_list[0]
    fig, axs = plt.subplots(nrows=2,ncols=6,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(22,8)) #width, height
    axs = axs.flatten()

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat

    for i in np.arange(len(ds_list)):
        #print(i)
        cs=axs[i].contourf(ds_list[i]['lon'],ds_list[i]['lat'],ds_list[i],levels,
                              transform = ccrs.PlateCarree(),
                              cmap=cmap, extend='both')   #cmap options: coolwarm,

        # Draw the coastines for each subplot
        axs[i].coastlines()
        axs[i].add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
        if subplot_title is not None:
            axs[i].set_title(subplot_title[i], pad = 2)
        
        axs[i].set_ylim(ylim)
        axs[i].set_xlim(xlim)
        
        if i in [6, 7, 8, 9, 10, 11]:
            # Longitude labels
            axs[i].set_xticks(xticks, crs=ccrs.PlateCarree())
            lon_formatter = cticker.LongitudeFormatter()
            axs[i].xaxis.set_major_formatter(lon_formatter)

        if i in [0, 6]:
            # Latitude labels
            axs[i].set_yticks(yticks, crs=ccrs.PlateCarree())
            lat_formatter = cticker.LatitudeFormatter()
            axs[i].yaxis.set_major_formatter(lat_formatter)

    # Delete the unwanted axes
    if unwanted_axes is not None:
        for i in unwanted_axes:
            fig.delaxes(axs[i])

    # # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.3, 0.15, 0.4, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)

    if main_title is not None:
        plt.suptitle(main_title)
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')
   
import geopandas as gpd

def draw_spatial_plot_addsh(ds, cmap, levels, main_title=None, out_dir=None, out_figname=None, add_cbar = True, cbar_label = '', shapefile=None):
    fig, axs = plt.subplots(nrows=1,ncols=1,
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(12,9)) #width, height

    xlim = [ds['lon'].values.min(), ds['lon'].values.max()]
    ylim = [ds['lat'].values.min(), ds['lat'].values.max()]

    xticks = np.arange(115,155,10)  #lon
    yticks = np.arange(-40,-10,5)   #lat
    
    cs=axs.contourf(ds['lon'],ds['lat'],ds,levels,
                          transform = ccrs.PlateCarree(),
                          cmap=cmap) #,extend='both')   #cmap options: coolwarm,

    # Draw the coastines for each subplot
    axs.coastlines()
    axs.add_feature(cfeature.OCEAN, zorder=2, edgecolor='k', facecolor='w')
    if main_title is not None:
        axs.set_title(main_title)

    # xticks = np.arange(113,156,1)  #lon
    # yticks = np.arange(-40,-10,1)   #lat
    # gl = axs.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    # gl.xlabels_top = False
    # gl.ylabels_left = False
    # gl.ylabels_right=True
    # gl.xlines = True
    # gl.ylines = True
    # gl.xlocator = mticker.FixedLocator(xticks)
    # gl.ylocator = mticker.FixedLocator(yticks)
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
 
########################################################
    # Longitude labels
    axs.set_xticks(xticks, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs.xaxis.set_major_formatter(lon_formatter)
    axs.set_xlim(xlim)

    # Latitude labels
    axs.set_yticks(yticks, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs.yaxis.set_major_formatter(lat_formatter)
    axs.set_ylim(ylim)
    
    if shapefile is not None:
        shdf = gpd.read_file(shapefile)
        #shdf = shdf.to_crs("EPSG:4326")
        projection = ccrs.PlateCarree() #central_longitude=0) 
        axs.add_geometries(shdf.geometry,
                  projection,
                  facecolor='none',
                  edgecolor='red', zorder=3)
        
        

    # Delete the unwanted axes
    # for i in [5]:
    #     fig.delaxes(axs[i])

    # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.95, left=0.05, right=0.95,
                        wspace=0.1, hspace=0.08)

    if add_cbar:
        # Add a colorbar axis at the bottom of the graph
        cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.03])

        # Draw the colorbar
        cbar=fig.colorbar(cs, cax=cbar_ax,orientation='horizontal', label=cbar_label)
        
    if out_dir is not None:
        if out_figname is not None:
            plt.savefig(out_dir + out_figname)
        else:
            plt.savefig(out_dir + 'figure.png')