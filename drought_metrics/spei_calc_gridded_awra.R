# read data
#===============================================================================

setwd("/Users/z3387998/OneDrive - UNSW/recent_drought/SPEI")

library(ncdf4) # package for netcdf manipulation
library(ncdf4.helpers)
library(SPEI)

# outfile = "SPEI3_gleam_monthly_1980_2020.nc"
# infile = "PminusPET_gleam_monthly_1980_2020.nc"

infile_name = "PminusPET_awra_monthly_1911_2020.nc"

infile <- nc_open(infile_name)

varname = "PminusPET"
# variable = ncvar_get(infile, varname)
time = ncvar_get(infile, "time")
time_atts = ncatt_get(infile, "time")
lat = ncvar_get(infile, "lat")
lon = ncvar_get(infile, "lon")
Time_in_Date_format <- nc.get.time.series(f = infile,
                                          time.dim.name = "time")
nc_close(infile)

# calculate SPEI3
#===============================================================================

year = 1900 + as.POSIXlt(Time_in_Date_format)$year
mon = 1 + as.POSIXlt(Time_in_Date_format)$mon

# # dividing the data into four sets spatially for calculation
# half_lat = round(length(lat)/2)
# # lat_divisions = list(lat[1:half_lat], lat[half_lat+1, ])
# lat_start = list(1, half_lat+1)
# lat_count = list(half_lat, -1)
# 
# half_lon = round(length(lon)/2)
# # lon_divisions = list(lon[1:half_lon], lat[half_lon+1, ])
# lon_start = list(1, half_lon+1)
# lon_count = list(half_lon, -1)
# 
# count = 0
# for (iLatDiv in 1:length(lat_start)) {
#   for (iLonDiv in 1:length(lon_start)) {
#     count = count + 1
#     outfile = paste0("SPEI3_awra_monthly_1981_2020_set", count, ".nc")
#     
#     # get the subsetted variable
#     infile <- nc_open(infile_name)
#     variable = ncvar_get(nc = infile, varid = varname, start = c(lon_start[[iLonDiv]], lat_start[[iLatDiv]], 1), count = c(lon_count[[iLonDiv]], lat_count[[iLatDiv]], -1))
#     lon_sel = ncvar_get(nc = infile, varid = "lon", start = c(lon_start[[iLonDiv]]), count = c(lon_count[[iLonDiv]]))
#     lat_sel = ncvar_get(nc = infile, varid = "lat", start = c(lat_start[[iLatDiv]]), count = c(lat_count[[iLatDiv]]))
#     nc_close(infile)
#     
#     spei3_gridded = array(NA, dim(variable))
#     for (lon_i in 1:length(lon_sel)) {
#       for (lat_i in 1:length(lat_sel)) {
#         if (all(!is.na(variable[lon_i, lat_i, ]))){
#           ts_pt <- ts(variable[lon_i, lat_i, ], start=c(year[1], mon[1]), end=c(year[length(year)], mon[length(mon)]), frequency=12)
#           spei3_pt = SPEI::spei(ts_pt, scale = 3, ref.start=c(1981,1), ref.end=c(2020,5))
#           spei3_gridded[lon_i, lat_i, ] = spei3_pt$fitted
#         }
#       }
#     }
#     
#     # write to output netcdf file
#     #===============================================================================
#     #----------------
#     # Make dimensions
#     #----------------
#     londim <- ncdim_def( name = 'lon', units = 'degrees_east', vals = lon_sel, longname = "Longitude")
#     latdim <- ncdim_def( name = 'lat', units = 'degrees_north', vals = lat_sel, longname = "Latitude" )
#     timedim <- ncdim_def( name = 'time', units = time_atts$units, calendar = time_atts$calendar, vals = time)
#     
#     #---------
#     # Make var
#     #---------
#     spei3_ncvar <- ncvar_def( name = 'SPEI3', units = 'standardised units', list(londim, latdim, timedim), missval = NA, prec="double")
#     
#     #---------------------
#     # Make new output file
#     #---------------------
#     ncid_out <- nc_create( outfile, vars = spei3_ncvar) #, force_v4 = TRUE)
#     
#     #-------------------------------
#     # Put some test data in the file
#     #-------------------------------
#     ncvar_put( nc = ncid_out, varid = spei3_ncvar, vals = spei3_gridded, start=NA, count=NA, verbose=TRUE )
#     nc_close(ncid_out)
#     rm(variable, spei3_gridded, spei3_ncvar)
#   }
# }

lat_start = which(lat == -44)
lat_end = which(lat == -22)
lat_count = lat_end - lat_start + 1

lon_start = which(lon == 135)
lon_end = which(lon == 154)
lon_count = lon_end - lon_start + 1

# set the baseline: 1980 to 2016 unless theres a reason to deviate from this
baseline = c(1980, 2016)

ref.start = c(baseline[1],1)
ref.end = c(baseline[2],12)

outfile = paste0("SPEI3_awra_monthly_baseline_", baseline[0], "to",  baseline[1], ".nc")

# get the subsetted variable
infile <- nc_open(infile_name)
variable = ncvar_get(nc = infile, varid = varname, start = c(lon_start, lat_start, 1), count = c(lon_count, lat_count, -1))
lon_sel = ncvar_get(nc = infile, varid = "lon", start = c(lon_start), count = c(lon_count))
lat_sel = ncvar_get(nc = infile, varid = "lat", start = c(lat_start), count = c(lat_count))
nc_close(infile)

spei3_gridded = array(NA, dim(variable))
for (lon_i in 1:length(lon_sel)) {
  for (lat_i in 1:length(lat_sel)) {
    if (all(!is.na(variable[lon_i, lat_i, ]))){
      ts_pt <- ts(variable[lon_i, lat_i, ], start=c(year[1], mon[1]), end=c(year[length(year)], mon[length(mon)]), frequency=12)
      spei3_pt = SPEI::spei(ts_pt, scale = 3, ref.start=ref.start, ref.end=ref.end)
      spei3_gridded[lon_i, lat_i, ] = spei3_pt$fitted
    }
  }
}

# write to output netcdf file
#===============================================================================
#----------------
# Make dimensions
#----------------
londim <- ncdim_def( name = 'lon', units = 'degrees_east', vals = lon_sel, longname = "Longitude")
latdim <- ncdim_def( name = 'lat', units = 'degrees_north', vals = lat_sel, longname = "Latitude" )
timedim <- ncdim_def( name = 'time', units = time_atts$units, calendar = time_atts$calendar, vals = time)

#---------
# Make var
#---------
spei3_ncvar <- ncvar_def( name = 'SPEI3', units = 'standardised units', list(londim, latdim, timedim), missval = NA, prec="double")

#---------------------
# Make new output file
#---------------------
ncid_out <- nc_create( outfile, vars = spei3_ncvar) #, force_v4 = TRUE)

#-------------------------------
# Put data in the file
#-------------------------------
ncvar_put( nc = ncid_out, varid = spei3_ncvar, vals = spei3_gridded, start=NA, count=NA, verbose=TRUE )
nc_close(ncid_out)
rm(variable, spei3_gridded, spei3_ncvar)

  
  

  

  

  






