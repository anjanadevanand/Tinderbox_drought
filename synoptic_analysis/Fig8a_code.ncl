
;;********************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

;;********************************

begin
;;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

;; Output PDF plot? -set to True or False
;;---------------

doplot = True

;;***************
;; CLIMATOLOGY 1980-2016
;;***************

;; Read in rainfall climatology data 1980-2016
;;---------------

infile2  = "Daily_rainfall_area_1980_2016.nc"
inf2     = addfile(infile2,"r")
precip_2 = inf2->precip

time_2  = inf2->time
date_2  = cd_calendar(time_2,0)
year_2  = date_2(:,0)
month_2 = date_2(:,1)

lat      = inf2->lat
lat_dims = dimsizes(lat)

lon      = inf2->lon
lon_dims = dimsizes(lon)

;; Seasonal rainfall total accumulations for climatology
;;---------------

djf_ind_2 = ind(month_2.eq.12.or.month_2.eq.1.or.month_2.eq.2)
mam_ind_2 = ind(month_2.eq.3.or.month_2.eq.4.or.month_2.eq.5)
jja_ind_2 = ind(month_2.eq.6.or.month_2.eq.7.or.month_2.eq.8)
son_ind_2 = ind(month_2.eq.9.or.month_2.eq.10.or.month_2.eq.11)

precip_djf_2 = precip_2(djf_ind_2,:,:)
precip_mam_2 = precip_2(mam_ind_2,:,:)
precip_jja_2 = precip_2(jja_ind_2,:,:)
precip_son_2 = precip_2(son_ind_2,:,:)

precip_clim_accum_djf_2 = dim_sum_n_Wrap(precip_djf_2,0)
precip_clim_accum_mam_2 = dim_sum_n_Wrap(precip_mam_2,0)
precip_clim_accum_jja_2 = dim_sum_n_Wrap(precip_jja_2,0)
precip_clim_accum_son_2 = dim_sum_n_Wrap(precip_son_2,0)

;; 90th percentile at each grid point for days >= 0.01 mm/day for climatology
;;---------------

precip_90 = precip_2(0,:,:)
precip_90 = precip_2@_FillValue

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp      = precip_2(:,ii,jj)
    if (.not.all(ismissing(temp))) then
      valid_ind = ind(.not.ismissing(temp))
      temp_1    = temp(valid_ind)
      rain_ind  = ind(temp_1.ge.0.01)
      temp_2    = temp_1(rain_ind)
      opt       = False
      statx     = stat_dispersion(temp_2,opt)
      precip_90(ii,jj) = (/statx(13)/) ; upper decile
      delete([/temp_1,valid_ind,rain_ind,temp_2,statx/])
    end if
    delete([/temp/])
  end do
end do

;; Accumulated rainfall total for days > 90th percentile by season for climatology
;;---------------

precip_bins_accum_djf_2 = precip_2(0,:,:)
precip_bins_accum_djf_2 = 0.0
precip_bins_accum_mam_2 = precip_bins_accum_djf_2
precip_bins_accum_jja_2 = precip_bins_accum_djf_2
precip_bins_accum_son_2 = precip_bins_accum_djf_2

;; DJF
do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf_2(:,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; MAM
do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam_2(:,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; JJA
do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja_2(:,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; SON
do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son_2(:,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; Percentage contribution to accumulation of >90p events for 1980-2016
;;---------------

precip_clim_accum_contr_djf_2 = precip_clim_accum_djf_2
precip_clim_accum_contr_djf_2 = (precip_bins_accum_djf_2 /precip_clim_accum_djf_2) * 100.0
precip_clim_accum_contr_mam_2 = precip_clim_accum_mam_2
precip_clim_accum_contr_mam_2 = (precip_bins_accum_mam_2 /precip_clim_accum_mam_2) * 100.0
precip_clim_accum_contr_jja_2 = precip_clim_accum_jja_2
precip_clim_accum_contr_jja_2 = (precip_bins_accum_jja_2 /precip_clim_accum_jja_2) * 100.0
precip_clim_accum_contr_son_2 = precip_clim_accum_son_2
precip_clim_accum_contr_son_2 = (precip_bins_accum_son_2 /precip_clim_accum_son_2) * 100.0

;; Statitics of >90p events contribution
;;---------------

accum_contr_stats_djf_2 = new(6,float)
accum_contr_stats_mam_2 = accum_contr_stats_djf_2
accum_contr_stats_jja_2 = accum_contr_stats_djf_2
accum_contr_stats_son_2 = accum_contr_stats_djf_2

;; DJF
temp_djf    = precip_clim_accum_contr_djf_2
temp_djf_1d = ndtooned(temp_djf)
opt = False
statx = stat_dispersion(temp_djf_1d,opt)
accum_contr_stats_djf_2(0) = (/statx(24)/) ; 5th perc
accum_contr_stats_djf_2(1) = (/statx(6)/)  ; 25th perc
accum_contr_stats_djf_2(2) = (/statx(8)/)  ; median
accum_contr_stats_djf_2(3) = (/statx(10)/) ; 75th perc
accum_contr_stats_djf_2(4) = (/statx(25)/) ; 95th perc
accum_contr_stats_djf_2(5) = (/statx(0)/)  ; mean
delete([/statx/])

;; MAM
temp_mam    = precip_clim_accum_contr_mam_2
temp_mam_1d = ndtooned(temp_mam)
opt = False
statx = stat_dispersion(temp_mam_1d,opt)
accum_contr_stats_mam_2(0) = (/statx(24)/)
accum_contr_stats_mam_2(1) = (/statx(6)/)
accum_contr_stats_mam_2(2) = (/statx(8)/)
accum_contr_stats_mam_2(3) = (/statx(10)/)
accum_contr_stats_mam_2(4) = (/statx(25)/)
accum_contr_stats_mam_2(5) = (/statx(0)/)
delete([/statx/])

;; JJA
temp_jja    = precip_clim_accum_contr_jja_2
temp_jja_1d = ndtooned(temp_jja)
opt = False
statx = stat_dispersion(temp_jja_1d,opt)
accum_contr_stats_jja_2(0) = (/statx(24)/)
accum_contr_stats_jja_2(1) = (/statx(6)/)
accum_contr_stats_jja_2(2) = (/statx(8)/) 
accum_contr_stats_jja_2(3) = (/statx(10)/)
accum_contr_stats_jja_2(4) = (/statx(25)/)
accum_contr_stats_jja_2(5) = (/statx(0)/)
delete([/statx/])

;; SON
temp_son    = precip_clim_accum_contr_son_2
temp_son_1d = ndtooned(temp_son)
opt = False
statx = stat_dispersion(temp_son_1d,opt)
accum_contr_stats_son_2(0) = (/statx(24)/) ; 5th perc
accum_contr_stats_son_2(1) = (/statx(6)/)  ; 25th perc
accum_contr_stats_son_2(2) = (/statx(8)/)  ; median
accum_contr_stats_son_2(3) = (/statx(10)/) ; 75th perc
accum_contr_stats_son_2(4) = (/statx(25)/) ; 95th perc
accum_contr_stats_son_2(5) = (/statx(0)/)  ; mean
delete([/statx/])


;;***************
;; CASE STUDY 2015-2020
;;***************

;; Read in rainfall data 2015-2020
;;---------------

infile1 = "Daily_rainfall_area_2015_2020.nc"
inf1    = addfile(infile1,"r")
precip  = inf1->precip

time  = inf1->time 
date  = cd_calendar(time,0)
year  = date(:,0)
month = date(:,1)

;; Seasonal rainfall by year
;;---------------

djf_ind = ind(month.eq.12.or.month.eq.1.or.month.eq.2)
mam_ind = ind(month.eq.3.or.month.eq.4.or.month.eq.5)
jja_ind = ind(month.eq.6.or.month.eq.7.or.month.eq.8)
son_ind = ind(month.eq.9.or.month.eq.10.or.month.eq.11)

precip_djf = precip(djf_ind,:,:)
precip_mam = precip(mam_ind,:,:)
precip_jja = precip(jja_ind,:,:)
precip_son = precip(son_ind,:,:)

;; DJF
time_djf = precip_djf&time
date_djf = cd_calendar(time_djf,-2)
case_ind_djf_2016 = ind(date_djf.ge.20151201.and.date_djf.le.20160229)
case_ind_djf_2017 = ind(date_djf.ge.20161201.and.date_djf.le.20170228)
case_ind_djf_2018 = ind(date_djf.ge.20171201.and.date_djf.le.20180228)
case_ind_djf_2019 = ind(date_djf.ge.20181201.and.date_djf.le.20190228)
case_ind_djf_2020 = ind(date_djf.ge.20191201.and.date_djf.le.20200229)

;; MAM
time_mam = precip_mam&time
date_mam = cd_calendar(time_mam,-2)
case_ind_mam_2016 = ind(date_mam.ge.20160301.and.date_mam.le.20160531)
case_ind_mam_2017 = ind(date_mam.ge.20170301.and.date_mam.le.20170531)
case_ind_mam_2018 = ind(date_mam.ge.20180301.and.date_mam.le.20180531)
case_ind_mam_2019 = ind(date_mam.ge.20190301.and.date_mam.le.20190531)
case_ind_mam_2020 = ind(date_mam.ge.20200301.and.date_mam.le.20200531)

;; JJA
time_jja = precip_jja&time
date_jja = cd_calendar(time_jja,-2)
case_ind_jja_2016 = ind(date_jja.ge.20160601.and.date_jja.le.20160831)
case_ind_jja_2017 = ind(date_jja.ge.20170601.and.date_jja.le.20170831)
case_ind_jja_2018 = ind(date_jja.ge.20180601.and.date_jja.le.20180831)
case_ind_jja_2019 = ind(date_jja.ge.20190601.and.date_jja.le.20190831)
case_ind_jja_2020 = ind(date_jja.ge.20200601.and.date_jja.le.20200831)

;; SON
time_son = precip_son&time
date_son = cd_calendar(time_son,-2)
case_ind_son_2016 = ind(date_son.ge.20160901.and.date_son.le.20161130)
case_ind_son_2017 = ind(date_son.ge.20170901.and.date_son.le.20171130)
case_ind_son_2018 = ind(date_son.ge.20180901.and.date_son.le.20181130)
case_ind_son_2019 = ind(date_son.ge.20190901.and.date_son.le.20191130)
case_ind_son_2020 = ind(date_son.ge.20200901.and.date_son.le.20201130)

;; Seasonal rainfall accumulation by year
;;---------------

;; DJF
precip_accum_djf_2016 = dim_sum_n_Wrap(precip_djf(case_ind_djf_2016,:,:),0)
precip_accum_djf_2017 = dim_sum_n_Wrap(precip_djf(case_ind_djf_2017,:,:),0)
precip_accum_djf_2018 = dim_sum_n_Wrap(precip_djf(case_ind_djf_2018,:,:),0)
precip_accum_djf_2019 = dim_sum_n_Wrap(precip_djf(case_ind_djf_2019,:,:),0)
precip_accum_djf_2020 = dim_sum_n_Wrap(precip_djf(case_ind_djf_2020,:,:),0)

;; MAM
precip_accum_mam_2016 = dim_sum_n_Wrap(precip_mam(case_ind_mam_2016,:,:),0)
precip_accum_mam_2017 = dim_sum_n_Wrap(precip_mam(case_ind_mam_2017,:,:),0)
precip_accum_mam_2018 = dim_sum_n_Wrap(precip_mam(case_ind_mam_2018,:,:),0)
precip_accum_mam_2019 = dim_sum_n_Wrap(precip_mam(case_ind_mam_2019,:,:),0)
precip_accum_mam_2020 = dim_sum_n_Wrap(precip_mam(case_ind_mam_2020,:,:),0)

;; JJA
precip_accum_jja_2016 = dim_sum_n_Wrap(precip_jja(case_ind_jja_2016,:,:),0)
precip_accum_jja_2017 = dim_sum_n_Wrap(precip_jja(case_ind_jja_2017,:,:),0)
precip_accum_jja_2018 = dim_sum_n_Wrap(precip_jja(case_ind_jja_2018,:,:),0)
precip_accum_jja_2019 = dim_sum_n_Wrap(precip_jja(case_ind_jja_2019,:,:),0)
precip_accum_jja_2020 = dim_sum_n_Wrap(precip_jja(case_ind_jja_2020,:,:),0)

;; SON
precip_accum_son_2016 = dim_sum_n_Wrap(precip_son(case_ind_son_2016,:,:),0)
precip_accum_son_2017 = dim_sum_n_Wrap(precip_son(case_ind_son_2017,:,:),0)
precip_accum_son_2018 = dim_sum_n_Wrap(precip_son(case_ind_son_2018,:,:),0)
precip_accum_son_2019 = dim_sum_n_Wrap(precip_son(case_ind_son_2019,:,:),0)
precip_accum_son_2020 = dim_sum_n_Wrap(precip_son(case_ind_son_2020,:,:),0)

;; Accumulated rainfall total for days > 90th percentile by season
;;---------------

;; DJF
precip_bins_accum_djf_2016 = precip(0,:,:)
precip_bins_accum_djf_2016 = 0.0
precip_bins_accum_djf_2017 = precip_bins_accum_djf_2016
precip_bins_accum_djf_2018 = precip_bins_accum_djf_2016
precip_bins_accum_djf_2019 = precip_bins_accum_djf_2016
precip_bins_accum_djf_2020 = precip_bins_accum_djf_2016

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf(case_ind_djf_2016,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2016(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf(case_ind_djf_2017,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2017(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf(case_ind_djf_2018,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2018(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf(case_ind_djf_2019,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2019(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_djf(case_ind_djf_2020,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_djf_2020(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; MAM
precip_bins_accum_mam_2016 = precip(0,:,:)
precip_bins_accum_mam_2016 = 0.0
precip_bins_accum_mam_2017 = precip_bins_accum_mam_2016
precip_bins_accum_mam_2018 = precip_bins_accum_mam_2016
precip_bins_accum_mam_2019 = precip_bins_accum_mam_2016
precip_bins_accum_mam_2020 = precip_bins_accum_mam_2016

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam(case_ind_mam_2016,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2016(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam(case_ind_mam_2017,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2017(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam(case_ind_mam_2018,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2018(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam(case_ind_mam_2019,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2019(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_mam(case_ind_mam_2020,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_mam_2020(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; JJA
precip_bins_accum_jja_2016 = precip(0,:,:)
precip_bins_accum_jja_2016 = 0.0
precip_bins_accum_jja_2017 = precip_bins_accum_jja_2016
precip_bins_accum_jja_2018 = precip_bins_accum_jja_2016
precip_bins_accum_jja_2019 = precip_bins_accum_jja_2016
precip_bins_accum_jja_2020 = precip_bins_accum_jja_2016

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja(case_ind_jja_2016,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2016(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja(case_ind_jja_2017,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2017(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja(case_ind_jja_2018,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2018(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja(case_ind_jja_2019,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2019(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_jja(case_ind_jja_2020,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_jja_2020(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; SON
precip_bins_accum_son_2016 = precip(0,:,:)
precip_bins_accum_son_2016 = 0.0
precip_bins_accum_son_2017 = precip_bins_accum_son_2016
precip_bins_accum_son_2018 = precip_bins_accum_son_2016
precip_bins_accum_son_2019 = precip_bins_accum_son_2016
precip_bins_accum_son_2020 = precip_bins_accum_son_2016

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son(case_ind_son_2016,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2016(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son(case_ind_son_2017,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2017(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son(case_ind_son_2018,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2018(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son(case_ind_son_2019,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2019(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

do ii = 0,lat_dims-1
  do jj = 0,lon_dims-1
    temp = precip_son(case_ind_son_2020,ii,jj)
    ind1 := ind(temp.ge.precip_90(ii,jj))
    if (.not.all(ismissing(ind1))) then
      precip_bins_accum_son_2020(ii,jj) = sum(temp(ind1))
    end if
    delete([/temp/])
  end do
end do

;; Percentage contribution to accumulation of >90p events
;;---------------

;; DJF
precip_bins_accum_contr_djf_2016 = precip_bins_accum_djf_2016
precip_bins_accum_contr_djf_2016 = (precip_bins_accum_djf_2016/precip_accum_djf_2016) * 100.0
precip_bins_accum_contr_djf_2017 = precip_bins_accum_djf_2017
precip_bins_accum_contr_djf_2017 = (precip_bins_accum_djf_2017/precip_accum_djf_2017) * 100.0
precip_bins_accum_contr_djf_2018 = precip_bins_accum_djf_2018
precip_bins_accum_contr_djf_2018 = (precip_bins_accum_djf_2018/precip_accum_djf_2018) * 100.0
precip_bins_accum_contr_djf_2019 = precip_bins_accum_djf_2019
precip_bins_accum_contr_djf_2019 = (precip_bins_accum_djf_2019/precip_accum_djf_2019) * 100.0
precip_bins_accum_contr_djf_2020 = precip_bins_accum_djf_2020
precip_bins_accum_contr_djf_2020 = (precip_bins_accum_djf_2020/precip_accum_djf_2020) * 100.0

;; MAM
precip_bins_accum_contr_mam_2016 = precip_bins_accum_mam_2016
precip_bins_accum_contr_mam_2016 = (precip_bins_accum_mam_2016/precip_accum_mam_2016) * 100.0
precip_bins_accum_contr_mam_2017 = precip_bins_accum_mam_2017
precip_bins_accum_contr_mam_2017 = (precip_bins_accum_mam_2017/precip_accum_mam_2017) * 100.0
precip_bins_accum_contr_mam_2018 = precip_bins_accum_mam_2018
precip_bins_accum_contr_mam_2018 = (precip_bins_accum_mam_2018/precip_accum_mam_2018) * 100.0
precip_bins_accum_contr_mam_2019 = precip_bins_accum_mam_2019
precip_bins_accum_contr_mam_2019 = (precip_bins_accum_mam_2019/precip_accum_mam_2019) * 100.0
precip_bins_accum_contr_mam_2020 = precip_bins_accum_mam_2020
precip_bins_accum_contr_mam_2020 = (precip_bins_accum_mam_2020/precip_accum_mam_2020) * 100.0

;; JJA
precip_bins_accum_contr_jja_2016 = precip_bins_accum_jja_2016
precip_bins_accum_contr_jja_2016 = (precip_bins_accum_jja_2016/precip_accum_jja_2016) * 100.0
precip_bins_accum_contr_jja_2017 = precip_bins_accum_jja_2017
precip_bins_accum_contr_jja_2017 = (precip_bins_accum_jja_2017/precip_accum_jja_2017) * 100.0
precip_bins_accum_contr_jja_2018 = precip_bins_accum_jja_2018
precip_bins_accum_contr_jja_2018 = (precip_bins_accum_jja_2018/precip_accum_jja_2018) * 100.0
precip_bins_accum_contr_jja_2019 = precip_bins_accum_jja_2019
precip_bins_accum_contr_jja_2019 = (precip_bins_accum_jja_2019/precip_accum_jja_2019) * 100.0
precip_bins_accum_contr_jja_2020 = precip_bins_accum_jja_2020
precip_bins_accum_contr_jja_2020 = (precip_bins_accum_jja_2020/precip_accum_jja_2020) * 100.0

;; SON
precip_bins_accum_contr_son_2016 = precip_bins_accum_son_2016
precip_bins_accum_contr_son_2016 = (precip_bins_accum_son_2016/precip_accum_son_2016) * 100.0
precip_bins_accum_contr_son_2017 = precip_bins_accum_son_2017
precip_bins_accum_contr_son_2017 = (precip_bins_accum_son_2017/precip_accum_son_2017) * 100.0
precip_bins_accum_contr_son_2018 = precip_bins_accum_son_2018
precip_bins_accum_contr_son_2018 = (precip_bins_accum_son_2018/precip_accum_son_2018) * 100.0
precip_bins_accum_contr_son_2019 = precip_bins_accum_son_2019
precip_bins_accum_contr_son_2019 = (precip_bins_accum_son_2019/precip_accum_son_2019) * 100.0
precip_bins_accum_contr_son_2020 = precip_bins_accum_son_2020
precip_bins_accum_contr_son_2020 = (precip_bins_accum_son_2020/precip_accum_son_2020) * 100.0

;; Accumulation anomalies - event minus climatology
;;---------------

;; DJF
precip_bins_accum_anom_djf_2016 = precip_bins_accum_contr_djf_2016
precip_bins_accum_anom_djf_2016 = precip_bins_accum_contr_djf_2016 - precip_clim_accum_contr_djf_2
precip_bins_accum_anom_djf_2017 = precip_bins_accum_contr_djf_2017
precip_bins_accum_anom_djf_2017 = precip_bins_accum_contr_djf_2017 - precip_clim_accum_contr_djf_2
precip_bins_accum_anom_djf_2018 = precip_bins_accum_contr_djf_2018
precip_bins_accum_anom_djf_2018 = precip_bins_accum_contr_djf_2018 - precip_clim_accum_contr_djf_2
precip_bins_accum_anom_djf_2019 = precip_bins_accum_contr_djf_2019
precip_bins_accum_anom_djf_2019 = precip_bins_accum_contr_djf_2019 - precip_clim_accum_contr_djf_2
precip_bins_accum_anom_djf_2020 = precip_bins_accum_contr_djf_2020
precip_bins_accum_anom_djf_2020 = precip_bins_accum_contr_djf_2020 - precip_clim_accum_contr_djf_2

;; MAM
precip_bins_accum_anom_mam_2016 = precip_bins_accum_contr_mam_2016
precip_bins_accum_anom_mam_2016 = precip_bins_accum_contr_mam_2016 - precip_clim_accum_contr_mam_2
precip_bins_accum_anom_mam_2017 = precip_bins_accum_contr_mam_2017
precip_bins_accum_anom_mam_2017 = precip_bins_accum_contr_mam_2017 - precip_clim_accum_contr_mam_2
precip_bins_accum_anom_mam_2018 = precip_bins_accum_contr_mam_2018
precip_bins_accum_anom_mam_2018 = precip_bins_accum_contr_mam_2018 - precip_clim_accum_contr_mam_2
precip_bins_accum_anom_mam_2019 = precip_bins_accum_contr_mam_2019
precip_bins_accum_anom_mam_2019 = precip_bins_accum_contr_mam_2019 - precip_clim_accum_contr_mam_2
precip_bins_accum_anom_mam_2020 = precip_bins_accum_contr_mam_2020
precip_bins_accum_anom_mam_2020 = precip_bins_accum_contr_mam_2020 - precip_clim_accum_contr_mam_2

;; JJA
precip_bins_accum_anom_jja_2016 = precip_bins_accum_contr_jja_2016
precip_bins_accum_anom_jja_2016 = precip_bins_accum_contr_jja_2016 - precip_clim_accum_contr_jja_2
precip_bins_accum_anom_jja_2017 = precip_bins_accum_contr_jja_2017
precip_bins_accum_anom_jja_2017 = precip_bins_accum_contr_jja_2017 - precip_clim_accum_contr_jja_2
precip_bins_accum_anom_jja_2018 = precip_bins_accum_contr_jja_2018
precip_bins_accum_anom_jja_2018 = precip_bins_accum_contr_jja_2018 - precip_clim_accum_contr_jja_2
precip_bins_accum_anom_jja_2019 = precip_bins_accum_contr_jja_2019
precip_bins_accum_anom_jja_2019 = precip_bins_accum_contr_jja_2019 - precip_clim_accum_contr_jja_2
precip_bins_accum_anom_jja_2020 = precip_bins_accum_contr_jja_2020
precip_bins_accum_anom_jja_2020 = precip_bins_accum_contr_jja_2020 - precip_clim_accum_contr_jja_2

;; SON
precip_bins_accum_anom_son_2016 = precip_bins_accum_contr_son_2016
precip_bins_accum_anom_son_2016 = precip_bins_accum_contr_son_2016 - precip_clim_accum_contr_son_2
precip_bins_accum_anom_son_2017 = precip_bins_accum_contr_son_2017
precip_bins_accum_anom_son_2017 = precip_bins_accum_contr_son_2017 - precip_clim_accum_contr_son_2
precip_bins_accum_anom_son_2018 = precip_bins_accum_contr_son_2018
precip_bins_accum_anom_son_2018 = precip_bins_accum_contr_son_2018 - precip_clim_accum_contr_son_2
precip_bins_accum_anom_son_2019 = precip_bins_accum_contr_son_2019
precip_bins_accum_anom_son_2019 = precip_bins_accum_contr_son_2019 - precip_clim_accum_contr_son_2
precip_bins_accum_anom_son_2020 = precip_bins_accum_contr_son_2020
precip_bins_accum_anom_son_2020 = precip_bins_accum_contr_son_2020 - precip_clim_accum_contr_son_2

;; Statitics of >90p accumulation anomalies
;;---------------

accum_anom_stats = new((/6,20/),float)
accum_anom_stats!0 = "stat"
accum_anom_stats!1 = "season"
stat = (/"5th","25th","median","75th","95th","mean"/)
accum_anom_stats&stat = stat
season = (/"DJF_2015/16","MAM_2016","JJA_2016","SON_2016","DJF_2016/17","MAM_2017","JJA_2017","SON_2017","DJF_2017/18","MAM_2018","JJA_2018","SON_2018","DJF_2018/19","MAM_2019","JJA_2019","SON_2019","DJF_2019/20","MAM_2020","JJA_2020","SON_2020"/)
accum_anom_stats&season = season

;; 2016
;; DJF
tempb    = precip_bins_accum_anom_djf_2016(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,0) = (/statx(24)/) ; 5th perc
accum_anom_stats(1,0) = (/statx(6)/)  ; 25th perc
accum_anom_stats(2,0) = (/statx(8)/)  ; median
accum_anom_stats(3,0) = (/statx(10)/) ; 75th perc
accum_anom_stats(4,0) = (/statx(25)/) ; 95th perc
accum_anom_stats(5,0) = (/statx(0)/)  ; mean
delete([/tempb,tempb_1d,statx/])

;; MAM
tempb    = precip_bins_accum_anom_mam_2016(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,1) = (/statx(24)/)
accum_anom_stats(1,1) = (/statx(6)/) 
accum_anom_stats(2,1) = (/statx(8)/) 
accum_anom_stats(3,1) = (/statx(10)/) 
accum_anom_stats(4,1) = (/statx(25)/)
accum_anom_stats(5,1) = (/statx(0)/) 
delete([/tempb,tempb_1d,statx/])

;; JJA
tempb    = precip_bins_accum_anom_jja_2016(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,2) = (/statx(24)/)
accum_anom_stats(1,2) = (/statx(6)/) 
accum_anom_stats(2,2) = (/statx(8)/) 
accum_anom_stats(3,2) = (/statx(10)/)
accum_anom_stats(4,2) = (/statx(25)/)
accum_anom_stats(5,2) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; SON
tempb    = precip_bins_accum_anom_son_2016(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,3) = (/statx(24)/) 
accum_anom_stats(1,3) = (/statx(6)/)  
accum_anom_stats(2,3) = (/statx(8)/) 
accum_anom_stats(3,3) = (/statx(10)/)
accum_anom_stats(4,3) = (/statx(25)/) 
accum_anom_stats(5,3) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; 2017
;; DJF
tempb    = precip_bins_accum_anom_djf_2017(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,4) = (/statx(24)/) 
accum_anom_stats(1,4) = (/statx(6)/)  
accum_anom_stats(2,4) = (/statx(8)/)  
accum_anom_stats(3,4) = (/statx(10)/) 
accum_anom_stats(4,4) = (/statx(25)/)
accum_anom_stats(5,4) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; MAM
tempb    = precip_bins_accum_anom_mam_2017(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,5) = (/statx(24)/) 
accum_anom_stats(1,5) = (/statx(6)/)  
accum_anom_stats(2,5) = (/statx(8)/) 
accum_anom_stats(3,5) = (/statx(10)/) 
accum_anom_stats(4,5) = (/statx(25)/) 
accum_anom_stats(5,5) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; JJA
tempb    = precip_bins_accum_anom_jja_2017(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,6) = (/statx(24)/) 
accum_anom_stats(1,6) = (/statx(6)/)  
accum_anom_stats(2,6) = (/statx(8)/)  
accum_anom_stats(3,6) = (/statx(10)/) 
accum_anom_stats(4,6) = (/statx(25)/) 
accum_anom_stats(5,6) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; SON
tempb    = precip_bins_accum_anom_son_2017(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,7) = (/statx(24)/) 
accum_anom_stats(1,7) = (/statx(6)/)  
accum_anom_stats(2,7) = (/statx(8)/)  
accum_anom_stats(3,7) = (/statx(10)/)
accum_anom_stats(4,7) = (/statx(25)/) 
accum_anom_stats(5,7) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; 2018
;; DJF
tempb    = precip_bins_accum_anom_djf_2018(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,8) = (/statx(24)/) 
accum_anom_stats(1,8) = (/statx(6)/)  
accum_anom_stats(2,8) = (/statx(8)/)  
accum_anom_stats(3,8) = (/statx(10)/) 
accum_anom_stats(4,8) = (/statx(25)/) 
accum_anom_stats(5,8) = (/statx(0)/) 
delete([/tempb,tempb_1d,statx/])

;; MAM
tempb    = precip_bins_accum_anom_mam_2018(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,9) = (/statx(24)/) 
accum_anom_stats(1,9) = (/statx(6)/)  
accum_anom_stats(2,9) = (/statx(8)/)  
accum_anom_stats(3,9) = (/statx(10)/) 
accum_anom_stats(4,9) = (/statx(25)/) 
accum_anom_stats(5,9) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; JJA
tempb    = precip_bins_accum_anom_jja_2018(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,10) = (/statx(24)/) 
accum_anom_stats(1,10) = (/statx(6)/) 
accum_anom_stats(2,10) = (/statx(8)/) 
accum_anom_stats(3,10) = (/statx(10)/)
accum_anom_stats(4,10) = (/statx(25)/) 
accum_anom_stats(5,10) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; SON
tempb    = precip_bins_accum_anom_son_2018(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,11) = (/statx(24)/)
accum_anom_stats(1,11) = (/statx(6)/)  
accum_anom_stats(2,11) = (/statx(8)/) 
accum_anom_stats(3,11) = (/statx(10)/) 
accum_anom_stats(4,11) = (/statx(25)/) 
accum_anom_stats(5,11) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; 2019
;; DJF
tempb    = precip_bins_accum_anom_djf_2019(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,12) = (/statx(24)/) 
accum_anom_stats(1,12) = (/statx(6)/) 
accum_anom_stats(2,12) = (/statx(8)/)  
accum_anom_stats(3,12) = (/statx(10)/) 
accum_anom_stats(4,12) = (/statx(25)/) 
accum_anom_stats(5,12) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; MAM
tempb    = precip_bins_accum_anom_mam_2019(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,13) = (/statx(24)/) 
accum_anom_stats(1,13) = (/statx(6)/)  
accum_anom_stats(2,13) = (/statx(8)/)  
accum_anom_stats(3,13) = (/statx(10)/) 
accum_anom_stats(4,13) = (/statx(25)/) 
accum_anom_stats(5,13) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; JJA
tempb    = precip_bins_accum_anom_jja_2019(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,14) = (/statx(24)/) 
accum_anom_stats(1,14) = (/statx(6)/)  
accum_anom_stats(2,14) = (/statx(8)/)  
accum_anom_stats(3,14) = (/statx(10)/) 
accum_anom_stats(4,14) = (/statx(25)/)
accum_anom_stats(5,14) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; SON
tempb    = precip_bins_accum_anom_son_2019(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,15) = (/statx(24)/) 
accum_anom_stats(1,15) = (/statx(6)/)  
accum_anom_stats(2,15) = (/statx(8)/) 
accum_anom_stats(3,15) = (/statx(10)/) 
accum_anom_stats(4,15) = (/statx(25)/) 
accum_anom_stats(5,15) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; 2020
;; DJF
tempb    = precip_bins_accum_anom_djf_2020(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,16) = (/statx(24)/) 
accum_anom_stats(1,16) = (/statx(6)/)  
accum_anom_stats(2,16) = (/statx(8)/)  
accum_anom_stats(3,16) = (/statx(10)/) 
accum_anom_stats(4,16) = (/statx(25)/) 
accum_anom_stats(5,16) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; MAM
tempb    = precip_bins_accum_anom_mam_2020(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,17) = (/statx(24)/) 
accum_anom_stats(1,17) = (/statx(6)/)  
accum_anom_stats(2,17) = (/statx(8)/) 
accum_anom_stats(3,17) = (/statx(10)/) 
accum_anom_stats(4,17) = (/statx(25)/) 
accum_anom_stats(5,17) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; JJA
tempb    = precip_bins_accum_anom_jja_2020(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,18) = (/statx(24)/) 
accum_anom_stats(1,18) = (/statx(6)/)  
accum_anom_stats(2,18) = (/statx(8)/)  
accum_anom_stats(3,18) = (/statx(10)/) 
accum_anom_stats(4,18) = (/statx(25)/) 
accum_anom_stats(5,18) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;; SON
tempb    = precip_bins_accum_anom_son_2020(:,:)
tempb_1d = ndtooned(tempb)
opt     = False
statx   = stat_dispersion(tempb_1d,opt)
accum_anom_stats(0,19) = (/statx(24)/) 
accum_anom_stats(1,19) = (/statx(6)/)  
accum_anom_stats(2,19) = (/statx(8)/)  
accum_anom_stats(3,19) = (/statx(10)/)
accum_anom_stats(4,19) = (/statx(25)/) 
accum_anom_stats(5,19) = (/statx(0)/)  
delete([/tempb,tempb_1d,statx/])

;;***************
;; PLOT BOX AND WHISKER PLOTS FOR ANOMALY STATISTICS AS PER FIG 8A  
;;***************

if (doplot) then

  wks = gsn_open_wks("pdf","Anomaly_statistics_8a")

  gsn_define_colormap(wks, "ncl_default")

  panel_plot = new(1,graphic)

  x1 = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/)

  season_str = (/"DJF","MAM","JJA","SON","DJF","MAM","JJA","SON","DJF","MAM","JJA","SON","DJF","MAM","JJA","SON","DJF","MAM","JJA","SON"/)

  res = True

  res@gsnDraw     = False
  res@gsnFrame    = False
  res@gsnMaximize = False ;True

  res@trXMinF = 0.0
  res@trXMaxF = 21.0
  res@trYMinF = -100.0
  res@trYMaxF = 100.0

  res@tiYAxisString         = "Anomalous proportion [%]"
  res@tiYAxisFontHeightF    = 0.01
  res@tiXAxisString         = ""
  res@tiXAxisFontHeightF    = 0.01

  res@tmXBLabelDeltaF      = 1.0
  res@tmXBLabelFontHeightF = 0.01
  res@tmXBMode             = "Explicit"
  res@tmXBLabels           = season_str
  res@tmXBLabelAngleF      = 90
  res@tmYLLabelFontHeightF = 0.01
  res@tmXBOn               = True
  res@tmXTOn               = False
  res@tmYLOn               = True
  res@tmYROn               = False
  res@tmXTBorderOn         = True
  res@tmYRBorderOn         = True

  res@vpXF      = 0.15
  res@vpYF      = 0.9
  res@vpWidthF  = 0.8
  res@vpHeightF = 0.3


  llres                  = True
  llres@gsLineThicknessF = 1.1 ;2.5

  opti                   = True
  opti@boxWidth          = 0.4 ;100.0 ;50.0 ;2.0 ;0.25
  opti@boxColors         = "blue"

  y_0a = accum_anom_stats(0:4,:)
  y_0 = y_0a(season|:,stat|:)

  plot_0 = boxplot(wks,x1,y_0,opti,res,llres)

  mres = True
  mres@gsMarkerIndex = 1
  mres@gsMarkerSizeF = 0.009
  dummy_0b = new(20,graphic)
  do ww = 0,19
    mres@gsMarkerColor = "blue"
    dummy_0b(ww) = gsn_add_polymarker(wks,plot_0,x1(ww),accum_anom_stats(5,ww),mres)
  end do

  polyres = True
  polyres@gsLineThicknessF  = 1.0
  polyres@gsLineColor       = "gray50"
  polyres@gsLineDashPattern = 1
  dummy_0 = gsn_add_polyline(wks,plot_0,(/0,res@trXMaxF/),(/0,0/),polyres)

  dummy_0c = new(4,graphic)
  dummy_0c(0) = gsn_add_polyline(wks,plot_0,(/4.5,4.5/),(/res@trYMinF,res@trYMaxF/),polyres)
  dummy_0c(1) = gsn_add_polyline(wks,plot_0,(/8.5,8.5/),(/res@trYMinF,res@trYMaxF/),polyres)
  dummy_0c(2) = gsn_add_polyline(wks,plot_0,(/12.5,12.5/),(/res@trYMinF,res@trYMaxF/),polyres)
  dummy_0c(3) = gsn_add_polyline(wks,plot_0,(/16.5,16.5/),(/res@trYMinF,res@trYMaxF/),polyres)

  txid = new(6,graphic)
  txres = True 
  txres@txFontHeightF = 0.01
  txid(0) = gsn_add_text(wks,plot_0,"2016",2.5,75,txres)
  txid(1) = gsn_add_text(wks,plot_0,"2017",6.5,75,txres)
  txid(2) = gsn_add_text(wks,plot_0,"2018",10.5,75,txres)
  txid(3) = gsn_add_text(wks,plot_0,"2019",14.5,75,txres)
  txid(4) = gsn_add_text(wks,plot_0,"2020",18.5,75,txres)

  txid(5) = gsn_add_text(wks,plot_0,"A",20.5,75,txres)

  panel_plot(0) = plot_0

  resP = True
  resP@gsn_maximize = True
  gsn_panel(wks,panel_plot,(/2,1/),resP)

end if

;;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
end