setwd(file.path(here::here()))
pacman::p_load(tidyverse,data.table,lubridate,stars,tictoc,
               patchwork,viridis,scico,future,plantecophys)

vp <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_vapourpres_h15_mean_r005_monthly_1971_2020.nc", proxy=F)
vp2 <- stars::read_stars("/home/sami/scratch/agcd_addendum/agcd_addendum/agcd_2020to22/monthly/agcd_vapourpres_h15_monmean_2020_2022.nc",proxy=F)

j1 <- vp %>% filter(time==ymd_hms("2020-01-16 00:00:00 UTC")) %>% st_as_stars()
j2 <- vp2 %>% filter(time==ymd_hms("2020-01-16 09:00:00 UTC")) %>% st_as_stars()





zone <- sf::read_sf("drought_focusArea/smooth_polygon_drought_focusArea.shp")
st_crs(zone)
zone <- sf::st_transform(zone, st_crs(vp))

svp <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_vapourpres_h15_mean_r005_monthly_1971_2020.nc")[zone] %>% 
  filter(time >= ymd("1980-01-01")) %>% 
  filter(time < ymd("2020-01-01")) %>% 
  st_as_stars()
svp <- svp %>% 
  units::drop_units() %>% 
  set_names('vp')

stmax <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_tmax_mean_r005_monthly_1910_2020.nc")[zone] %>% 
  filter(time >= ymd("1980-01-01")) %>% 
  filter(time < ymd("2020-01-01")) %>% 
  st_as_stars()
stmax <- stmax %>% 
  units::drop_units() %>% 
  set_names('tmax')

clim <- merge(as.data.table(stmax),as.data.table(svp),by=c("x","y","time"))
clim <- clim %>% 
  rename(date=time) %>% 
  drop_na()
clim[,`:=`(year=year(date),month=month(date))]
clim[,`:=`(vp = vp/10)]
clim[,`:=`(esat = bigleaf::Esat.slope(Tair=tmax)$Esat)]
clim[,`:=`(vpd = esat-vp)]

clim$vpd %>% hist(100)

ref <- clim[date>=ymd("2002-01-01")][date<=ymd("2016-12-31")][,
                                                              .(vpd_u = mean(vpd,na.rm=T)), 
                                                              by=.(x,y,month)]

clim <- merge(clim,ref,by=c("x","y","month"))
clim[,`:=`(vpd_anom = vpd-vpd_u)]



# Plot =====================================================
clim[,.(val = median(vpd,na.rm=T), 
        val_lo = quantile(vpd,0.05,na.rm=T),
        val_hi = quantile(vpd,0.95,na.rm=T)),by=date] %>% 
  .[order(date)] %>% 
  .[,`:=`(vpd12 = frollmean(val,12,align = 'right'), 
          vpd12_lo = frollmean(val_lo,12,align='right'),
          vpd12_hi = frollmean(val_hi,12,align='right'))] %>% 
  ggplot(data=.,aes(date,vpd12))+
  # geom_rect(aes(xmin=ymd("2017-01-01",tz = 'UTC'), 
  #               xmax=ymd("2019-12-31",tz='UTC'),
  #               ymin=0,ymax=4.15), 
  #           fill='grey90')+
  geom_ribbon(aes(x=date,ymin=vpd12_lo,ymax=vpd12_hi), 
              lty=0, fill="#cf0000",alpha=0.25)+
  geom_line()+
  geom_line(data=. %>% 
              filter(date>=ymd("2002-01-01",tz='UTC')) %>% 
              filter(date<=ymd("2010-01-01",tz='UTC')), 
            color='#000000',lwd=1.5)+
  geom_line(data=. %>% 
              filter(date>=ymd("2002-01-01",tz='UTC')) %>% 
              filter(date<=ymd("2010-01-01",tz='UTC')), 
            color='#cf0000',lwd=0.75)+
  geom_line(data=. %>% 
              filter(date>=ymd("2017-01-01",tz='UTC')), 
            color='#000000',lwd=1.5)+
  geom_line(data=. %>% 
              filter(date>=ymd("2017-01-01",tz='UTC')), 
            color='#cf0000',lwd=0.75)+
  coord_cartesian(xlim=c(ymd("1980-01-01",tz='UTC'),
                         ymd("2020-09-01",tz='UTC')), 
                  ylim=c(0.99,4.15), 
                  expand=F)+
  labs(x=NULL,y=expression(paste(VPD["12 mo."]~(kPa))))+
  theme_linedraw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        text = element_text(size=14))
  
ggsave("figures/figure_timeseries_absolute-vpd12-range_1981_2020.png", 
       width=22,
       height=8,
       units='cm',
       dpi=350)

