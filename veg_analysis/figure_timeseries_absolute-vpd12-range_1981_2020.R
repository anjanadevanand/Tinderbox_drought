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


vpd_ma <- ref[,.(vpd_ma = mean(vpd_u))]


# Plot =====================================================
clim[,.(val = median(vpd,na.rm=T), 
          val_lo1 = quantile(vpd,0.25,na.rm=T),
        val_hi1 = quantile(vpd,0.75,na.rm=T),
        val_lo2 = quantile(vpd,0.05,na.rm=T),
        val_hi2 = quantile(vpd,0.95,na.rm=T)),by=date] %>% 
  .[order(date)] %>% 
  .[,`:=`(vpd12 = frollmean(val,12,align = 'right'),
          vpd12_lo1 = frollmean(val_lo1,12,align='right'),
          vpd12_hi1 = frollmean(val_hi1,12,align='right'),
          vpd12_lo2 = frollmean(val_lo2,12,align='right'),
          vpd12_hi2 = frollmean(val_hi2,12,align='right'))] %>% 
  ggplot(data=.,aes(date,vpd12))+
  # geom_rect(aes(xmin=ymd("2017-01-01",tz = 'UTC'), 
  #               xmax=ymd("2019-12-31",tz='UTC'),
  #               ymin=0,ymax=4.15), 
  #           fill='grey90')+
  geom_ribbon(aes(x=date,ymin=vpd12_lo2,ymax=vpd12_hi2), 
              lty=0, fill="#ff0000",alpha=0.15)+
  geom_ribbon(aes(x=date,ymin=vpd12_lo1,ymax=vpd12_hi1), 
              lty=0, fill="#ff0000",alpha=0.3)+
  geom_hline(aes(yintercept=vpd_ma$vpd_ma),lty=2,col='grey30')+
  geom_hline(aes(yintercept=vpd_ma$vpd_ma*0.9),lty=3,col='grey30')+
  geom_hline(aes(yintercept=vpd_ma$vpd_ma*1.1),lty=3,col='grey30')+
  geom_line()+
  geom_line(data=. %>% 
              filter(date>=ymd("2002-01-01",tz='UTC')) %>% 
              filter(date<=ymd("2010-01-01",tz='UTC')), 
            color='#000000',lwd=1.5)+
  geom_line(data=. %>% 
              filter(date>=ymd("2002-01-01",tz='UTC')) %>% 
              filter(date<=ymd("2010-01-01",tz='UTC')), 
            color='#b30404',lwd=0.75)+
  geom_line(data=. %>% 
              filter(date>=ymd("2017-01-01",tz='UTC')), 
            color='#000000',lwd=1.5)+
  geom_line(data=. %>% 
              filter(date>=ymd("2017-01-01",tz='UTC')), 
            color='#ff6f00',lwd=0.75)+
  coord_cartesian(xlim=c(ymd("1980-01-01",tz='UTC'),
                         ymd("2020-09-01",tz='UTC')), 
                  ylim=c(0.99,4.15), 
                  expand=F)+
  labs(x=NULL,y=expression(paste(VPD["12 mo."]~(kPa))))+
    # manual legend for mean and 10% anom
    annotate('segment', 
    x=ymd("1981-01-01",tz="UTC"),
    xend=ymd("1984-01-01",tz="UTC"),
    y=4,yend=4,
    linetype=c(2),
    color='grey30')+
  annotate('segment', 
    x=ymd("1981-01-01",tz="UTC"),
    xend=ymd("1984-01-01",tz="UTC"),
    y=3.75,yend=3.75,
    linetype=c(3),
    color='grey30')+
  annotate('text', 
    x=ymd("1984-06-01",tz="UTC"),
    # xend=ymd("1989-01-01",tz="UTC"),
    y=4,yend=4,
    hjust=0,
    label='Mean')+
  annotate('text', 
    x=ymd("1984-06-01",tz="UTC"),
    # xend=ymd("1989-01-01",tz="UTC"),
    y=3.75,yend=3.75,
    hjust=0,
    label='10% Relative Anomaly')+
  
  # manual legend for MD and Tinderbox droughts
  annotate('segment', 
    x=ymd("2000-01-01",tz="UTC"),
    xend=ymd("2005-01-01",tz="UTC"),
    y=4,yend=4,
    linetype=c(1),
    lwd=1,
    color='#b30404')+
  annotate('segment', 
    x=ymd("2000-01-01",tz="UTC"),
    xend=ymd("2005-01-01",tz="UTC"),
    y=3.75,yend=3.75,
    linetype=c(1),
    lwd=1,
    color='#ff6f00')+
  annotate('text', 
    x=ymd("2006-01-01",tz="UTC"),
    # xend=ymd("1989-01-01",tz="UTC"),
    y=4,yend=4,
    hjust=0,
    label='Millenium Drought')+
  annotate('text', 
    x=ymd("2006-01-01",tz="UTC"),
    # xend=ymd("1989-01-01",tz="UTC"),
    y=3.75,yend=3.75,
    hjust=0,
    label='Tinderbox Drought')+
  theme_linedraw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        text = element_text(size=14))
  
ggsave("figures/supplement_figure_timeseries_absolute-vpd12-range_1981_2020.png", 
       width=22,
       height=8,
       units='cm',
       dpi=350)

