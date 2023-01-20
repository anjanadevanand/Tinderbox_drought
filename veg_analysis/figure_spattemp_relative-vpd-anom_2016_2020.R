setwd(file.path(here::here()))
pacman::p_load(tidyverse,data.table,lubridate,stars,tictoc,
  patchwork,viridis,scico,future,plantecophys)

# Load datasets ===================================================
vp <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_vapourpres_h15_mean_r005_monthly_1971_2020.nc", proxy=T)

zone <- sf::read_sf("drought_focusArea/smooth_polygon_drought_focusArea.shp")
st_crs(zone)
zone <- sf::st_transform(zone, st_crs(vp))
zone_bb <- st_bbox(c("xmin"=135,"ymin"=-39,"xmax"=154,"ymax"=-23)) %>% st_as_sfc() %>% 
  st_set_crs(4326)
zone_bb <- sf::st_transform(zone_bb,st_crs(zone))

oz_land_poly <- rnaturalearth::ne_states(
  country='Australia',
  returnclass = 'sf') %>% 
  select(admin) %>% 
  st_transform(., st_crs(zone_bb)) %>% 
  st_crop(., zone_bb) %>% 
  st_union()


svp <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_vapourpres_h15_mean_r005_monthly_1971_2020.nc")[oz_land_poly] %>% 
  filter(time >= ymd("1980-01-01")) %>% 
  filter(time < ymd("2020-01-01")) %>% 
  st_as_stars()
svp <- svp %>% 
  units::drop_units() %>% 
  set_names('vp')

svp2 <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/monthly/monthly/agcd_v1_vapourpres_h15_mean_r005_monthly_2020.nc") %>% 
  st_set_crs(st_crs(zone)) %>% 
  .[oz_land_poly] %>% 
  filter(time >= ymd("2020-01-01")) %>% 
  filter(time < ymd("2021-01-01")) %>% 
  st_as_stars() %>% 
  units::drop_units() %>% 
  set_names('vp')

svp <- c(svp,svp2)

stmax <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/agcd_v1_tmax_mean_r005_monthly_1910_2020.nc")[oz_land_poly] %>% 
  filter(time >= ymd("1980-01-01")) %>% 
  filter(time < ymd("2020-01-01")) %>% 
  st_as_stars()
stmax <- stmax %>% 
  units::drop_units() %>% 
  set_names('tmax')

stmax2 <- stars::read_stars("/media/sami/srifai-2tb/data_2tb/clim_grid/tmp_agcd/monthly/monthly/agcd_v1_tmax_mean_r005_monthly_2020.nc") %>% 
  st_set_crs(st_crs(zone)) %>% 
  .[oz_land_poly] %>% 
  filter(time >= ymd("1980-01-01")) %>% 
  filter(time < ymd("2021-01-01")) %>% 
  st_as_stars() %>% 
  units::drop_units() %>% 
  set_names('tmax')

stmax <- c(stmax,stmax2)

# merges and calc VPD ======================================
clim <- merge(as.data.table(stmax),as.data.table(svp),by=c("x","y","time"))
clim <- clim %>% 
  rename(date=time) %>% 
  drop_na()
clim[,`:=`(year=year(date),month=month(date))]
clim[,`:=`(vp = vp/10)]
clim[,`:=`(esat = bigleaf::Esat.slope(Tair=tmax)$Esat)]
clim[,`:=`(vpd = esat-vp)] # vpd in kPa

clim$vpd %>% hist(100)

# reference period selectd to be consistent with 
# the other RS based anoms
ref <- clim[date>=ymd("2002-01-01")][date<=ymd("2016-12-31")][,
  .(vpd_u = mean(vpd,na.rm=T)), 
  by=.(x,y,month)]

clim <- merge(clim,ref,by=c("x","y","month"))
clim[,`:=`(vpd_anom = vpd-vpd_u)]

# PLOT ============================================= 
oz_poly <- rnaturalearth::ne_states(
  country='Australia',
  returnclass = 'sf') %>% 
  select(admin) 

oz_poly <- st_crop(oz_poly, 
                   st_as_stars(unique(ref[,.(x,y)])) %>% 
                     st_bbox())

rel_anom <- clim[year%in%2016:2020][,
  .(rel_anom = 100*vpd_anom/vpd_u),by=.(year,x,y)]

library(cols4all)
# clim[year%in%c(2016:2019)][,.(vpd_anom = mean(vpd_anom)),by=.(x,y,year)] %>% 

p_out <- rel_anom %>% 
  ggplot(data=., aes(x,y,fill=rel_anom))+
  geom_sf(data=oz_poly, 
          inherit.aes = F,
          color="grey40")+
  geom_tile()+
  geom_sf(data=zone, 
          inherit.aes = F, 
          lwd=0.5,
          color='black', 
          fill='transparent')+
  coord_sf(expand = F, 
    crs=st_crs(4326))+
  scale_x_continuous(breaks=c(seq(135,152,by=4)))+
  # scale_fill_viridis_b(option='B', 
  #                      limits=c(quantile(anoms_max$lst_anom,0.01), 
  #                               quantile(anoms_max$lst_anom,0.99)), 
  #                      oob=scales::squish, 
  #                      n.breaks=7)+
  scale_fill_continuous_c4a_div(palette='okeeffe1', 
    mid=0, 
    limits=c(-50,50),
    oob=scales::squish,
    breaks=c(-50,-25,0,25,50),
    labels=c("≤ -50","-25","0","25","≥ 50"),
    reverse=T)+
  # scale_fill_scico(palette='vik',midpoint = 0, direction = 1,
  #                      # limits=c(quantile(clim$vpd_anom,0.01),
  #                      #          quantile(clim$vpd_anom,0.99)),
  #                  oob=scales::squish)+
  labs(x=NULL,y=NULL, 
       fill="Annual Mean Relative 15:00 VPD Anomaly (%)     ")+
  facet_grid(~year)+
  theme_linedraw()+
  theme(panel.grid=element_blank(), 
        legend.position = 'bottom', 
        legend.key.width = unit(1.75,'cm'), 
        legend.key.height = unit(0.15,'cm'), 
        strip.background = element_blank(),
        strip.text = element_text(color='black',face='bold',size=13))
ggsave(p_out, 
  filename="figures/figure_spattemp_relative-vpd-anom_2016_2020.png", 
       width=22*(1.25),
       height=8,
       units='cm',
  device=grDevices::png,
       dpi=350)
# END PLOT =================================================

# plot(oz_poly)
# oz_poly %>% sf::st_union() %>% plot
