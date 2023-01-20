# ****************************************************************************
# VOD ---------------------------------------------------------------------
# ****************************************************************************
pacman::p_load(tidyverse,data.table,lubridate,stars,scico,viridis,cols4all)

out2 <- arrow::read_parquet("../../data_general/AMSRE/se_oz_vod_2002_2016.parquet")
out3 <- arrow::read_parquet("../../data_general/AMSRE/se_oz_vod_2017_2020.parquet")

# [x>140 & y< -28 & y> -39]
out2 <- out2[x>135 & y< -22.5 & y> -39][band=='vod'&qa==0][,`:=`(hour=hour(date))][hour==13]
out3 <- out3[x>135 & y< -22.5 & y> -39][band=='vod'&qa==0][,`:=`(hour=hour(date))][hour==13]




out2[,`:=`(
  hour=hour(date),
  month=month(date))][,`:=`(season = case_when(
    month%in%c(12,1,2)~'DJF',
    month%in%c(3,4,5)~'MAM',
    month%in%c(6,7,8)~'JJA',
    month%in%c(9,10,11)~'SON'))][
      ,`:=`(season=factor(season,levels = c("DJF","MAM","JJA","SON"),ordered = T))
    ]

out3[,`:=`(
  hour=hour(date),
  month=month(date))][,`:=`(season = case_when(
    month%in%c(12,1,2)~'DJF',
    month%in%c(3,4,5)~'MAM',
    month%in%c(6,7,8)~'JJA',
    month%in%c(9,10,11)~'SON'))][
      ,`:=`(season=factor(season,levels = c("DJF","MAM","JJA","SON"),ordered = T))
    ]

# 2002-2016 reference period
ref <- out2[band=='vod'&qa==0][,
                               .(vod_mean = mean(value,na.rm=T), 
                                 vod_min = min(value,na.rm=T),
                                 vod_50 = median(value,na.rm=T),
                                 vod_05 = quantile(value,0.05,na.rm=T)),
                               by=.(x,y,season)]

rbindlist(list(out2[date >= ymd("2016-01-01")][date <= ymd("2019-11-30")],
out3[date >= ymd("2016-01-01")][date <= ymd("2019-11-30")]))

# big dry 2017- Aug 31 2019
bd <- rbindlist(
  list(out2[date >= ymd("2016-01-01")][date <= ymd("2021-01-01")],
      out3[date >= ymd("2016-01-01")][date <= ymd("2021-01-01")])) %>% 
.[band=='vod'&qa==0] %>% 
  .[,`:=`(hour=hour(date),
                                                                 month=month(date),
                                                                 year=year(date))] %>% .[,
                                                                                   .(vod_mean = mean(value,na.rm=T), 
                                                                                     vod_min = min(value,na.rm=T),
                                                                                     vod_50 = median(value,na.rm=T),
                                                                                     vod_05 = quantile(value,0.05,na.rm=T)),
                                                                                   by=.(x,y,season,year)]

# tmp1 <- merge(bd,ref,by=c("x","y","season"),suffixes = c("",".ref"))

rel_diff <- merge(bd,ref,by=c("x","y","season"), suffix=c("_bd","_ref")) %>% 
  .[,`:=`(vod_mean_diff = 100*(vod_mean_bd - vod_mean_ref)/vod_mean_ref  )]


# plotting ----------------------------------------------------------------
zone <- sf::read_sf("drought_focusArea/smooth_polygon_drought_focusArea.shp")
zone_bb <- st_bbox(c("xmin"=135,"ymin"=-39,"xmax"=154,"ymax"=-23)) %>% st_as_sfc() %>% 
  st_set_crs(4326)
oz_poly <- rnaturalearth::ne_states(
  country='Australia',
  returnclass = 'sf') %>% 
  select(admin) 

# oz_poly <- st_crop(oz_poly, 
#                    st_as_stars(unique(rel_diff[,.(x,y)])) %>% 
#                      st_bbox())
oz_poly <- st_crop(oz_poly, 
                     zone_bb %>% st_bbox())

oz_poly

rel_diff_son <- rel_diff[season=="SON"] %>% 
  st_as_sf(coords=c("x","y"))
st_crs(rel_diff_son) <- st_crs(4326)
rel_diff_son$x <- rel_diff[season=="SON"]$x
rel_diff_son$y <- rel_diff[season=="SON"]$y
rel_diff_son <- st_intersection(rel_diff_son,zone_bb)
# rel_diff_son <- rel_diff_son %>%
#   mutate(x=st_coordinates(.)[,'X'],
#     y=st_coordinates(.)[,'Y'])

rel_diff[season=="SON"]$vod_mean_diff %>% quantile(., c(0.01,0.99))

sf::sf_use_s2(F)

p_out <- rel_diff_son %>% 
  ggplot(data=., aes(x,y,fill=vod_mean_diff))+
  geom_sf(data=oz_poly,
          inherit.aes = F,
          color="grey40")+
  geom_tile()+
  geom_sf(data=zone, 
          inherit.aes = F, 
          lwd=0.5,
          color='black', 
          fill='transparent')+
  scale_fill_distiller(type='div',direction = 1,
                       limits=c(-50,50), 
                       oob=scales::squish)+
  coord_sf(expand=F,
           crs=4326)+
  scale_x_continuous(breaks=c(135,140,145,150))+
  labs(x=NULL,
       y=NULL,
       fill='Sept-Nov Vegetation Optical Depth Relative Anomaly (%)      ')+
  
  facet_wrap(~year,drop = T,nrow=1)+
  theme_linedraw()+
  theme(panel.grid=element_blank(), 
        legend.position = 'bottom', 
        legend.key.width = unit(1.75,'cm'), 
        legend.key.height = unit(0.15,'cm'), 
        strip.background = element_blank(),
        strip.text = element_text(color='black',face='bold',size=13))
ggsave(p_out, 
  filename = "figures/figure_relative-median-VOD-anom_SON_2016-2020.png", 
       width=22*(5/4),
       height=8,
       units='cm',
  device=grDevices::png,
       dpi=350)
