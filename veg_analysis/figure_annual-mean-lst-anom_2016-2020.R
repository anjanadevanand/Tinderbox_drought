pacman::p_load(tidyverse,data.table,lubridate,stars,scico,viridis,cols4all)

fl_ref <- list.files("../../data_general/proc_oz_npv/","medianLST_norms",full.names = T)
fl <- list.files("../../data_general/proc_oz_npv/","medianLST_2016_2020",full.names = T)

r <- stars::read_stars(fl_ref) %>% 
  set_names('lst_u') %>% 
  st_set_dimensions(., 3, values=1:12, names='month')

s <- stars::read_stars(fl) %>% 
  set_names('lst') %>% 
  st_set_dimensions(., 3, 
  values = seq(ymd("2016-01-01"),ymd("2020-12-01"),by='1 months'), 
  names = 'date')

zone <- sf::read_sf("drought_focusArea/smooth_polygon_drought_focusArea.shp")
zone_bb <- st_bbox(c("xmin"=135,"ymin"=-39,"xmax"=154,"ymax"=-23)) %>% st_as_sfc() %>% 
  st_set_crs(4326)

s <- s[zone_bb]
r <- r[zone_bb]


s <- s %>% 
  as.data.table()
s[,`:=`(year=year(date+months(0)), 
        month=month(date))]
s[,`:=`(season = case_when(
  month%in%c(12,1,2)~"DJF",
  month%in%c(3,4,5)~"MAM",
  month%in%c(6,7,8)~"JJA",
  month%in%c(9,10,11)~"SON"
))]
# s <- s[date<ymd("2019-12-01")][,`:=`(year=ifelse(year==2016,2017,year))]

r <- r %>% 
  as.data.table()
r[,`:=`(season = case_when(
  month%in%c(12,1,2)~"DJF",
  month%in%c(3,4,5)~"MAM",
  month%in%c(6,7,8)~"JJA",
  month%in%c(9,10,11)~"SON"
))]

dat <- merge(s,r,by=c("x","y","month","season"))
dat <- dat[is.na(lst_u)==F]

# PLOT ============================================= 
dat[,`:=`(idx = .GRP), by=.(x,y)]
dat_anoms <- dat[,`:=`(lst_anom = lst-lst_u)]
# anoms_max <- dat_anoms[dat_anoms[,(lst_anom==max(lst_anom)),by=.(x,y,year)]$V1]
lst_rel_anoms <- dat_anoms[,.(rel_lst_anom = 100*mean(lst_anom,na.rm = T)/lst_u, 
                          lst_anom = mean(lst_anom,na.rm=T)), 
                       by=.(x,y,year)]

oz_poly <- rnaturalearth::ne_states(
  country='Australia',
  returnclass = 'sf') %>% 
  select(admin) 

oz_poly <- st_crop(oz_poly, 
  zone_bb %>% 
                     st_bbox())

p_out <- lst_rel_anoms %>% 
  ggplot(data=., aes(x,y,fill=lst_anom))+
  geom_sf(data=oz_poly, 
          inherit.aes = F,
          color="grey40")+
  geom_raster()+
  geom_sf(data=zone, 
          inherit.aes = F, 
          lwd=0.5,
          color='black', 
          fill='transparent')+
  geom_sf(data=oz_poly,inherit.aes = F, 
          fill=NA)+
  coord_sf(expand = F, 
    crs = st_crs(4326))+
  scale_x_continuous(breaks=c(seq(138,152,by=4)))+
    scale_fill_continuous_c4a_div(palette='vik', 
    mid=0, 
    limits=c(-6,6),
    oob=scales::squish,
    breaks=c(-6,-3,0,3,6),
    labels=c("≤ -6","-3","0","3","≥ 6"),
    reverse=F)+
  # scale_fill_viridis_b(option='B', 
  #                      limits=c(quantile(anoms_max$lst_anom,0.01), 
  #                               quantile(anoms_max$lst_anom,0.99)), 
  #                      oob=scales::squish, 
  #                      n.breaks=7)+
  # scale_fill_scico(palette='vik',midpoint = 0, direction = 1,
  #                      # limits=c(quantile(anoms_max$lst_anom,0.01),
  #                      #          quantile(anoms_max$lst_anom,0.99)),
  #                  oob=scales::squish)+
  labs(x=NULL,y=NULL, 
       fill="(Annual)\n13:30 LST\nAnom. (°C)")+
       # fill="Annual Mean \n13:30 LST\nAnom. (°C)")+
  facet_wrap(~year,nrow = 1)+
  theme_linedraw()+
  theme(panel.grid=element_blank(), 
        legend.position = 'right', 
        legend.key.width = unit(0.2,'cm'), 
        legend.key.height = unit(0.6,'cm'), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=12),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color='black',face='bold',size=13))
p_out
ggsave(p_out, 
  filename = "figures/v2/figure_annual-mean-lst-anom_2016-2020.png", 
       width=22*(5/4),
       height=8*0.75,
       units='cm',
       device=grDevices::png,
       dpi=350)
ggsave(p_out, 
  filename = "figures/v2/figure_annual-mean-lst-anom_2016-2020.svg", 
       width=22*(5/4),
       height=8*0.75,
       units='cm',
       device=grDevices::svg,
       dpi=350)

# END PLOT =================================================


# r <- r[,.(lst_u = mean(lst_u,na.rm=T)), by=.(x,y,season)]
# 
# dat <- merge(s,r,by=c("x","y","season"))
# dat[,`:=`(lst_max_anom = lst-lst_u)]
# dat <- dat[x>135 & y< -23.5 & y> -39]
# dat[,`:=`(season = factor(season,levels = c("DJF","MAM","JJA","SON"),ordered = T))]
# 
# oz_poly <- rnaturalearth::ne_states(
#   country='Australia',
#   returnclass = 'sf') %>% 
#   select(admin) 
# 
# oz_poly <- st_crop(oz_poly, 
#   st_as_stars(unique(dat[,.(x,y)])) %>% 
#   st_bbox())
# 
# dat$lst_max_anom %>% quantile(., c(0.01,0.99),na.rm=T)
# 
# sf::sf_use_s2(F)
# dat[year<2020][is.na(lst_max_anom)==F] %>% 
#   ggplot(data=.,aes(x,y,
#     fill= lst_max_anom ))+
#   geom_sf(data=oz_poly,inherit.aes = F, 
#     fill='grey')+
#   geom_tile()+
#   geom_sf(data=oz_poly,inherit.aes = F, 
#     fill=NA)+
#   scale_fill_distiller(
#     type='div',direction = -1,
#     palette = 5,
#     limits=c(-15,15),
#     oob=scales::squish)+
#   coord_sf(expand=F,
#     crs=4326)+
#   scale_x_continuous(breaks=c(135,140,145,150))+
#   labs(x=NULL,
#     y=NULL,
#     title='MODIS AQUA Max LST Anomaly',
#     fill="(K)")+
#   facet_grid(season~year,drop = T)+
#   theme_linedraw()+
#   theme(panel.grid = element_blank(),
#     panel.background = element_rect(fill='lightblue'),
#     legend.key.height = unit(3.33,'cm'))
# ggsave("figs/fig_spattemp-LST-Max-anom_seasonal_2017-2019.png", 
#   width=20,
#   height=20,
#   units='cm',
#   dpi=350)
