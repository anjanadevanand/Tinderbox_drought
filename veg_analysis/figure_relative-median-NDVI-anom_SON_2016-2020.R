pacman::p_load(tidyverse,data.table,lubridate,stars,scico,viridis,cols4all)
# dat <- arrow::read_parquet("/media/sami/srifai-ssd/data-ssd/AMSRE/se_oz_vod_2017_2020.parquet")
setwd("/home/sami/srifai@gmail.com/work/research/CLEX/CLEX_drought_synthesis/")


# fire mask ---------------------------------------------------------------
fm <- stars::read_stars("../../data_general/proc_oz_npv/blackSummerFiresMask.tif") %>% 
  set_names('bs_burned') %>% 
  as.data.table()

# Drought zone -----------------------------------------------------------
zone <- sf::read_sf("drought_focusArea/smooth_polygon_drought_focusArea.shp")
zone_bb <- st_bbox(c("xmin"=135,"ymin"=-39,"xmax"=153.65,"ymax"=-23)) %>% st_as_sfc() %>% 
  st_set_crs(4326)


ndvi_son_ref <- stars::read_stars("../../data_general/proc_oz_npv/myd13a2_c61_mean-annual-mean_ndvi_son_2002-2016.tif") %>% 
  # st_set_dimensions(., 3, values=2003:2016,names='year') %>% 
  # st_apply(., 1:2,mean,na.rm=T) %>% 
  set_names('ndvi') %>% 
  .[zone_bb] %>% 
  as.data.table() %>% 
  .[,`:=`(season='SON')] %>% 
  .[is.na(ndvi)==F]


ndvi_son <- stars::read_stars("../../data_general/proc_oz_npv/myd13a2_c61_mean_ndvi_son_2016-2020.tif") %>% 
  st_set_dimensions(., 3, values=2016:2020,names='year') %>% 
  set_names('ndvi') %>% 
  .[zone_bb] %>% 
  as.data.table() %>% 
  .[,`:=`(season='SON')] %>% 
  .[is.na(ndvi)==F] %>% 
  .[year>=2016]

dat_ndvi <- merge(ndvi_son_ref %>% rename(ndvi_u=ndvi),
ndvi_son, 
by=c("x","y","season"))

dat_ndvi[,`:=`(rel_anom = 100*(ndvi-ndvi_u)/ndvi_u)]
dat_ndvi[year>=2016 & year<=2019]$rel_anom %>% 
  quantile(.,c(0.05,0.95),na.rm=T)
dat_ndvi[year==2019][ndvi_u > 0] %>% summarize(val = mean(rel_anom))
dat_ndvi[year==2019][ndvi_u > 0] %>% summarize(val = sum(rel_anom < 0)/n())
dat_ndvi[year==2019][ndvi_u > 0.05]$rel_anom %>% summary

# dndvi <- dndvi[x>140 & y< -28 & y> -39]
# dndvi <- dndvi[x>135 & y< -23.5 & y> -39]
# dndvi[,`:=`(season = factor(season,levels = c("DJF","MAM","JJA","SON"),ordered = T))]

# Plotting ================================================
oz_poly <- rnaturalearth::ne_states(
  country='Australia',
  returnclass = 'sf') %>% 
  select(admin) 

# oz_poly <- st_crop(oz_poly, 
#                    st_as_stars(unique(dat_ndvi[,.(x,y)])) %>% 
#                      st_bbox())
oz_poly <- st_crop(oz_poly, zone_bb %>% st_bbox())


sf::sf_use_s2(F)
library(cols4all)
p_out <- merge(dat_ndvi,fm,by=c("x","y"))[year>=2016 & year<=2020] %>% 
  # .[!(year==2019 & season=="SON")] %>% 
  mutate(bs_burned = case_when((season=='SON')&(year==2019)&(bs_burned)==1 ~ 1, 
                               TRUE ~ NA_real_)) %>% #pull(bs_burned) %>% table
  ggplot(data=.,aes(x,y,
                    fill= 100*(ndvi-ndvi_u)/ndvi_u ))+
  geom_sf(data=oz_poly, 
          inherit.aes = F,
          color="grey40")+
  geom_tile()+
  geom_sf(data=zone, 
          inherit.aes = F, 
          lwd=0.5,
          color='black', 
          fill='transparent')+
  geom_point(inherit.aes = F,
             data=. %>% filter(bs_burned==1) %>% sample_frac(0.1),
             aes(x,y), 
             size=2,
             shape="x",
             col='#ff0000')+
  geom_sf(data=oz_poly,inherit.aes = F, 
          fill=NA)+
  scale_fill_continuous_c4a_div(palette='bam', 
    mid=0, 
    limits=c(-50,50),
    oob=scales::squish,
    breaks=c(-50,-25,0,25,50),
    labels=c("≤ -50","-25","0","25","≥ 50"),
    reverse=F)+
  # scale_fill_distiller(
  #   type='div',direction = 1, 
  #   limits=c(-60,60), 
  #   breaks = c(-60,-30,0,30,60),
  #   labels=c("< -60","-30","0","30","> 60"),
  #   oob=scales::squish)+
  coord_sf(expand=F,
           crs=4326)+
  scale_x_continuous(breaks=c(135,140,145,150))+
  labs(x=NULL,
       y=NULL,
       fill='Sept-Nov NDVI Relative Anomaly (%)      ')+
  facet_wrap(~year,drop = T,nrow=1)+
  theme_linedraw()+
  theme_linedraw()+
  theme(panel.grid=element_blank(), 
        legend.position = 'bottom', 
        legend.key.width = unit(1.75,'cm'), 
        legend.key.height = unit(0.15,'cm'), 
        strip.background = element_blank(),
        strip.text = element_text(color='black',face='bold',size=13))
gc()
ggsave(p_out,
  filename="figures/figure_relative-median-ndvi-anom_SON_2016-2020.png", 
       width=22*(5/4),
       height=8,
       units='cm',
  device=grDevices::png,
       dpi=350)
