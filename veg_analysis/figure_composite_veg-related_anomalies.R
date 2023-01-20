# pacman::p_load(tidyverse,data.table,lubridate,stars,scico,viridis,cols4all)
library(magick)

p_vpd <- magick::image_read("figures/figure_spattemp_relative-vpd-anom_2016_2020.png") %>%
  image_annotate(.,text = '(a)',size = 80)
p_lst <- magick::image_read("figures/figure_annual-mean-lst-anom_2016-2020.png") %>%
  image_annotate(.,text = '(d)',size = 80)
p_vod <- magick::image_read("figures/figure_relative-median-VOD-anom_SON_2016-2020.png") %>%
  image_annotate(.,text = '(b)',size = 80)
p_ndvi <- magick::image_read("figures/figure_relative-median-ndvi-anom_SON_2016-2020.png") %>%
  image_annotate(.,text = '(c)',size = 80)


trim_y <- 95
p_out <- image_append(c(
  magick::image_crop(p_vpd,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_vod,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_ndvi,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_lst,geometry_size_percent(width=100,height=trim_y))),
 stack=T)


image_write(p_out,"figures/figure_spattemp_veg-related-anoms_wNDVI_2016_2020.png")




p_vpd <- magick::image_read("figures/figure_spattemp_relative-vpd-anom_2016_2020.png")
p_lst <- magick::image_read("figures/figure_annual-mean-lst-anom_2016-2020.png")
p_vod <- magick::image_read("figures/figure_relative-median-VOD-anom_SON_2016-2020.png")
# p_ndvi <- magick::image_read("figures/figure_relative-median-ndvi-anom_SON_2016-2020.png")
p_lai <- magick::image_read("figures/figure_relative-median-LAI-anom_SON_2016-2020.png")

trim_y <- 95
p_out <- image_append(c(
  magick::image_crop(p_vpd,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_vod,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_lai,geometry_size_percent(width=100,height=trim_y)),
  magick::image_crop(p_lst,geometry_size_percent(width=100,height=trim_y))),
 stack=T)


image_write(p_out,"figures/figure_spattemp_veg-related-anoms_wLAI_2016_2020.png")


