#function for printing levelplot of cross-section velocities

print_crossection = function(my_data, interpolace, oko){
  grid_size = oko
  mesh = sit_sample(my_data, grid_size)
  SIT_hl = HlouSvislic_sit(mesh$POLYGONCS, mesh$SIT)
  if (interpolace== 1) {JOHA=  arith_mean(my_data, SIT_hl)
      naz = 'The arithmetical mean interpolation' }
  if (interpolace == 2) {JOHA=  inv_dis_w(my_data, SIT_hl)
      naz = 'The inverse distance weighted interpolation'}
  if (interpolace == 3) {JOHA=  ord_kriging(my_data, SIT_hl)
  naz = 'The ordinary kriging interpolation'}
  JOHA = data.frame(JOHA)

  library(latticeExtra)
  rampa <- colorRampPalette(c('red', 'blue', 'cyan'))(20)
  levelplot(JOHA$vysl~ JOHA$stan * JOHA$hlouBodu, 	panel = panel.levelplot.points, type = c("p"),cex = 1.5,
				col = 0,
				col.regions=rampa,
				prepanel = prepanel.default.xyplot,
				colorkey = TRUE,
				main = naz, ylab = "Depth [m]", xlab = "Width [m]"
    )
}
