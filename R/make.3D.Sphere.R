make.3D.Sphere <-
function(position.matrix, center, radius, alpha = 0.5){
  #scatterplot3d(x = position.matrix[,4], y = position.matrix[,5], z = position.matrix[,6], type = 'o', color = 'blue',
  #             cex.symbols=c(1,radius))
  
  rgl.open()
  rgl.bg(color = c("white"))
  plot3d(x=position.matrix[,4], y=position.matrix[,5], z=position.matrix[,6], type = 's', size =0.5, add= FALSE, xlab = '', ylab = '', zlab = '')
  plot3d(x=position.matrix[,4], y=position.matrix[,5], z=position.matrix[,6], type = 'l', add= TRUE, col = 'blue')
  decorate3d(main = "Protein Structure", axes= TRUE, xlab = "x-axis", ylab = "y-axis", zlab = "z-axis")
  index <- which(position.matrix[,2]==center)
  spheres3d(x = position.matrix[index,4], y= position.matrix[index,5], z = position.matrix[index,6], radius = radius, color = 'red',
            alpha= alpha)
  #text3d(position.matrix[,4], position.matrix[,5], position.matrix[,6], text = position.matrix[,2], add = TRUE)
  
}
