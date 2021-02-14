
#--------------------------------------#
### colScale and grayScale functions ###
#--------------------------------------#
colScale <- function(min, max, nb, epsilon)
{
  if(min*max < 0)
  {
    colorTable <- designer.colors(3*nb-1, c( "blue","white", "red"))
    M <- max(abs(min), abs(max))
    brks<- c(seq(-M, -epsilon, length.out = nb+1)[-(nb+1)],
             seq(-epsilon, epsilon,length.out = nb+1)[-(nb+1)],
             seq( epsilon, M,length.out = nb))
    if(-M < min)
    {
      rk <- sum(brks < min)
      colorTable <- colorTable[rk:(3*nb-1)]
      brks <- brks[rk:(3*nb)]
    }
    else
    {
      rk <- sum(brks > max)
      colorTable <- colorTable[1:(3*nb-1-rk)]
      brks <- brks[1:(3*nb-rk)]
    }
  }
  if(min*max >= 0)
  {
    if(min < 0)
    {
      colorTable <- designer.colors(nb-1, c( "blue","white"))
      brks <- seq(min, 0, length.out = nb)
    }
    else
    {
      colorTable <- designer.colors(nb-1, c( "white","red"))
      brks<- seq(0, max, length.out = nb)
    }

  }
  return(list(col = colorTable, breaks = brks))
}


grayScale <- function(max, nb)
{
  colorTable <- designer.colors(nb-1, c("white", "black"))
  brks<- seq(0, max, length.out = nb)
  return(list(col = colorTable, breaks = brks))
}



#------------------------#
### Plot simu function ###
#------------------------#

plotSimu <- function(z1, w1, z2, w2, z3, w3,
                     u1 = u1,
                     v1 = v1,
                     u2 = u2,
                     v2 = v2,
                     u3 = u3,
                     v3 = v3,
                     .phi = phi,
                     .logOmega2 = logOmega2,
                     .myscale = myscale,
                     .positions = positions
)

{
  par(mfrow=c(2,3))
  #par(mar = c(3,4,2,3))
  par(mar = c(4,4,3,3))###change

  #####sd_Eta
  image.plot(.phi, .logOmega2, z1, breaks=u1$breaks, col=u1$col, axis.args=list(cex.axis=2),
             main = expression(E((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  #####sd_Nu
  image.plot(.phi, .logOmega2,z2,breaks=u2$breaks, col=u2$col ,axis.args=list(cex.axis=2),
             main = expression(E((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)


  #####phi
  image.plot(.phi, .logOmega2,z3, breaks=u3$breaks, col=u3$col, axis.args=list(cex.axis=2),
             main = expression(E (hat(phi) - phi) ),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ), mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)



  image.plot(.phi, .logOmega2,w1,breaks=v1$breaks, col=v1$col,axis.args=list(cex.axis=2),
             main = expression(SD((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  image.plot(.phi, .logOmega2,w2,breaks=v2$breaks, col=v2$col,axis.args=list(cex.axis=2),
             main = expression(SD((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  image.plot(.phi, .logOmega2,w3,breaks=v3$breaks, col=v3$col,axis.args=list(cex.axis=2),
             main = expression(SD (hat(phi) - phi)),
             xlab = expression(paste("      ",  phi ) ),
             ylab = expression( omega^2 ),  mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)
}



#-------------------------#
###      FIGURE 10      ###
#-------------------------#

u1 <- colScale(min = min(z1_2,z1_3), max = max(z1_2,z1_3),
               epsilon = abs(max(z1_2,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1_2,w1_3),  nb = 20)
u2 <- colScale(min = min(z2_2,z2_3), max = max(z2_2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2_2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3_2,z3_3), max = max(z3_2,z3_3), epsilon = abs(max(z3_2,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3_2,w3_3),  nb = 20)

u1 <- colScale(min = min(z1_1,z1_2,z1_3), max = max(z1_1,z1_2,z1_3),
               epsilon = abs(max(z1_1,z1_2,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1_1,w1_2,w1_3),  nb = 20)
u2 <- colScale(min = min(z2_1,z2_2,z2_3), max = max(z2_1,z2_2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2_1,w2_2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3_1,z3_2,z3_3), max = max(z3_1,z3_2,z3_3), epsilon = abs(max(z3_1,z3_2,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3_1,w3_2,w3_3),  nb = 20)
plotSimu(z1_1, w1_1, z2_1, w2_1, z3_1, w3_1, u1, v1, u2, v2, u3, v3)
plotSimu(z1_2, w1_2, z2_2, w2_2, z3_2, w3_2, u1, v1, u2, v2, u3, v3)
plotSimu(z1_3, w1_3, z2_3, w2_3, z3_3, w3_3, u1, v1, u2, v2, u3, v3)


#-------------------------#
###      FIGURE 11      ###
#-------------------------#

u1 <- colScale(min = min(z1,z1_3), max = max(z1,z1_3),
               epsilon = abs(max(z1,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1,w1_3),  nb = 20)
u2 <- colScale(min = min(z2,z2_3), max = max(z2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3,z3_3), max = max(z3,z3_3), epsilon = abs(max(z3,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3,w3_3),  nb = 20)

plotSimu(z1_3, w1_3, z2_3, w2_3, z3_3, w3_3, u1, v1, u2, v2, u3, v3)
plotSimu(z1, w1, z2, w2, z3, w3, u1, v1, u2, v2, u3, v3)
