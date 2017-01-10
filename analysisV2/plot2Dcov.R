

### Plotting function plot_2Dcov 

# Plots a genome-wide distribution for two variables in greyscale
### The number of bins is given by nbin
### The 95% prediction ellipse for all data is plotted in black

# Plots up to 4 groups of SNPs in different colors (orange, blue, green, or golden)
### The 95% prediction ellipse for each group is plotted in the respective color

# If the absolute value of association between x and y (PE) 
### is greater than 0.1, the function will colors the quadrant 
### for the expected direction based on the x,y association. 
### This version also draws the 95% CI on the covariance matrix.

plot2Dcov <- function(x,y, xlab, ylab, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), 
                      nbin, 
                      x_sub_orange=NULL, y_sub_orange=NULL, 
                      x_sub_green=NULL, y_sub_green=NULL, 
                      x_sub_blue=NULL,  y_sub_blue=NULL,
                      x_sub_golden=NULL,y_sub_golden=NULL, 
                      PE=0, PElab = c("V1", "V2")){

  data1 <- cbind(x, y)
  data1b <- data1[complete.cases(data1),]
  if(length(xlim)==0){
    xlim_up <- max(x, na.rm=TRUE)+0.15*max(x, na.rm=TRUE)
    xlim_lower <- min(x, na.rm=TRUE)-0.15*max(x, na.rm=TRUE)
  }else{
    xlim_lower <- xlim[1]; xlim_up <- xlim[2]
  }
  if(length(ylim)==0){
    ylim_up <- max(y, na.rm=TRUE)+0.15*max(y, na.rm=TRUE)
    ylim_lower <- min(y, na.rm=TRUE)-0.15*max(y, na.rm=TRUE)
  }else{
    ylim_lower <- ylim[1]; ylim_up <- ylim[2]
  }
        binned <- bin2(data1b,
                 matrix(c(xlim_lower,xlim_up,ylim_lower,ylim_up), 2,2, byrow=TRUE),
                 nbin=c(nbin,nbin))
    binned$nc[binned$nc==0]=NA

    ### Start plotting 
    xlocs <- seq(xlim_lower,xlim_up,length.out = nbin)
    ylocs <- seq(ylim_lower,ylim_up, length.out=nbin)
    image(xlocs,  ylocs,
               z=matrix(0, nrow=length(ylocs), ncol=length(xlocs)),
          col=(rgb(0,0,0,0)),
             xlab=xlab, ylab=ylab, add=FALSE, bty="n")#, 
             #xaxt="n", yaxt="n", bty="n")
    
    if(PE >= 0.1){
      polygon(c(0,1,1,0,0), c(0,0,1,1,0),col=rgb(0,0,0,0.05))
      polygon(c(0,-1,-1,0,0), c(0,0,-1,-1,0), col=rgb(0,0,0,0.05))
      }
    if(PE <= -0.1){
      polygon(c(0,-1,-1,0,0), c(0,0,1,1,0),col=rgb(0,0,0,0.05))
      polygon(c(0,1,1,0,0), c(0,0,-1,-1,0),col=rgb(0,0,0,0.05))
    }
    lines(c(-100,100), c(0,0), col="black") #add x axis
    lines(c(0,0), c(-100,100), col="black") # add y axis
    
    image(seq(xlim_lower,xlim_up,length.out = nbin), 
               seq(ylim_lower,ylim_up, length.out=nbin),
               binned$nc,
             xlab=xlab, ylab=ylab, add=TRUE, 
             col=grey.colors(75, 0.8,0.1))
             #xaxt="l", yaxt="l", bty="n")
    
    ### Function for plotting covariance matrix
    plot_mycov <- function(data1, CI = 0.95, color="black", ...){
      # data1 is a data frame with x, y in columns
      C.ls2 <- cov(data1[complete.cases(data1),])
      m.ls2 <- colMeans(data1[complete.cases(data1),])
      d2.95 <- qchisq(CI, df = 2)
      #d2.95 <- qnorm(0.99)
      lines(ellipsoidPoints(C.ls2, d2.95, loc=m.ls2), lwd=2, col=color, ...)
    }
    
    ### Make cov plot for all data
    plot_mycov(data1)
    
    ### Make function for plotting sub points
    plot_mycov_sub <- function(x_sub, y_sub, outlinecolor, bgcolor, linecolor,  mylty, mypch, ...){
      if(length(x_sub)>0){
      points(x_sub, y_sub, cex=1, col=outlinecolor, bg=bgcolor, pch=mypch)
      plot_mycov(data.frame(x_sub, y_sub), col=linecolor, lty=mylty,...)
      }
    }
    
    ### Make cov plot for x_sub_orange, y_sub_orange
    plot_mycov_sub(x_sub_orange, y_sub_orange, outlinecolor=adjustcolor("darkorange",0.7), bgcolor=adjustcolor("orange",0.7), mypch=24, mylty=2, linecolor="#CC6600")
    
    plot_mycov_sub(x_sub_green, y_sub_green, outlinecolor=adjustcolor("darkgreen",0.7), bgcolor="lightgreen", mypch=21, mylty=3, linecolor="darkgreen")
    
    plot_mycov_sub(x_sub_blue, y_sub_blue, outlinecolor=adjustcolor("blue",0.7), bgcolor=adjustcolor("lightblue",0.7), mypch=22, mylty=4, linecolor="darkblue")
    
    plot_mycov_sub(x_sub_golden, y_sub_golden, linecolor="darkred", bgcolor="gold", mypch=23, mylty=6)

    t <- bquote(~rho ~ "(" ~ .(PElab[1]) ~ "," ~ .(PElab[2]) ~ ") =" ~ .(PE))
    text(xlim_up*0.5,ylim_lower*0.95, t, cex=1.2)
} # end plot_2D.b