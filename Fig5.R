##Load all required packages
library(latticeExtra)
library(RColorBrewer)
library(Orcs)

#Function to make individual panels
panels <- function(j){
  k <- ifelse(j<3,1,2) #Type of differentiation
  data <- expand.grid(x=x,y=y) #Combinations of x and y values
  data$z <- sapply(seq(nrow(data)), function(i){ #Function to calculate error
    d1 <- data$y[i]; d2 <- data$x[i]
    if(j%%2==1) p2 <- 0*d2
    if(j%%2==0) p2 <- 0.8*d2
    temp <- ifelse(d1==d2,(1-d1)*exp(-d1),(d2*exp(-d2)-d1*exp(-d1))/(d2-d1))
    p_star_1 <- d2*exp(-d2) - (d2-p2)*temp/k
    err <- p_star_1/d2
    return(err)
  })
  
  #Make individual panels
  P <- levelplot(z~x*y,data,#aspect='iso',
                 xlab='True turnover rate of the POI (relative to the labelling period)',
                 ylab='True turnover rate of the precursors (relative to the labelling period)',
                 scales=list(y=list(log=T),x=list(log=T),tck=c(1,0)),
                 xscale.components=xscale.components.log10ticks,
                 yscale.components=yscale.components.log10ticks,
                 col.regions=colorRampPalette(brewer.pal(8,'RdBu'))(25),
                 panel=panel.levelplot.points,cex=1)
  return(P)
}

##Main block
lbx <- 0.01; ubx <- 10; nx <- 100 #Choose number of partitions on x-axis
x <- exp(seq(log(lbx),log(ubx),by=(log(ubx)-log(lbx))/nx)) #Uniform partitions on log-scale
lby <- 0.01; uby <- 10; ny <- 100 #Choose number of partitions on y-axis
y <- exp(seq(log(lby),log(uby),by=(log(uby)-log(lby))/ny)) #Uniform partitions on log-scale

#Define panels
panel_names <- c('(a)','(b)','(c)','(d)')
plots <- vector('list',length=4)
for(i in seq(4)){
  plots[[i]] <- panels(i) #Make individual panels
}

#Plot the figure
print(latticeCombineGrid(plots,layout=c(2,2),between=list(y=1.5,x=1.5)))
