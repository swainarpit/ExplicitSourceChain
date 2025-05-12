##Load all required packages
library(pdftools)
library(ggplot2)

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the deuterium labelling availability in plasma
plsma <- function(t,state,parms){
  parms <- 10^parms
  with(as.list(c(state,parms)),{
    if(t<tau){
      P <- f*(1-exp(-del*t))+beta*exp(-del*t)
    }else{
      P <- (f*(1-exp(-del*tau))+beta*exp(-del*tau))*exp(-del*(t-tau))
    }
    return(list(c(P)))
  })
}

#Function defining a kinetically heterogeneous population
chain1 <- function(t,state,parms){
  parms <- 10^parms
  with(as.list(c(state,parms)),{
    if(t<tau){
      P <- f*(1-exp(-del*t))+beta*exp(-del*t)
    }else{
      P <- (f*(1-exp(-del*tau))+beta*exp(-del*tau))*exp(-del*(t-tau))
    }
    P <- 3.5*P
    dl11 <- d11*(P-l11)
    dl12 <- d12*(P-l12)
    dl1 <- (alp*d11+(1-alp)*d12)*P - (alp*d11*l11+(1-alp)*d12*l12)
    return(list(c(dl11,dl12,dl1)))
  })
}

#Function defining a chain of 2 populations, where the precursors are kinetically heterogeneous
chain2 <- function(t,state,parms){
  parms <- 10^parms
  with(as.list(c(state,parms)),{
    if(t<tau){
      P <- f*(1-exp(-del*t))+beta*exp(-del*t)
    }else{
      P <- (f*(1-exp(-del*tau))+beta*exp(-del*tau))*exp(-del*(t-tau))
    }
    P <- 3.5*P
    p2 <- alp2*d2
    dl11 <- d11*(P-l11)
    dl12 <- d12*(P-l12)
    dl1 <- (alp*d11+(1-alp)*d12)*P - (alp*d11*l11+(1-alp)*d12*l12)
    dl2 <- (d2*(l1+(k-1)*P-k*l2)+p2*(P-l1))/k
    return(list(c(dl11,dl12,dl1,dl2)))
  })
}

#Function defining a chain of 3 populations, where the precursors are kinetically heterogeneous
chain3 <- function(t,state,parms){
  parms <- 10^parms
  with(as.list(c(state,parms)),{
    if(t<tau){
      P <- f*(1-exp(-del*t))+beta*exp(-del*t)
    }else{
      P <- (f*(1-exp(-del*tau))+beta*exp(-del*tau))*exp(-del*(t-tau))
    }
    P <- 3.5*P
    p1 <- alp1*d1; p2 <- alp2*d2
    dl01 <- d01*(P-l01)
    dl02 <- d02*(P-l02)
    dl0 <- (alp*d01+(1-alp)*d02)*P - (alp*d01*l01+(1-alp)*d02*l02)
    dl1 <- (d1*(l0+(k-1)*P-k*l1)+p1*(P-l0))/k
    dl2 <- (d2*(l1+(k-1)*P-k*l2)+p2*(P-l1))/k
    return(list(c(dl01,dl02,dl0,dl1,dl2)))
  })
}

#Function to transform data before fitting
#This helps to reduce the bias against data points with low absolute values
func <- function(x){
  return(asin(sqrt(x)))
}

#Function for fitting all models
fitt <- function(datax,x){
  #Fit plasma data
  dataP <- datax[datax$subset=='Plasma',c(1,2)]
  names(dataP) <- c('time','P')
  
  s <- c(P=0); p <- c(tau=49,f=0.015,del=0.11,beta=0.01)
  free <- names(p)[-1]; lower <- -5; upper <- c(1,1,5)
  p <- log10(p); upper <- log10(upper)
  
  print('Plasma')
  fP <- fit(dataP,arrest="tau",free=free,hmax=0.01,solution=T,
            odes=plsma,state=s,parms=p,ymax=0.05,timeplot=F,
            lower=lower,upper=upper,#method="Pseudo",
            xlab="Time (in days)",ylab="Fraction labelled")
  xP <- round(10^fP$par,5)
  
  #Fit all models
  k <<- 1
  
  #Fit the kinetic heterogeneity model
  dataC4p <- datax[datax$subset=='CD4+CD57+',]
  dataC4 <- data.frame(time=dataC4p$time,l1=dataC4p$APE)
  ggData <<- rbind(ggData,data.frame(time=dataC4$time,l=dataC4$l1,
                                     sett='CD57+',modd=modnames[1]))
  
  s <- c(l11=0,l12=0,l1=0); p <- c(tau=49,10^fP$par,d11=0.1,d12=0.1,alp=0.1)
  free <- names(p)[-c(1,2,3,4)]
  lower <- -5; upper <- 1
  p <- log10(p); upper <- log10(upper)
  
  fC41p <- fit(dataC4,arrest="tau",free=free,hmax=0.01,fun=func,
               odes=chain1,ymax=0.05,state=s,parms=p,timeplot=F,
               lower=lower,upper=upper,#method="Pseudo",
               xlab="Time (in days)",ylab="Fraction labelled")
  xC41p <- round(10^fC41p$par,5); p[free] <- fC41p$par
  dC41p <- run(130,arrest="tau",hmax=0.01,table=T,timeplot=F,
               odes=chain1,state=s,parms=p)
  ggSim <<- rbind(ggSim,data.frame(time=dC41p$time,l=dC41p$l1,
                                   sett='CD57+',modd=modnames[1]))
  
  #Fit the 2-population chain model
  dataC4n <- datax[datax$subset=='CD4+CD57-',]
  dataC4p <- datax[datax$subset=='CD4+CD57+',]
  tempdataC4n <- dataC4n; tempdataC4n$APE <- NA
  tempdataC4p <- dataC4p; tempdataC4p$APE <- NA
  dataC4 <- data.frame(time=c(dataC4n$time,dataC4p$time),
                       l1=c(dataC4n$APE,tempdataC4p$APE),
                       l2=c(tempdataC4n$APE,dataC4p$APE))
  ggData <<- rbind(ggData,data.frame(time=dataC4$time,l=dataC4$l1,
                                     sett='CD57-',modd=modnames[2]))
  ggData <<- rbind(ggData,data.frame(time=dataC4$time,l=dataC4$l2,
                                     sett='CD57+',modd=modnames[2]))
  
  s <- c(l11=0,l12=0,l1=0,l2=0); p <- c(tau=49,10^fP$par,d11=0.1,d12=0.1,alp=0.1,d2=0.1,alp2=0.1)
  free <- names(p)[-c(1,2,3,4)]
  lower <- -5; upper <- 1; p <- log10(p); upper <- log10(upper)
  
  fC42 <- fit(dataC4,arrest="tau",free=free,hmax=0.01,fun=func,
              odes=chain2,ymax=0.05,state=s,parms=p,timeplot=F,
              lower=lower,upper=upper,#method="Pseudo",
              xlab="Time (in days)",ylab="Fraction labelled")
  xC42 <- round(10^fC42$par,5); p[free] <- fC42$par
  dC42 <- run(130,arrest="tau",hmax=0.01,table=T,timeplot=F,
              odes=chain2,state=s,parms=p)
  ggSim <<- rbind(ggSim,data.frame(time=dC42$time,l=dC42$l1,
                                   sett='CD57-',modd=modnames[2]))
  ggSim <<- rbind(ggSim,data.frame(time=dC42$time,l=dC42$l2,
                                   sett='CD57+',modd=modnames[2]))
  
  #Fit the 3-population chain model
  ggData <<- rbind(ggData,data.frame(time=dataC4$time,l=dataC4$l1,
                                     sett='CD57-',modd=modnames[3]))
  ggData <<- rbind(ggData,data.frame(time=dataC4$time,l=dataC4$l2,
                                     sett='CD57+',modd=modnames[3]))
  ggData <<- rbind(ggData,data.frame(time=NA,l=NA,
                                     sett='Precursors',modd=modnames[3]))
  
  p <- c(tau=49,10^fP$par,d01=0.1,d02=0.1,alp=0.1,d1=0.1,alp1=0.1,d2=0.1,alp2=0.1)
  s <- c(l01=0,l02=0,l0=0,l1=0,l2=0); free <- names(p)[-c(1,2,3,4)]
  lower <- -5; upper <- 1
  p <- log10(p); upper <- log10(upper)
  
  fC43 <- fit(dataC4,arrest="tau",free=free,hmax=0.01,fun=func,
              state=s,parms=p,odes=chain3,ymax=0.05,timeplot=F,
              lower=lower,upper=upper,#method="Pseudo",
              xlab="Time (in days)",ylab="Fraction labelled")
  xC43 <- round(10^fC43$par,5); p[free] <- fC43$par
  dC43 <- run(130,arrest="tau",hmax=0.01,table=T,timeplot=F,
              odes=chain3,state=s,parms=p)
  ggSim <<- rbind(ggSim,data.frame(time=dC43$time,l=dC43$l1,
                                   sett='CD57-',modd=modnames[3]))
  ggSim <<- rbind(ggSim,data.frame(time=dC43$time,l=dC43$l2,
                                   sett='CD57+',modd=modnames[3]))
  ggSim <<- rbind(ggSim,data.frame(time=dC43$time,l=dC43$l0,
                                   sett='Precursors',modd=modnames[3]))
  
  ggData$sett <<- factor(ggData$sett,levels=c('Precursors','CD57-','CD57+'))
  ggSim$sett <<- factor(ggSim$sett,levels=c('Precursors','CD57-','CD57+'))
  
  return(c(xP,xC41p,xC42,xC43))
}

##Main block
#Read data
data <- read.table('CD57/Becca/AllData.txt',header=T)
data[data<0] <- 0; data$time <- round(data$time)

ggData <- data.frame(); ggSim <- data.frame() #Define dataframe to be used to plot
modnames <- c('(a) Single population',
              '(b) Two populations','(c) Three populations')
datai <- data[data$ID=='DW02',c(1,2,3)]; ests <- fitt(datai,'DW02')
print('Plasma'); print(ests[1:3])
print('a'); print(ests[4:6])
print('b'); print(ests[7:11]); print('c'); print(ests[12:18])

#Make the ggplot
par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1))
print(ggplot(ggData,aes(x=time,y=l,color=sett))+geom_point()+
        facet_wrap(~modd)+theme_test()+xlim(0,130)+ylim(0,0.05)+
        geom_line(data=ggSim,aes(x=time,y=l,color=sett,group=sett,linetype=sett))+
        geom_vline(aes(xintercept=49),size=0.4,linetype='dashed')+
        labs(y='Deuterium enrichment (APE)',x='Time (in days)')+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        scale_linetype_manual(name='Subsets',values=c(2,4,1),
                              labels=c("Precursors","CD57-","CD57+"))+
        scale_color_manual(name='Subsets',values=c("orange","darkred","dodgerblue"),
                           labels=c("Precursors","CD57-","CD57+"),
                           guide=guide_legend(override.aes=list(shape=NA))))

# save.image('Fig6.RData')
# ggsave('Fig6.pdf',width=8,height=5)
# ggsave('Fig6.svg',width=8,height=5)
# ggsave('Fig6.png',width=8,height=5)
