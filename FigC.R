##Load all required packages
library(ggplot2) #For plotting

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the explicit source model
esm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    if(t<delta){
      U <- 0
    }else{
      if(t<=tau){
        U <- f*(1-exp(-del*t))+b*exp(-del*t)
      }else{
        U <- (f*(1-exp(-del*tau))+b*exp(-del*tau))*exp(-del*(t-tau))
      }
    }
    U <- U/f
    
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

#Function to generate the time trajectory and slopes
data <- function(p){
  s <- c(l1=0,l2=0) #State vector: precursors and POI are both unlabelled at the start
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=esm,table=T,
               parms=p,state=s,timeplot=F) #Simulate the explicit source model
  Ut <- sapply(data1$time, function(i){
    t <- i
    if(t<delta){
      U <- 0
    }else{
      if(t<=tau){
        U <- f*(1-exp(-del*t))+b*exp(-del*t)
      }else{
        U <- (f*(1-exp(-del*tau))+b*exp(-del*tau))*exp(-del*(t-tau))
      }
    }
    U <- U/f
    return(U)
  })
  
  #Calculate the slopes from the simulated trajectory
  hgh1 <- data1$l2[data1$time==tau] #POI at the end of the labelling period
  upp10 <- (p['d']*(k-1)+p['p'])/k #Calculcated p*(0)
  upp11 <- (p['d']*(p['d1']+k-1)+p['p']*(1-p['d1']))/(k*(1+p['d'])) #Calculcated p*(1)
  dwn1 <- p['d']-(p['d']-p['p'])*p['d1']/(k*upp11) #Calculcated d*(1+e)
  
  #Saving data in a certain way to make ggplot
  ltim <- length(data1$time); tim <- rep(data1$time,2)
  varC <- c(data1$l1,data1$l2); varnC <- c('1l1','1l2')
  vxn <- c(0,1,1,0); vxx <- c(ts,1-ts,1+ts,ts)
  vyn <- c(0,hgh1,hgh1,0)
  vyx <- c(upp10*ts,hgh1-upp11*ts,hgh1-hgh1*dwn1*ts,p['d']*ts)
  varnS <- c('3p*(0)','4p*(1)','5d*(1)','6d2')
  
  #Create a dataframe
  tempC <- data.frame(tim,varC,rep(varnC,each=ltim),rep(mn,2*ltim))
  names(tempC) <- c('time','var','varn','panel')
  tempS <- data.frame(vxn,vxx,vyn,vyx,varnS,rep(mn,4))
  names(tempS) <- c('xmin','xmax','ymin','ymax','varn','panel')
  tempU <- data.frame(data1$time,Ut,'0Ut',rep(mn,ltim))
  names(tempU) <- c('time','var','varn','panel')
  
  #Save dataframe to the global dataframe
  dataC <<- rbind(dataC,tempC)
  dataS <<- rbind(dataS,tempS)
  dataU <<- rbind(dataU,tempU)
  
  print(c(upp10,upp11,dwn1))
  return()
}

##Main block
par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1)) #Define the graphics of the plot
delta <- 0; f <- 0.024; del <- 0.1*28; b <- 0.015; amp <- 3.093 #Define the parameters of the deuterated water labelling
ts <- 0.4; tstp <- 0.01; tau <- 1; tmax <- 5 #Define the timeline of the simulation
dataC <- data.frame(); dataS <- data.frame(); dataU <- data.frame() #Define a global dataframe

##Simulate all 4 panels of the figure
#Fast precursors
#k=1
k <- 1; p <- c(d1=0.2,p=0.25,d=0.5); mn <- '1'; data(p)
p['p'] <- 0; mn <- '2'; data(p)

#k=2
k <- 2; p['p'] <- 0.25; mn <- '3'; data(p)
p['p'] <- 0; mn <- '4'; data(p)

#Panel names
panells <- c('(a)','(b)','(c)','(d)')
dataC$panel <- factor(dataC$panel,labels=panells)
dataS$panel <- factor(dataS$panel,labels=panells)
dataU$panel <- factor(dataU$panel,labels=panells)

#Legend names
labells <- c(expression(D*(t)),"Precursors","POI",
             expression(p*'*'*(0)),expression(p*'*'*(1)),
             expression(d*'*'*(1+epsilon)*l[2]),expression(d[2]))

#Offset predicted slopes
dataS[,3:4] <- dataS[,3:4] + 0.01*rep(c(-1,-1,1,1),4)

#ggplot with facets
print(ggplot(data=dataC,aes(group=varn,
                            linetype=varn,size=varn,col=varn))+
        geom_line(data=dataC,aes(x=time,y=var))+
        geom_line(data=dataU,aes(x=time,y=0.5*var))+
        geom_segment(data=dataS,aes(x=xmin,xend=xmax,
                                    y=ymin,yend=ymax,col=varn),show.legend=F)+
        geom_point(data=dataS,aes(x=xmin,y=ymin,col=varn),size=2,show.legend=F)+
        facet_wrap(~panel,labeller=label_parsed)+theme_test()+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        scale_y_continuous(name="Fraction labelled",sec.axis=sec_axis(~.*2,name='Body Deuterium Concentration'))+
        scale_color_manual(name='',values=c('red',rep("#000000",2),"#CD0BBC",
                                            "orange","darkred","dodgerblue"),labels=labells)+
        scale_linetype_manual(name='',values=c("dotted","dashed","solid",
                                               rep("solid",4)),labels=labells)+
        scale_size_manual(name='',values=c(0.8,0.5,1,rep(1,4)),labels=labells)+
        geom_vline(xintercept=1,colour='grey',linetype='dotted')+
        xlab('Time (in scaled time units)'))

# ggsave('FigC.pdf',width=8,height=5)
# ggsave('FigC.svg',width=8,height=5)
# ggsave('FigC.png',width=8,height=5)
