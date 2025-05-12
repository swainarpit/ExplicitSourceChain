##Load all required packages
library(ggplot2) #For plotting

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the differential equation of explicit source model
esm <- function(t, state, parms) { #Use current time, variable and parameter values
  with(as.list(c(state,parms)), { #To be able to use the variable and parameter names
    U <- ifelse(t<tau,1,0) #Deuterium availability function
    dl1 <- d1*(U-l1) #Precursors
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2 #POI
    return(list(c(dl1,dl2)))
  })
}

#Function to generate the time trajectory and slopes
data <- function(p){
  s <- c(l1=0,l2=0) #State vector: precursors and POI are both unlabelled at the start
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=esm,table=T,
               parms=p,state=s,timeplot=F) #Simulate the explicit source model
  
  #Calculate the slopes from the simulated trajectory
  hgh1 <- data1$l2[data1$time==tau] #POI at the end of the labelling period
  upp10 <- (p['d']*(k-1)+p['p'])/k #Calculated p*(0)
  upp11 <- (p['d']*(p['d1']+k-1)+p['p']*(1-p['d1']))/(k*(1+p['d'])) #Calculated p*(1)
  dwn1 <- p['d']-(p['d']-p['p'])*p['d1']/(k*upp11) #Calculated d*(1+e)
  
  #Saving data in a certain way to make ggplot
  ltim <- length(data1$time); tim <- rep(data1$time,2)
  varC <- c(data1$l1,data1$l2); varnC <- c('1l1','1l2')
  vxn <- c(0,1,1,0); vyn <- c(0,hgh1,hgh1,0)
  if(mn=='2'){
    vxx <- c(ts*0.7,1-ts*0.5,1+ts*0.5,ts*0.1)
    vyx <- c(upp10*ts*0.7,
             hgh1-upp11*ts*0.5,hgh1-hgh1*dwn1*ts*0.5,p['d']*ts*0.1)
  }else{
    vxx <- c(ts*0.2,1-ts*0.5,1+ts*0.1,ts*0.1)
    vyx <- c(upp10*ts*0.2,
             hgh1-upp11*ts*0.5,hgh1-hgh1*dwn1*ts*0.1,p['d']*ts*0.1)
  }
  varnS <- c('3p*(0)','4p*(1)','5d*(1)','6d2')
  
  #Create a dataframe
  tempC <- data.frame(tim,varC,rep(varnC,each=ltim),rep(mn,2*ltim))
  names(tempC) <- c('time','var','varn','panel')
  tempS <- data.frame(vxn,vxx,vyn,vyx,varnS,rep(mn,4))
  names(tempS) <- c('xmin','xmax','ymin','ymax','varn','panel')
  
  #Save dataframe to the global dataframe
  dataC <<- rbind(dataC,tempC)
  dataS <<- rbind(dataS,tempS)
  
  print(c(upp10,upp11,dwn1))
  return()
}

##Main block
par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1)) #Define the graphics of the plot
ts <- 0.4; tstp <- 0.01; tau <- 1; tmax <- 5 #Define the timeline of the simulation
dataC <- data.frame(); dataS <- data.frame() #Define a global dataframe

##Simulate all 4 panels of the figure
#Very fast POI
#k=1
k <- 1; p <- c(d1=0.5,p=5,d=10); mn <- '1'; data(p)
p['p'] <- 0; mn <- '2'; data(p)

#k=2
k <- 2; p['p'] <- 5; mn <- '3'; data(p)
p['p'] <- 0; mn <- '4'; data(p)

#Panel names
panells <- c('(a)','(b)','(c)','(d)')
dataC$panel <- factor(dataC$panel,labels=panells)
dataS$panel <- factor(dataS$panel,labels=panells)

#Legend names
labells <- c("Precursors","POI",
             expression(p*'*'*(0)),expression(p*'*'*(1)),
             expression(d*'*'*(1+epsilon)*l[2]),expression(d[2]))

#Offset predicted slopes
dataS[,3:4] <- dataS[,3:4] + 0.01*rep(c(-1,-1,1,1),4)

#ggplot with facets
print(ggplot(dataC,aes(x=time,y=var,group=varn,
                       linetype=varn,size=varn,col=varn))+
        geom_line()+geom_segment(data=dataS,aes(x=xmin,xend=xmax,
        y=ymin,yend=ymax,col=varn),show.legend=F)+theme_test()+
        geom_point(data=dataS,aes(x=xmin,y=ymin,col=varn),size=2,show.legend=F)+
        facet_wrap(~panel,labeller=label_parsed)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        scale_color_manual(name='',values=c(rep("#000000",2),"#CD0BBC",
              "orange","darkred","dodgerblue"),labels=labells)+
        scale_linetype_manual(name='',values=c("dashed","solid",
                        rep("solid",4)),labels=labells)+
        scale_size_manual(name='',values=c(0.5,1,rep(1,4)),labels=labells)+
        geom_vline(xintercept=1,linetype='dotted')+
        xlab('Time (in scaled time units)')+ylab('Fraction labelled'))

# ggsave('FigA.pdf',width=8,height=5)
# ggsave('FigA.svg',width=8,height=5)
# ggsave('FigA.png',width=8,height=5)
