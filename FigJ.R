##Load all required packages
library(ggplot2) #For plotting

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the explicit source model
esm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

#Function defining the implicit kinetic heterogeneity model
pd <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U  <- ifelse(t<tau,1,0)
    dl <- ps*U-ds*l
    return(list(c(dl)))
  })
}

#Function to generate the data and fit it
fn <- function(p){
  s <- c(l1=0,l2=0)
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=esm,table=T,
               parms=p,state=s,timeplot=F)
  data1[which(data1$time==tau),] <- NA
  data <- na.omit(data.frame(time=data1$time,l=data1$l2))
  
  tsstp <- 0.1
  ss <- c(l=0); pp <- c(ps=0.1,ds=0.1)
  fit1 <- fit(data,tstep=tsstp,hmax=tsstp,arrest=tau,
              state=ss,parms=pp,free=names(pp),weight="mean",
              odes=pd,lower=0,upper=5,method="Pseudo",timeplot=F)
  
  data2 <- run(tmax,tsstp,hmax=tsstp,arrest=tau,odes=pd,table=T,
               parms=fit1$par,state=ss,timeplot=F)
  
  #Saving data to the global dataframe
  tempD <- data.frame(data$time,data$l,x)
  names(tempD) <- c('time','var','panel')
  dataD <<- rbind(dataD,tempD)
  
  tempF <- data.frame(data2$time,data2$l,x)
  names(tempF) <- c('time','var','panel')
  dataF <<- rbind(dataF,tempF)
  
  d1_v <<- c(d1_v,p[1]); p2_v <<- c(p2_v,p[2]); d2_v <<- c(d2_v,p[3])
  ps_v <<- c(ps_v,round(fit1$par[1],2)); ds_v <<- c(ds_v,round(fit1$par[2],2))
  ssr_v <<- c(ssr_v,round(fit1$ssr,4))
  
  return()
}

##Main block
#load('Fig4.RData')
d1_v <- c(); p2_v <- c(); d2_v <- c(); ps_v <- c(); ds_v <- c(); ssr_v <- c()
dataD <- data.frame(); dataF <- data.frame()

tau <- 1; k <- 1; tmax <- 5
tstp <- 0.1
x <- 1; temp <- which(ssrv<quantile(ssrv,0.8)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))
x <- 2; temp <- which(ssrv<quantile(ssrv,0.5)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))
x <- 3; temp <- which(ssrv<quantile(ssrv,0.25)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))

tstp <- 0.4
x <- 4; temp <- which(ssrv<quantile(ssrv,0.8)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))
x <- 5; temp <- which(ssrv<quantile(ssrv,0.5)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))
x <- 6; temp <- which(ssrv<quantile(ssrv,0.25)); y <- which(ssrv==max(ssrv[temp]))
fn(round(c(d1=d1v[y],p=p2v[y],d=d2v[y]),5))

#Panel names
panells <- c('(a)','(b)','(c)','(d)','(e)','(f)')
dataF$panel <- factor(dataF$panel,labels=panells)
dataD$panel <- factor(dataD$panel,labels=panells)

#ggplot with facets
d1_v <- unname(d1_v); p2_v <- unname(p2_v); d2_v <- unname(d2_v)

dat_txt1 <- data.frame(label=as.character(as.expression(sapply(seq(6), function(i) 
  bquote(d[1]~" = "~.(d1_v[i])~", "~p[2]~" = "~.(p2_v[i])~", "~d[2]~" = "~.(d2_v[i]))))),
  panel=panells,x=3,y=0.8)
dat_txt2 <- data.frame(label=paste('p* = ',ps_v,', d* = ',ds_v,sep=''),
                       panel=panells,x=3,y=0.6)
dat_txt3 <- data.frame(label=paste('SSR: ',ssr_v,sep=''),panel=panells,x=3,y=0.4)

print(ggplot(dataD,aes(x=time,y=var))+geom_point()+
        geom_line(data=dataF,aes(x=time,y=var),col=2,size=1)+theme_test()+
        facet_wrap(~panel,labeller=label_parsed)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        geom_vline(xintercept=1,linetype='dotted')+
        xlab('Time (in scaled time units)')+ylab('Fraction labelled')+
        geom_text(data=dat_txt1,mapping=aes(x=x,y=y,label=label),size=3,parse=T)+
        geom_text(data=dat_txt2,mapping=aes(x=x,y=y,label=label),size=3)+
        geom_text(data=dat_txt3,mapping=aes(x=x,y=y,label=label),size=3))

# ggsave('FigJ.pdf',width=8,height=5)
# ggsave('FigJ.svg',width=8,height=5)
# ggsave('FigJ.png',width=8,height=5)
