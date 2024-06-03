library(ggplot2)

nsm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

pds <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl <- pp*U-dd*l
    return(list(c(dl)))
  })
}

data <- function(p){
  s <- c(l1=0,l2=0)
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=nsm,table=T,
               parms=p,state=s,timeplot=F)
  
  dataf <- data.frame(time=data1$time,l=data1$l2)
  ss <- c(l=0); ps <- c(pp=0.1,dd=0.1)
  fit1 <- fit(dataf,tstp,hmax=tstp,arrest=tau,odes=pds,
              free=names(ps),parms=ps,state=ss,timeplot=F,
              method="Pseudo",lower=0,upper=1)
  
  data2 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=pds,table=T,
               parms=fit1$par,state=ss,timeplot=F)
  
  #Making data
  tempD <- data.frame(dataf$time,dataf$l,mn)
  names(tempD) <- c('time','var','panel')
  dataD <<- rbind(dataD,tempD)
  
  tempF <- data.frame(data2$time,data2$l,mn)
  names(tempF) <- c('time','var','panel')
  dataF <<- rbind(dataF,tempF)
  
  return()
}

par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1))
tstp <- 0.1; tau <- 1; tmax <- 5
dataD <- data.frame(); dataF <- data.frame()

#Fast precursors
#k=1
k <- 1; p <- c(d1=0.2,p=0.25,d=0.5); mn <- '1'; data(p)
p['p'] <- 0; mn <- '2'; data(p)

#k=2
k <- 2; p['p'] <- 0.25; mn <- '3'; data(p)
p['p'] <- 0; mn <- '4'; data(p)

#Panel names
#panells <- c(expression('(a) '*p[2]*' = 0.5'*d[2]*', k = 1'),
#             expression('(b) '*p[2]*' = 0, k = 1'),
#             expression('(c) '*p[2]*' = 0.5'*d[2]*', k = 2'),
#             expression('(d) '*p[2]*' = 0, k = 2'))
panells <- c('(a)','(b)','(c)','(d)')
dataF$panel <- factor(dataF$panel,labels=panells)
dataD$panel <- factor(dataD$panel,labels=panells)

#Plot
print(ggplot(dataD,aes(x=time,y=var))+geom_point()+
        geom_line(data=dataF,aes(x=time,y=var),col=2,size=1)+theme_test()+
        facet_wrap(~panel,labeller=label_parsed)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        geom_vline(xintercept=1,linetype='dotted')+
        xlab('Time (in scaled time units)')+ylab('Fraction labelled'))

ggsave('FigureS3.pdf',width=8,height=5)
ggsave('FigureS3.svg',width=8,height=5)
ggsave('FigureS3.png',width=8,height=5)
