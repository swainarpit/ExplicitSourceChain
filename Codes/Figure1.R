library(ggplot2)

nsm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

data <- function(p){
  s <- c(l1=0,l2=0)
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=nsm,table=T,
               parms=p,state=s,timeplot=F)
  
  hgh1 <- data1$l2[data1$time==tau]
  upp10 <- (p['d']*(k-1)+p['p'])/k
  upp11 <- (p['d']*(p['d1']+k-1)+p['p']*(1-p['d1']))/(k*(1+p['d']))
  dwn1 <- p['d']-(p['d']-p['p'])*p['d1']/(k*upp11)
  
  #Making data
  ltim <- length(data1$time); tim <- rep(data1$time,2)
  varC <- c(data1$l1,data1$l2); varnC <- c('1l1','1l2')
  vxn <- c(0,1,1,0); vxx <- c(ts,1-ts,1+ts,ts)
  vyn <- c(0,hgh1,hgh1,0)
  vyx <- c(upp10*ts,hgh1-upp11*ts,hgh1-hgh1*dwn1*ts,p['d']*ts)
  varnS <- c('3p*(0)','4p*(1)','5d*(1)','6d2')
  
  #Total data
  tempC <- data.frame(tim,varC,rep(varnC,each=ltim),rep(mn,2*ltim))
  names(tempC) <- c('time','var','varn','panel')
  tempS <- data.frame(vxn,vxx,vyn,vyx,varnS,rep(mn,4))
  names(tempS) <- c('xmin','xmax','ymin','ymax','varn','panel')
  
  dataC <<- rbind(dataC,tempC)
  dataS <<- rbind(dataS,tempS)
  
  print(c(upp10,upp11,dwn1))
  return()
}

par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1))
ts <- 0.4; tstp <- 0.01; tau <- 1; tmax <- 5
dataC <- data.frame(); dataS <- data.frame()

#Fast precursors
#k=1
k <- 1; p <- c(d1=0.5,p=0.1,d=0.2); mn <- '1'; data(p)
p['p'] <- 0; mn <- '2'; data(p)

#k=2
k <- 2; p['p'] <- 0.1; mn <- '3'; data(p)
p['p'] <- 0; mn <- '4'; data(p)

#Panel names
#panells <- c(expression('(a) '*p[2]*' = 0.5'*d[2]*', k = 1'),
#             expression('(b) '*p[2]*' = 0, k = 1'),
#             expression('(c) '*p[2]*' = 0.5'*d[2]*', k = 2'),
#             expression('(d) '*p[2]*' = 0, k = 2'))
panells <- c('(a)','(b)','(c)','(d)')
dataC$panel <- factor(dataC$panel,labels=panells)
dataS$panel <- factor(dataS$panel,labels=panells)

#Legend names
labells <- c(expression(l[1]),expression(l[2]),
             expression(p*'*'*(0)),expression(p*'*'*(1)),
             expression(d*'*'*(1+epsilon)*l[2]),expression(d[2]))

#Plot
print(ggplot(dataC,aes(x=time,y=var,group=varn,
                       linetype=varn,size=varn,col=varn))+
        geom_line()+geom_segment(data=dataS,aes(x=xmin,xend=xmax,
        y=ymin,yend=ymax,col=varn),arrow=arrow(length=unit(0.2,'cm')),
        show.legend=F)+theme_test()+
        facet_wrap(~panel,labeller=label_parsed)+
        theme(plot.title=element_text(hjust=0.5),legend.text.align=0.3,
              strip.background=element_blank())+
        scale_color_manual(name='',values=c(rep("black",2),"red",
              "orange","green","blue"),labels=labells)+
        scale_linetype_manual(name='',values=c("dashed","solid",
              rep("solid",4)),labels=labells)+
        scale_size_manual(name='',values=c(0.5,1,rep(1,4)),labels=labells)+
        geom_vline(xintercept=1,linetype='dotted')+
        xlab('Time (in scaled time units)')+ylab('Fraction labelled'))

ggsave('Figure1.pdf',width=8,height=5)
ggsave('Figure1.svg',width=8,height=5)
ggsave('Figure1.png',width=8,height=5)
