##Load all required packages
library(ggplot2) #For plotting

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the explicit source model
#k=1; p=0
esm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    p1 <- aa1*d1; p2 <- aa2*d2
    dn <- d0*(U-n)
    dl1 <- ((d1-p1)/k1)*((k1-1)*U+n)+p1*U-d1*l1
    dl2 <- ((d2-p2)/k2)*((k2-1)*U+l1)+p2*U-d2*l2
    return(list(c(dn,dl1,dl2)))
  })
}

#Function defining a waiting box model
wb <- function(t,state,parms){
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    tlag1 <- t-Delta1
    tlag2 <- tlag1-Delta2
    lags1 <- ifelse(tlag1<0,0,lagvalue(tlag1,1))
    lags2 <- ifelse(tlag2<0,0,lagvalue(tlag2,1))
    dn <- d0*(U-n)
    dl1 <- (n-lags1)/Delta1
    dl2 <- (lags1-lags2)/Delta2
    return(list(c(dn,dl1,dl2)))
  })
}

##Main block
par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,2))
tau <- 4; k1 <- 1; k2 <- 1
s <- c(n=0,l1=0,l2=0); p <- c(d0=0.1,aa1=0,d1=0.2,aa2=0,d2=0.25)
dnsm <- run(20,0.01,hmax=0.1,arrest=tau,odes=esm,table=T,
            timeplot=F)

s <- c(n=0,l1=0,l2=0); p <- c(d0=0.1,Delta1=5,Delta2=4)
dwb <- run(20,0.01,hmax=0.1,arrest=tau,odes=wb,delay=T,table=T,
           timeplot=F)

data <- data.frame(time=c(rep(dnsm$time,3),rep(dwb$time,3)),
                   l=c(dnsm$n,dnsm$l1,dnsm$l2,dwb$n,dwb$l1,dwb$l2),
                   var=c(rep(c('1st','2nd','3rd'),each=length(dnsm$time)),
                         rep(c('1st','2nd','3rd'),each=length(dwb$time))),
                   model=c(rep('(a)',length(dnsm$time)*3),rep('(b)',length(dwb$time)*3)))

#ggplot with facets
print(ggplot(data,aes(x=time,y=l,group=var,col=var))+
        geom_line(size=1)+facet_wrap(~model)+theme_test()+
        theme(plot.title=element_text(hjust=0.5),
              legend.title=element_blank(),
              legend.text.align=0.3,strip.background=element_blank())+
        scale_color_manual(name='',values=c("orange","darkred","dodgerblue"),
                           labels=c("1st","2nd","3rd"))+
        xlab('Time (in scaled time units)')+ylab('Fraction labelled'))

# ggsave('FigK.pdf',width=8,height=5)
# ggsave('FigK.svg',width=8,height=5)
# ggsave('FigK.png',width=8,height=5)
