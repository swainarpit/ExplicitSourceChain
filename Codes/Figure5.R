library(ggplot2)

availMCMV <- function(t,tau){
  if(t<delta){
    U <- 0
  }else{
    if(t<=tau){
      U <- f*(1-exp(-del*t))+b*exp(-del*t)
    }else{
      U <- (f*(1-exp(-del*tau))+b*exp(-del*tau))*exp(-del*(t-tau))
    }
  }
  U <- amp*U
  return(U)
}

kh <- function(t, state, parms){
  with(as.list(c(state, parms)),{
    U <- availMCMV(t,tau)
    dl1 <- d1*(U-l1)
    dl2 <- d2*(U-l2)
    return(list(c(dl1,dl2)))
  })
}

tweakkh <- function(parms, nsol){
  with(as.list(parms), {
    nsol$l <- a*nsol$l1+(1-a)*nsol$l2
    return(nsol)
  })
}

nsm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- availMCMV(t,tau)
    p <- aa*d
    dn <- d1*(U-n)
    dl <- ((d-p)/k)*((k-1)*U+n)+p*U-d*l
    return(list(c(dn,dl)))
  })
}

data <- read.delim("148/148_1wk.txt",sep="\t")
data1 <- data.frame(data$Time,data$M8)
data1[data1<0] <- 0; names(data1) <- c("time","l")

data <- read.delim("148/148_4wk.txt",sep="\t")
data2 <- data.frame(data$Time,data$M8)
data2[data2<0] <- 0; names(data2) <- c("time","l")

data <- read.delim("148/148_8wk.txt",sep="\t")
data3 <- data.frame(data$Time,data$M8)
data3[data3<0] <- 0; names(data3) <- c("time","l")

datadf <- data.frame(time=c(data1$time,data2$time,data3$time),
          l=c(data1$l,data2$l,data3$l),
          set=c(rep('(a) 1 week',length(data1$time)),
                rep('(b) 4 weeks',length(data2$time)),
                rep('(c) 8 weeks',length(data3$time))))

tend <- max(c(data1$time,data2$time,data3$time))
delta <- 0; f <- 0.024; del <- 0.261; b <- 0.015; amp <- 3.093
tstep <- 0.1; hmax <- 0.1; tauall <- c(7,28,56)
vln <- data.frame(tau=tauall,
                  set=c('(a) 1 week','(b) 4 weeks','(c) 8 weeks'))

par(mar=c(2.6,2.6,1.6,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,3))
print("KH")
s <- c(l1=0,l2=0); p <- c(d1=0.1,d2=0.1,a=0.1,tau=1)
free <- names(p)[1:3]; differ <- c(); fixed <- list(tau=tauall)
ff <- fit(list(data1,data2,data3),tstep=tstep,hmax=hmax,
          arrest=tauall,odes=kh,free=free,fixed=fixed,
          differ=differ,lower=0,upper=c(3,3,1),method="Pseudo",
          tweak="nsol<-tweakkh(parms,nsol)")
p[free] <- ff$par
p['tau'] <- tauall[1]
dataff1 <- run(tend,tstep=tstep,arrest="tau",odes=kh,table=T,
               tweak="nsol<-tweakkh(parms,nsol)",timeplot=F)
p['tau'] <- tauall[2]
dataff2 <- run(tend,tstep=tstep,arrest="tau",odes=kh,table=T,
               tweak="nsol<-tweakkh(parms,nsol)",timeplot=F)
p['tau'] <- tauall[3]
dataff3 <- run(tend,tstep=tstep,arrest="tau",odes=kh,table=T,
               tweak="nsol<-tweakkh(parms,nsol)",timeplot=F)

print("NSM")
s <- c(n=0,l=0); p <- c(d1=0.1,aa=0.1,d=0.1,tau=7)
free <- names(p)[1:3]; differ <- c(); fixed <- list(tau=tauall)
k <- 1
e1 <- fit(list(data1,data2,data3),tstep=tstep,hmax=hmax,
          arrest=tauall,odes=nsm,free=free,fixed=fixed,
          method="Pseudo",lower=0,upper=1)
p[free] <- e1$par
p['tau'] <- tauall[1]
datae11 <- run(tend,tstep=tstep,arrest="tau",odes=nsm,table=T,
               timeplot=F)
p['tau'] <- tauall[2]
datae12 <- run(tend,tstep=tstep,arrest="tau",odes=nsm,table=T,
               timeplot=F)
p['tau'] <- tauall[3]
datae13 <- run(tend,tstep=tstep,arrest="tau",odes=nsm,table=T,
               timeplot=F)

fitdf <- data.frame(time=c(dataff1$time,dataff2$time,dataff3$time,
          datae11$time,datae12$time,datae13$time),
          l=c(dataff1$l,dataff2$l,dataff3$l,datae11$l,datae12$l,
          datae13$l),
          set=rep(c(rep('(a) 1 week',length(dataff1$time)),
                    rep('(b) 4 weeks',length(dataff2$time)),
                    rep('(c) 8 weeks',length(dataff3$time))),2),
          Model=c(rep('KH',length(c(dataff1$time,dataff2$time,
          dataff3$time))),rep('ES',length(c(datae11$time,
          datae12$time,datae13$time)))))

print(ggplot(datadf,aes(x=time,y=l,group~set))+geom_point()+
        facet_wrap(~set)+ylim(0,0.05)+xlab('Time (in days)')+
        ylab('Fraction labelled')+theme_test()+
        geom_line(data=fitdf,aes(col=Model,linetype=Model,size=Model))+
        geom_vline(data=vln,aes(xintercept=tau),size=0.2,
                   linetype='dashed')+
        scale_color_manual(name='',values=c("blue","red"),
                           labels=c('KH','ES'))+
        scale_linetype_manual(name='',values=c("solid","dashed"),
                              labels=c('KH','ES'))+
        scale_size_manual(name='',values=c(1,1.25),
                          labels=c('KH','ES'))+
        theme(plot.title=element_text(hjust=0.5),
        legend.text.align=0.3,strip.background=element_blank()))

ggsave('Figure6.pdf',width=8,height=5)
ggsave('Figure6.svg',width=8,height=5)
ggsave('Figure6.png',width=8,height=5)
