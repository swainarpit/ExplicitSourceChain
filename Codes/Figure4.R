library(qpdf)
library(gtools)
library(ggplot2)
library(patchwork)
library(svglite)

nsm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

pd <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U  <- ifelse(t<tau,1,0)
    dl <- ps*U-ds*l
    return(list(c(dl)))
  })
}

iter <- function(p){
  s <- c(l1=0,l2=0)
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=nsm,table=T,
               parms=p,state=s,timeplot=F)
  data <- data.frame(time=data1$time,l=data1$l2)
  
  ss <- c(l=0); pp <- c(ps=0.1,ds=0.1)
  fit1 <- fit(data,tstep=tstp,hmax=tstp,
              state=ss,parms=pp,free=names(pp),weight="mean",
              arrest=tau,odes=pd,lower=0,upper=5,method="Pseudo")
  
  return(c(p,fit1$par,fit1$ssr))
}

tau <- 1; tstp <- 0.1; k <- 2; tmax <- 5
par_all <- expand.grid(list(d1vec=seq(0.5,5,0.5),
                            d2vec=seq(0.5,5,0.5),
                            pd2vec=seq(0,0.9,0.1)))
#set.seed(1); temp <- runif(30)
#par_all <- expand.grid(list(d1vec=temp[1:10]*2,
#                            d2vec=temp[11:20]*2,
#                            pd2vec=temp[21:30]))
crs <- ifelse(.Platform$OS.type=="windows",1,detectCores()-1)
all_ests <- pbmclapply(seq(dim(par_all)[1]), function(i){
  pdf(paste('Fit',i,'.pdf'))
  x <- unlist(par_all[i,]); x[3] <- x[2]*x[3]
  names(x) <- c('d1','d','p'); y <- iter(x)
  dev.off()
  y
},mc.cores=crs,ignore.interactive=T)

pdf_combine(mixedsort(list.files(pattern='.pdf')),output='FigureX2_Fits.pdf')
file.remove(tail(mixedsort(list.files(pattern='.pdf')),-1))

d1v <- sapply(all_ests,'[[',1); d2v <- sapply(all_ests,'[[',2)
p2v <- sapply(all_ests,'[[',3); psv <- sapply(all_ests,'[[',4)
dsv <- sapply(all_ests,'[[',5); ssrv <- sapply(all_ests,'[[',6)
ps0v <- p2v; ps1v <- (d2v*d1v+p2v*(1-d1v))/(1+d2v)
dsev <- d2v-(d2v-p2v)*d1v/ps1v

#pdata <- data.frame(x=c(d1v,d2v,p2v,ps0v,psv,ps0v,ps1v,dsv),
#                    y=c(rep(psv,3),ps1v,dsv,psv,psv,dsev),
#                    nm=rep(c('1. d1 vs p*','2. d2 vs p*','3. p2 vs p*',
#                             '4. p*(0) vs p*(1)','5. p* vs d*','6. p*(0) vs p*',
#                             '7. p*(1) vs p*','8. d* vs d*(1)'),
#                           each=length(psv)),
#                    ssr=rep(ssrv,8))
pdata <- data.frame(x=c(d1v,d2v,p2v,(p2v+(d2v-p2v)*(k-1)/k)),
                    y=c(rep(psv,4)),
                    nm=rep(c('(a)','(b)','(c)','(d)'),
                           each=length(psv)),
                    ssr=rep(ssrv,4))
pdata <- pdata[order(-pdata$ssr),]
#pdata[pdata$ssr>quantile(ssrv,0.25),c(2,4)] <- NA

pdata$coll <- '4'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.25))] <- '3'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.50))] <- '2'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.80))] <- '1'

pdata1 <- pdata[pdata$nm=='(a)',]; pdata2 <- pdata[pdata$nm=='(b)',]
pdata3 <- pdata[pdata$nm=='(c)',]; pdata4 <- pdata[pdata$nm=='(d)',]

plt1 <- ggplot(pdata1,aes(x,y))+geom_point(aes(colour=coll))+theme_test()+
  xlab('d1')+ylab('p*')+ggtitle('(a)')+xlim(0,NA)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,col='red')+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','yellow','orange','red'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt2 <- ggplot(pdata2,aes(x,y))+geom_point(aes(colour=coll))+theme_test()+
  xlab('d2')+ylab('p*')+ggtitle('(b)')+xlim(0,NA)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,col='red')+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','yellow','orange','red'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt3 <- ggplot(pdata3,aes(x,y))+geom_point(aes(colour=coll))+theme_test()+
  xlab('p2')+ylab('p*')+ggtitle('(c)')+xlim(0,NA)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,col='red')+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','yellow','orange','red'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt4 <- ggplot(pdata4,aes(x,y))+geom_point(aes(colour=coll))+theme_test()+
  xlab('p*(0)')+ylab('p*')+ggtitle('(d)')+xlim(0,NA)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,col='red')+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','yellow','orange','red'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt <- plt1+plt2+plt3+plt4+plot_layout(ncol=2,guides='collect') & ylab(NULL)# & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
print(wrap_elements(plt)+labs(tag='The estimated average turnover rate, p*')+
        theme(plot.tag=element_text(size=rel(1),angle=90),plot.tag.position='left'))

#save.image('Figure4.RData')
#ggsave('Figure4.pdf',width=8,height=5)
#ggsave('Figure4.svg',width=8,height=5)
#ggsave('Figure4.png',width=8,height=5)
