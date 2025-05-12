##Load all required packages
library(qpdf)
library(gtools)
library(ggplot2)
library(patchwork)

##Run the wrapper grind.R; used for data fitting
source('grind.R')

##Define all functions
#Function defining the differential equation of explicit source model
esm <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U <- ifelse(t<tau,1,0)
    dl1 <- d1*(U-l1)
    dl2 <- ((d-p)/k)*((k-1)*U+l1)+p*U-d*l2
    return(list(c(dl1,dl2)))
  })
}

#Function defining the differential equation of the implicit kinetic heterogeneity model
pd <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    U  <- ifelse(t<tau,1,0)
    dl <- ps*U-ds*l
    return(list(c(dl)))
  })
}

#Function to generate data and fit the same
iter <- function(p){
  s <- c(l1=0,l2=0) #Initialize the variables
  
  data1 <- run(tmax,tstp,hmax=tstp,arrest=tau,odes=esm,table=T,
               parms=p,state=s,timeplot=F) #Simulate the explicit source model
  data <- data.frame(time=data1$time,l=data1$l2) #Save the required data
  
  ss <- c(l=0); pp <- c(ps=0.1,ds=0.1) #Initialize the variable and initial parameter guesses
  fit1 <- fit(data,tstep=tstp,hmax=tstp,
              state=ss,parms=pp,free=names(pp),weight="mean",
              arrest=tau,odes=pd,lower=0,upper=5,method="Pseudo") #Fit the implicit kinetic heterogeneity model to simulated data
  
  return(c(p,fit1$par,fit1$ssr))
}

##Main block
tau <- 1; tstp <- 0.1; k <- 1; tmax <- 5 #Set timeline and type of differentiation
par_all <- expand.grid(list(d1vec=seq(0.5,5,0.5),
                            d2vec=seq(0.5,5,0.5),
                            pd2vec=seq(0,0.9,0.1))) #Generate a wide range of parameter combinations
#Run the fits in parallel
crs <- ifelse(.Platform$OS.type=="windows",1,detectCores()-1) #Set the number of cores
all_ests <- pbmclapply(seq(dim(par_all)[1]), function(i){
  pdf(paste('Fit',i,'.pdf')) #Save the fits in individual pdfs
  x <- unlist(par_all[i,]); x[3] <- x[2]*x[3]
  names(x) <- c('d1','d','p'); y <- iter(x)
  dev.off()
  y #Save the output of the 'iter' function in a list
},mc.cores=crs,ignore.interactive=T) #Call the function of interest in parallel

#Lines 62-63 work well if the current directory does not already have pdfs 
pdf_combine(mixedsort(list.files(pattern='.pdf')),output='Fig3_Fits.pdf') #Combine all the pdfs
file.remove(tail(mixedsort(list.files(pattern='.pdf')),-1)) #Delete the individual pdfs

#Extract the required quantities from the fits
d1v <- sapply(all_ests,'[[',1); d2v <- sapply(all_ests,'[[',2)
p2v <- sapply(all_ests,'[[',3); psv <- sapply(all_ests,'[[',4)
dsv <- sapply(all_ests,'[[',5); ssrv <- sapply(all_ests,'[[',6)
ps0v <- p2v; ps1v <- (d2v*d1v+p2v*(1-d1v))/(1+d2v)
dsev <- d2v-(d2v-p2v)*d1v/ps1v

#Create a dataframe that can be used to plot
pdata <- data.frame(x=c(d1v,d2v,p2v,(p2v+(d2v-p2v)*(k-1)/k)),
                    y=c(rep(psv,4)),
                    nm=rep(c('(a)','(b)','(c)','(d)'),
                           each=length(psv)),
                    ssr=rep(ssrv,4))
pdata <- pdata[order(-pdata$ssr),]

pdata$coll <- '4'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.25))] <- '3'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.50))] <- '2'
pdata$coll[which(pdata$ssr>quantile(ssrv,0.80))] <- '1'

pdata1 <- pdata[pdata$nm=='(a)',]; pdata2 <- pdata[pdata$nm=='(b)',]
pdata3 <- pdata[pdata$nm=='(c)',]; pdata4 <- pdata[pdata$nm=='(d)',]

#Make individual ggplot panels
plt1 <- ggplot(pdata1,aes(x,y))+geom_point(aes(colour=coll,shape=coll))+theme_test()+
  xlab('d1')+ylab('p*')+ggtitle('(a)')+xlim(0,5)+ylim(0,5)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,lty=2)+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_shape_manual(name='Quantile(SSR)',
                     values=c(15,16,17,18),
                     labels=c('>80%','50-80%','25-50%','<25%'))+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','dodgerblue','orange','darkred'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt2 <- ggplot(pdata2,aes(x,y))+geom_point(aes(colour=coll,shape=coll))+theme_test()+
  xlab('d2')+ylab('p*')+ggtitle('(b)')+xlim(0,5)+ylim(0,5)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,lty=2)+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_shape_manual(name='Quantile(SSR)',
                     values=c(15,16,17,18),
                     labels=c('>80%','50-80%','25-50%','<25%'))+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','dodgerblue','orange','darkred'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt3 <- ggplot(pdata3,aes(x,y))+geom_point(aes(colour=coll,shape=coll))+theme_test()+
  xlab('p2')+ylab('p*')+ggtitle('(c)')+xlim(0,5)+ylim(0,5)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,lty=2)+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_shape_manual(name='Quantile(SSR)',
                     values=c(15,16,17,18),
                     labels=c('>80%','50-80%','25-50%','<25%'))+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','dodgerblue','orange','darkred'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

plt4 <- ggplot(pdata4,aes(x,y))+geom_point(aes(colour=coll,shape=coll))+theme_test()+
  xlab('p*(0)')+ylab('p*')+ggtitle('(d)')+xlim(0,5)+ylim(0,5)+
  geom_abline(intercept=0,slope=1)+geom_abline(intercept=0,slope=2,lty=2)+
  theme(plot.title=element_text(hjust=0.5,size=10),legend.text.align=0.3,
        strip.background=element_blank(),strip.placement='outside')+
  scale_shape_manual(name='Quantile(SSR)',
                     values=c(15,16,17,18),
                     labels=c('>80%','50-80%','25-50%','<25%'))+
  scale_colour_manual(name='Quantile(SSR)',
                      values=c('grey','dodgerblue','orange','darkred'),
                      labels=c('>80%','50-80%','25-50%','<25%'))

#Combine the individual plots into one figure
plt <- plt1+plt2+plt3+plt4+plot_layout(ncol=2,guides='collect') & ylab(NULL)# & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
print(wrap_elements(plt)+labs(tag='The estimated average turnover rate, p*')+
        theme(plot.tag=element_text(size=rel(1),angle=90),plot.tag.position='left'))

# save.image('Fig3.RData')
# ggsave('Fig3.pdf',width=8,height=7)
# ggsave('Fig3.svg',width=8,height=7)
# ggsave('Fig3.png',width=8,height=7)
