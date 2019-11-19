## Model calculations

## === utilities
logit <- function(x) log(x/(1-x))
ilogit <- function(x) 1/(exp(-x)+1)
relogit <- function(x,RR) ilogit(RR*logit(x))
odds <- function(x) (x/(1-x))
iodds <- function(x) x/(x+1)
reodds <- function(x,RR) iodds(RR*odds(x))

bracket <- function(x,y,z,digits=0,spc=FALSE,sep='comma'){
  x <- format(round(x,digits = digits),big.mark=',')
  y <- format(round(y,digits = digits),big.mark=',')
  z <- format(round(z,digits = digits),big.mark=',')
  if(sep=='comma'){
    if(spc)
      ans <- paste0(x,' (',y,' , ',z,')')
    else
      ans <- paste0(x,' (',y,',',z,')')
  } else{
    if(spc)
      ans <- paste0(x,' (',y,' - ',z,')')
    else
      ans <- paste0(x,' (',y,'-',z,')')
  }
  ans
} 
## bracket(0.2,0.1,0.9,digits=1)

sq <- 0.674                             #75 percentile limits are this many SDs


## ==== model logic
## construct tree from org file
DT <- org2tree('LTBI.org')              #NB this needs data.tree

## ==== model parameters
## check parm
DT$Set(check=0)
DT$Set(check=1,filterFun=isLeaf)
## Counters
DT$Set(TBcount=0)                            #TB
DT$Set(TBcount=1,filterFun=function(x) x$name=='TB' )
DT$Set(AEcount=0)                            #AE
DT$Set(AEcount=1,filterFun=function(x) x$name=='AE' )
DT$Set(PTcount=0)                            #PT
DT$Set(PTcount=1,filterFun=function(x) x$name=='PT' )
DT$HHC$p <- 1
## TST prev
DT$HHC$TSTp$p <- "LTBI"
DT$HHC$TSTn$p <- "1-LTBI"
## AE probs
DT$Set(p='aeprob',filterFun=function(x) x$name=='AE' )
DT$Set(p='1-aeprob',filterFun=function(x) x$name=='noAE' )
## no TB probs
DT$Set(p='tbprob',filterFun=function(x) x$name=='TB' & rev(x$path)[2]=='noPT') #parent nm
DT$Set(p='1-tbprob',filterFun=function(x) x$name=='noTB' & rev(x$path)[2]=='noPT')
DT$Set(p='tbprob',filterFun=function(x) x$name=='TB' & rev(x$path)[2]=='AE') #AE=no PT
DT$Set(p='1-tbprob',filterFun=function(x) x$name=='noTB' & rev(x$path)[2]=='AE')
DT$Set(p='PTtbprob',filterFun=function(x) x$name=='TB' & rev(x$path)[2]=='noAE') #noAE=PT
DT$Set(p='1-PTtbprob',filterFun=function(x) x$name=='noTB' & rev(x$path)[2]=='noAE')
## PT cov
DT$HHC$TSTp$PT$p <- "posPTcov"
DT$HHC$TSTp$noPT$p <- "1-posPTcov"
DT$HHC$TSTn$PT$p <- "negPTcov"
DT$HHC$TSTn$noPT$p <- "1-negPTcov"

## inspect
print(DT,'p','TBcount','AEcount','PTcount','check')

## ## plot
## plotter(DT,varz = c('name'))

## ==== make model functions
## get functions for calculation
F <- makeTfuns(DT,c("check","PTcount","AEcount","TBcount"))
F


## nsim x regs x acat x (2 LTBI)
## base DF
regs <- c("none","6H","4R","3HP")
nsim <- 1e4                             #change N here
N <- length(regs) * length(acat) * 2
DF <- expand.grid(
  LTBI=c(0,1),
  regimen=regs,
  acat=acat
)
DF <- as.data.table(DF)
## DF
## extend by nsim
DF <- DF[rep(1:N,each=nsim)]
DF[,iter:=rep(1:nsim,N)]
## DF

## add other trivial parameters
DF[,c('posPTcov','negPTcov'):=1] 

## TB progression: depends age
DF[,tbprob:=.1*(1+1e-2*runif(nrow(DF)))] #safety
for(ag in acat){
  nn <- nrow(DF[acat==ag])
  DF[acat==ag,tbprob:=P[[paste0("P.",ag)]]$r(nn)]
  ## if(ag=="0-17"){
  ##   DF[acat==ag,
  ##      tbprob:=pmax(0, (1-exp(-rnorm(nn,0.7704387,0.1819483)) )/3+
  ##                      (1-exp(-rnorm(nn,0.2720844,0.06204593)))*2/3)]
  ## } else {
  ##   DF[acat==ag,
  ##      tbprob:=pmax(0, (1-exp(-rnorm(nn,0.01985839,0.009111896))))] #1yr
  ## }
}

## LTBI: make TST-ve less likely to progress
nn <- nrow(DF[LTBI==0])
DF[LTBI==0,tbprob:=tbprob * P[['RR.TST']]$r(nn)]


## efficacy  
DF[regimen=="none",PTtbprob:=tbprob]
for(reg in regs[regs!="none"]){
  nn <- nrow(DF[regimen==reg])
  DF[regimen==reg,PTtbprob:=reodds(tbprob,P[[paste("OR",reg,sep=".")]]$r(nn))]
}

## ## adjust PT efficacy of TST-ve (taken==1)
## nn <- nrow(DF[LTBI==0 & regimen=="none"])
## DF[LTBI==0 & regimen!="none",PTtbprob:=PTtbprob/P[["RR.PT.TST"]]$r(nn)]

## AE probs
for(reg in regs[regs!="none"]){
  for(age in acat){
    nn <- nrow(DF[regimen==reg & acat==age])
    DF[regimen==reg & acat==age,aeprob:=P[[paste("AE",reg,age,sep=".")]]$r(nn)]
  }
}
## DF[regimen=="none",aeprob:=1e-4*runif(nrow(DF[regimen=="none"]))]
DF[regimen=="none",aeprob:=0]

## DF

## Apply model
DF[,check:=F$checkfun(DF)]
DF[,PTcount:=F$PTcountfun(DF)]
DF[,AEcount:=F$AEcountfun(DF)]
DF[,TBcount:=F$TBcountfun(DF)]
## save(DF,file='DF.Rdata')                #save
DF[,summary(check)]

## summary table
DFS <- DF[,.(AE = 1e3*mean(AEcount),TB = 1e3*mean(TBcount),
             AE.lo = 1e3*quantile(AEcount,.25),TB.lo = 1e3*quantile(TBcount,.25),
             AE.hi = 1e3*quantile(AEcount,.75),TB.hi = 1e3*quantile(TBcount,.75)),
          by=.(LTBI,acat,regimen)]
DFS

DFT <- dcast(DFS,LTBI + acat ~ regimen,value.var = c('AE','TB'))
knitr::kable(DFT)

cbbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
DF$regimen <- factor(DF$regimen,levels=regs)

## summary plot
## print("Saving plot...")
GP <- ggplot(DF,aes(AEcount,TBcount,col=regimen)) +
  geom_point(alpha=.1,size=0.5) +
  scale_color_manual(values=cbbPal)+
  facet_grid(LTBI ~ acat)               #lots of points/time


## ggsave(GP,file='GP1.png',width=15)

## GP <- ggplot(DF,aes(AEcount,TBcount,col=regimen)) +
##   stat_ellipse() +
##   facet_grid(LTBI ~ acat)               #needs jitter on AE

## GP <- ggplot(DF,aes(AEcount,TBcount,col=regimen)) +
##   facet_grid(LTBI ~ acat) + geom_density_2d() #needs jitter on none AE...

## summarized data
DFR <- DF[,
          .(AEcount=mean(AEcount),
            AEcount.sd=sd(AEcount),
            TBcount=mean(TBcount),
            TBcount.sd=sd(TBcount)),
          by = .(LTBI,acat,regimen)]

DFR[,TST:=paste0('TST=',LTBI)]
DFR[TST=='TST=0',TST:='TST-negative']
DFR[TST=='TST=1',TST:='TST-positive']


DFR$regimen <- factor(DFR$regimen,levels=regs)
rotx <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot rexduced
GP <- ggplot(DFR,aes(x=AEcount*1e3,xmin=(AEcount-AEcount.sd*sq)*1e3,
                     xmax=(AEcount+AEcount.sd*sq)*1e3,
                     y=TBcount*1e3,ymin=(TBcount-TBcount.sd*sq)*1e3,
                     ymax=(TBcount+TBcount.sd*sq)*1e3,
                     col=regimen)) +
  geom_point(size=2) + scale_color_manual(values=cbbPal)+
  geom_errorbar(width=0) + geom_errorbarh(height=0)+
  ## scale_x_log10() + scale_y_log10() +
  facet_grid(TST ~ acat)

GP <- GP + xlab('Severe adverse event risk (per 1000)') +
  ylab('TB disease progression risk (per 1000)')
GP <- GP + rotx

Figure2 <- GP + theme_bw()

## ggsave(GP,file='GP2.png',width=15)


## looking at differences
DFW <- dcast(DF,LTBI + acat + iter ~ regimen,value.var = c('AEcount','TBcount'))
aecols0 <- grep('AE',names(DFW),value=TRUE)
aecols <- aecols0[aecols0!="AEcount_none"]
tbcols0 <- grep('TBcount',names(DFW),value=TRUE)
tbcols <- tbcols0[tbcols0!="TBcount_none"]
DFW[, paste0("d", aecols) := .SD - DFW$AEcount_none, .SDcols=aecols]
DFW[, paste0("d", tbcols) := .SD - DFW$TBcount_none, .SDcols=tbcols]
DFW[,c(aecols0,tbcols0):=NULL]
DFWS <- melt(DFW,id.vars = c('LTBI','acat','iter'))

DFWM <- DFWS[,.(mid=mean(value),lo=quantile(value,.25),hi=quantile(value,.75)),
             by=.(LTBI,acat,variable)]


DFWM[,regimen:= unlist(lapply(variable,
                              function(x) unlist(
                                            strsplit(as.character(x),split="_"))[2]))]

DFWM[,type:='TB']
DFWM[grepl('AE',variable),type:='AE']
DFWM[,variable:=NULL]
DFWM2 <- dcast(DFWM,LTBI + acat + regimen ~ type,value.var = c('mid','lo','hi'))
DFWM2[,TST:=paste0('TST=',LTBI)]

GP2 <- ggplot(DFWM2,aes(x=mid_AE*1e3,xmin=lo_AE*1e3,xmax=hi_AE*1e3,
                        y=mid_TB*1e3,ymin=lo_TB*1e3,ymax=hi_TB*1e3,
                        col=regimen)) +
  geom_errorbar(width=0) + geom_errorbarh(height=0) +
  geom_point() +
  scale_color_manual(values=cbbPal[-1])+
  facet_grid(TST~acat) +
  ylab('Change in TB cases (per 1000)') +
  xlab('Change in severe adverse event risk (per 1000)')


DFWM2[,NNH:=1/abs(mid_AE)];
## DFWM2[,NNHlo:=1/abs(hi_AE)]; DFWM2[,NNHhi:=1/abs(lo_AE)];
DFWM2[,NNT:=1/abs(mid_TB)];
## DFWM2[,NNThi:=1/abs(hi_TB)]; DFWM2[,NNTlo:=1/abs(lo_TB)];

GP3 <- ggplot(DFWM2,aes(x=NNH,
                        y=NNT,
                        col=regimen,group=acat)) +
  ## geom_errorbar(aes(ymin=NNTlo,ymax=NNThi),width=0,alpha=.3) +
  ## geom_errorbarh(aes(xmin=NNHlo,xmax=NNHhi),height=0,alpha=.3) +
  ylim(c(0,200)) + xlim(c(0,400))+
  geom_point() + geom_text_repel(aes(label=acat))+
  geom_line(col=1,lty=2,alpha=.5)+
  scale_color_manual(values=cbbPal[-1])+
  facet_grid(~TST) +
  geom_abline(intercept=0,slope = 1,col='grey',lty=2)+
  annotate(x=200,y=200-14,geom="text",label="NNT=NNH",col='grey')+
  geom_abline(intercept=0,slope = 2,col='grey',lty=2)+
  annotate(x=100,y=200-7,geom="text",label="NNT=2.NNH",col='grey')+
  geom_abline(intercept=0,slope = 3,col='grey',lty=2)+
  annotate(x=200/3,y=200,geom="text",label="NNT=3.NNH",col='grey')+
  xlab('NNH through an extra AE') + ylab('NNT to avert a TB case')
## GP3

## looking at difference LTBI+/-
DFLW <- dcast(DF,acat + iter ~ LTBI + regimen,value.var = c('AEcount','TBcount'))
aecols0 <- grep('AEcount_0',names(DFLW),value=TRUE)
aecols01 <- aecols0[aecols0!="AEcount_0_none"]
aecols1 <- grep('AEcount_1',names(DFLW),value=TRUE)
aecols11 <- aecols1[aecols1!="AEcount_1_none"]
tbcols0 <- grep('TBcount_0',names(DFLW),value=TRUE)
tbcols01 <- tbcols0[tbcols0!="TBcount_0_none"]
tbcols1 <- grep('TBcount_1',names(DFLW),value=TRUE)
tbcols11 <- tbcols1[tbcols1!="TBcount_1_none"]

DFLW[, paste0("d", aecols01) := .SD - DFLW$AEcount_0_none, .SDcols=aecols01]
DFLW[, paste0("d", aecols11) := .SD - DFLW$AEcount_1_none, .SDcols=aecols11]
DFLW[, paste0("d", tbcols01) := .SD - DFLW$TBcount_0_none, .SDcols=tbcols01]
DFLW[, paste0("d", tbcols11) := .SD - DFLW$TBcount_1_none, .SDcols=tbcols11]
DFLW[,c(aecols0,tbcols0,aecols1,tbcols1):=NULL]


## addin extra cols
DFL <- melt(DFLW,id.vars = c('acat','iter'))
DFL[,qty:='AE']
DFL[grepl('TB',variable),qty:='TB']
DFL[,LTBI:=1]
DFL[grepl('0',variable),LTBI:=0]
DFL[,regimen:= unlist(lapply(variable,
                             function(x) unlist(
                                           strsplit(as.character(x),split="_"))[3]))]
DFL[,variable:=NULL]
DFLW <- dcast(DFL,acat+iter+qty+regimen ~ LTBI,value.var = 'value')
## names(DFLW)[5:6] <- c('Lp','Ln')
names(DFLW)[names(DFLW)=='0'] <- 'Ln'
names(DFLW)[names(DFLW)=='1'] <- 'Lp'

## TS <- copy(DFLW)                        #for potential later use

## 10%,25%,50%,75% averages
DFLW[,pc10_n:=1e-1*Lp+(1-1e-1)*Ln]
DFLW[,pc10_y:=1e-1*Lp+0*(1-1e-1)*Ln]
DFLW[,pc25_n:=25e-2*Lp+(1-25e-2)*Ln]
DFLW[,pc25_y:=25e-2*Lp+0*(1-25e-2)*Ln]
DFLW[,pc50_n:=5e-1*Lp+(1-5e-1)*Ln]
DFLW[,pc50_y:=5e-1*Lp+0*(1-5e-1)*Ln]
DFLW[,pc75_n:=75e-2*Lp+(1-75e-2)*Ln]
DFLW[,pc75_y:=75e-2*Lp+0*(1-75e-2)*Ln]
DFLW[,c('Lp','Ln'):=NULL]

DLM <- melt(DFLW,id.vars = c('acat','iter','qty','regimen'))

DLM[,TST:=sub(".*_", "", as.character(variable))]
DLM[,LTBI:=unlist(str_extract_all(variable,"[0-9,.]+"))]
DLM[,LTBI:=paste0(LTBI," %")]
DLM[,variable:=NULL]

DLM <- DLM[,.(mid=mean(value),hi=quantile(value,.75),lo=quantile(value,.25)),
           by=.(acat,qty,regimen,TST,LTBI)]
DLM$LTBI <- factor(DLM$LTBI,levels=paste0(c(10,25,50,75)," %"),ordered = TRUE)

## differences in mean risks
GP4 <- ggplot(DLM,aes(x=TST,y=mid,ymin=lo,ymax=hi,
                      col=acat,shape=regimen,group=paste(acat,regimen)))+
  geom_point(size=2) +
  facet_grid(qty~LTBI,scales='free_y') +
  geom_line(lty=2) +
  geom_errorbar(alpha=.2,width=.1) +
  ylab('Mean risk difference in contacts')
GP4


GP4 <- ggplot(DLM[qty=='AE'],aes(x=LTBI,y=mid,ymin=lo,ymax=hi,
                      col=TST,shape=regimen,group=paste(acat,regimen,TST)))+
  geom_point(size=2) +
  facet_grid(qty~acat,scales='free_y') +
  geom_line(lty=2) +
  geom_errorbar(alpha=.2,width=.1) +
  ylab('Mean risk difference in contacts')
GP4

GP4 <- ggplot(DLM[qty=='TB'],aes(x=LTBI,y=mid,ymin=lo,ymax=hi,
                      col=TST,shape=regimen,group=paste(acat,regimen,TST)))+
  geom_point(size=2) +
  facet_grid(qty~acat,scales='free_y') +
  geom_line(lty=2) +
  geom_errorbar(alpha=.2,width=.1) +
  ylab('Mean risk difference in contacts')
GP4


## looking specifically at difference between using test and not
DLM3 <- copy(DFLW)
DLM3[,pc10:=pc10_y-pc10_n]
DLM3[,pc25:=pc25_y-pc25_n]
DLM3[,pc50:=pc50_y-pc50_n]
DLM3[,pc75:=pc75_y-pc75_n]

dp <- c('pc10','pc25','pc50','pc75')
DLM3[,c(paste(dp,'y',sep="_"),paste(dp,'n',sep="_")):=NULL]
DLM3 <- melt(DLM3,id.vars = names(DLM3)[1:4])
DLM3[,LTBI:="10 %"]
DLM3[grep('25',variable),LTBI:="25 %"]
DLM3[grep('50',variable),LTBI:="50 %"]
DLM3[grep('75',variable),LTBI:="75 %"]
DLM3[,variable:=NULL]
DLM4 <- DLM3[,.(mid=mean(value),lo=quantile(value,0.25),hi=quantile(value,0.75)),
             by=.(acat,qty,regimen,LTBI)]
DLM4$LTBI <- factor(DLM4$LTBI,levels=paste0(c(10,25,50,75)," %"),ordered = TRUE)

DLM4

DLM4a <- copy(DLM4)
## DLM4a[qty=='AE',mid:=-mid]
## DLM4a[qty=='AE',lo:=-lo]
## DLM4a[qty=='AE',hi:=-hi]
DLM4a[,mid:=-mid]
DLM4a[,lo:=-lo]
DLM4a[,hi:=-hi]


DLM4a[,QTY:='Tuberculosis cases']
DLM4a[qty=='AE',
      QTY:='Severe adverse events']
DLM4a$QTY <- factor(DLM4a$QTY,levels=sort(unique(DLM4a$QTY),dec=TRUE),ordered=TRUE)

GP4ab <- ggplot(DLM4a,
                aes(x=acat,y=mid*1e3,ymin=lo*1e3,ymax=hi*1e3,
                    col=LTBI,shape=regimen,group=paste(LTBI,regimen)))+
  geom_point(size=2) +
  facet_grid(regimen~QTY) +
  geom_line(lty=2) +
  geom_errorbar(alpha=.2,width=.1) +
  ylab('Change in outcomes per 1,000 contacts by not using TST') +
  xlab('Age group') +
  scale_color_viridis_d(name='TST positivity')+
  geom_hline(yintercept=0,col='grey',lty=1) + theme_bw()
## GP4ab

Figure3 <- GP4ab + rotx

GP4ab2 <- ggplot(DLM4a,
                aes(x=acat,y=mid*1e3,ymin=lo*1e3,ymax=hi*1e3,
                    col=LTBI,shape=regimen,group=paste(LTBI,regimen)))+
  geom_point(size=2) +
  facet_grid(QTY~regimen) +
  geom_line(lty=2) +
  geom_errorbar(alpha=.2,width=.1) +
  ylab('Change in outcomes per 1,000 contacts by not using TST') +
  xlab('Age group') +
  scale_color_viridis_d(name='TST positivity')+
  geom_hline(yintercept=0,col='grey',lty=1) + theme_bw()
## GP4ab2
Figure3v2 <- GP4ab2 + rotx

GP4abc <- ggplot(DLM4a,aes(x=acat,y=1/mid,
                      col=LTBI,shape=regimen,group=paste(LTBI,regimen)))+
  geom_point(size=2) +
  facet_grid(regimen~QTY) +
  geom_line(lty=2) +
  ylab('Number of contacts per event, from not using TST')
## GP4abc

fwrite(DLM4a[,.(acat,qty,regimen,LTBI,mid,lo,hi)],file=here('figures/TableS4a.csv'))
## flip hi/lo for negatives to make sense
DLM4a[lo<0 & hi<0, c('tmplo','tmphi'):=.(hi,lo)]
DLM4a[lo<0 & hi<0, c('lo','hi'):=.(tmplo,tmphi)]

fwrite(cbind(DLM4a[,.(acat,qty,regimen,LTBI)],
             value=bracket(1e3*DLM4a$mid,1e3*DLM4a$lo,1e3*DLM4a$hi,digits=1)),
       file=here('figures/TableS4b.csv'))


DLM5 <- dcast(DLM4,acat+regimen+LTBI~qty,value.var = c('mid','lo','hi'))

GP4b <- ggplot(DLM5,aes(x=-mid_AE*1e3,y=mid_TB*1e3,#ymin=lo,ymax=hi,
                        col=acat,shape=regimen,group=paste(acat,regimen),
                        alpha=LTBI))+
  geom_point(size=2)  +
  geom_line()  +
  xlab('Increased serious adverse events per 1000 contacts by not using TST') +
  ylab('Decreased tuberculosis cases per 1000 contacts by not using TST') +
  geom_abline(intercept=0,slope = 1,lty=2) +
  annotate(geom='text',x=10,y=10,label='TB ~ AE') +
  geom_abline(intercept=0,slope = 0.5,lty=2,col='grey') +
  annotate(geom='text',x=20,y=8,label='TB ~ 2 x AE',col='grey') +
  annotate(geom='text',x=0.00,y=0.0,label='+') + 
  theme_classic() +
  theme(legend.position = c(0.5,0.8), legend.box = "horizontal") +
  xlim(c(0,NA)) + ylim(c(0,NA))
GP4b


## DF
DFLW0 <- dcast(DF,acat + regimen + iter ~ LTBI,value.var = c('AEcount','TBcount'))
prev <- 0.25

## needs to compare to no PT
DFLW0 <- merge(DFLW0,DFLW0[regimen=='none',.(acat,iter,TBcount_0bc=TBcount_0)],
               by=c('acat','iter'),all.x = TRUE)

DFLW0[,aecount_all:=AEcount_0*1e3]
DFLW0[,aecount_tst:=prev*AEcount_0*1e3]
DFLW0[,tbcount_all:=((1-prev)*TBcount_0 + prev*TBcount_1)*1e3]
DFLW0[,tbcount_tst:=((1-prev)*TBcount_0bc + prev*TBcount_1)*1e3]

Table3 <- DFLW0[,.(all_ae.mid=mean(aecount_all),
                   all_ae.lo=quantile(aecount_all,.25),
                   all_ae.hi=quantile(aecount_all,.75),
                   all_tb.mid=mean(tbcount_all),
                   all_tb.lo=quantile(tbcount_all,.25),
                   all_tb.hi=quantile(tbcount_all,.75),
                   tst_ae.mid=mean(aecount_tst),
                   tst_ae.lo=quantile(aecount_tst,.25),
                   tst_ae.hi=quantile(aecount_tst,.75),
                   tst_tb.mid=mean(tbcount_tst),
                   tst_tb.lo=quantile(tbcount_tst,.25),
                   tst_tb.hi=quantile(tbcount_tst,.75)
                   ),by=.(regimen,acat)]
Table3 <- Table3[order(regimen,acat)]

Table3pretty <- cbind(Table3[,.(regimen,acat)],
                      all_tb=Table3[,bracket(all_tb.mid,all_tb.lo,all_tb.hi,digits=1)],
                      all_ae=Table3[,bracket(all_ae.mid,all_ae.lo,all_ae.hi,digits=1)],
                      tst_tb=Table3[,bracket(tst_tb.mid,tst_tb.lo,tst_tb.hi,digits=1)],
                      tst_ae=Table3[,bracket(tst_ae.mid,tst_ae.lo,tst_ae.hi,digits=1)]
                      )

fwrite(Table3pretty,file=here('figures/Table3.csv'))

## loop for supplementary tables (actually including 25%)
Table3X <- list()
for(pv in c(.1,.25,.5,.75)){
  tmp <- dcast(DF,acat + regimen + iter ~ LTBI,value.var = c('AEcount','TBcount'))
  tmp <- merge(tmp,tmp[regimen=='none',.(acat,iter,TBcount_0bc=TBcount_0)],
               by=c('acat','iter'),all.x = TRUE)
  tmp[,aecount_all:=AEcount_0*1e3]
  tmp[,aecount_tst:=pv*AEcount_0*1e3]
  tmp[,tbcount_all:=((1-pv)*TBcount_0 + pv*TBcount_1)*1e3]
  tmp[,tbcount_tst:=((1-pv)*TBcount_0bc + pv*TBcount_1)*1e3]
  tmptab <- tmp[,.(all_ae.mid=mean(aecount_all),
                   all_ae.lo=quantile(aecount_all,.25),
                   all_ae.hi=quantile(aecount_all,.75),
                   all_tb.mid=mean(tbcount_all),
                   all_tb.lo=quantile(tbcount_all,.25),
                   all_tb.hi=quantile(tbcount_all,.75),
                   tst_ae.mid=mean(aecount_tst),
                   tst_ae.lo=quantile(aecount_tst,.25),
                   tst_ae.hi=quantile(aecount_tst,.75),
                   tst_tb.mid=mean(tbcount_tst),
                   tst_tb.lo=quantile(tbcount_tst,.25),
                   tst_tb.hi=quantile(tbcount_tst,.75)
                   ),by=.(regimen,acat)]
  tmptab <- tmptab[order(regimen,acat)]
  tmptab <- cbind(tmptab[,.(regimen,acat)],
                  all_tb=tmptab[,bracket(all_tb.mid,all_tb.lo,all_tb.hi,digits=1)],
                  all_ae=tmptab[,bracket(all_ae.mid,all_ae.lo,all_ae.hi,digits=1)],
                  tst_tb=tmptab[,bracket(tst_tb.mid,tst_tb.lo,tst_tb.hi,digits=1)],
                  tst_ae=tmptab[,bracket(tst_ae.mid,tst_ae.lo,tst_ae.hi,digits=1)]
                  )
  Table3X[[as.character(100*pv)]] <- tmptab
}

## write out
for(pv in names(Table3X)) fwrite(Table3X[[pv]],
                                 file=here('figures',paste0('Table3X_',pv,'.csv')))


DFLW0 <- merge(DFLW0,CT,by='acat',all.x=TRUE)

## all-tst
DFLW0[,aecount_lb:=1e3*LBprop*AEcount_0*(1-LBtst)]
DFLW0[,aecount_hb:=1e3*Hbprop*AEcount_0*(1-HBtst)]
DFLW0[,tbcount_lb:=1e3*LBprop*(1-LBtst)*(TBcount_0 - TBcount_0bc )]
DFLW0[,tbcount_hb:=1e3*Hbprop*(1-HBtst)*(TBcount_0 - TBcount_0bc )]

## DFLW0[,aecount_all:=AEcount_0*1e3]
## DFLW0[,aecount_tst:=prev*AEcount_0*1e3]
## DFLW0[,tbcount_all:=((1-prev)*TBcount_0 + prev*TBcount_1)*1e3]
## DFLW0[,tbcount_tst:=((1-prev)*TBcount_0bc + prev*TBcount_1)*1e3]

TS <- DFLW0[regimen!='none']
TS <- TS[,.(aecount_lb=sum(aecount_lb),aecount_hb=sum(aecount_hb),
            tbcount_lb=sum(tbcount_lb),tbcount_hb=sum(tbcount_hb)),
         by=.(regimen,iter)]

Table4 <- TS[,.(lb.ae.mid=mean(aecount_lb),
                lb.ae.lo=quantile(aecount_lb,.25),
                lb.ae.hi=quantile(aecount_lb,.75),
                lb.tb.mid=mean(tbcount_lb),
                lb.tb.lo=quantile(tbcount_lb,.25),
                lb.tb.hi=quantile(tbcount_lb,.75),
                hb.ae.mid=mean(aecount_hb),
                hb.ae.lo=quantile(aecount_hb,.25),
                hb.ae.hi=quantile(aecount_hb,.75),
                hb.tb.mid=mean(tbcount_hb),
                hb.tb.lo=quantile(tbcount_hb,.25),
                hb.tb.hi=quantile(tbcount_hb,.75)),
             by=regimen]


s <- TRUE
Table4pretty <- cbind(Table4[,.(regimen)],
                      lb_tb=Table4[,bracket(lb.tb.mid,lb.tb.lo,lb.tb.hi,digits=1,s)],
                      lb_ae=Table4[,bracket(lb.ae.mid,lb.ae.lo,lb.ae.hi,digits=1,s)],
                      hb_tb=Table4[,bracket(hb.tb.mid,hb.tb.lo,hb.tb.hi,digits=1,s)],
                      hb_ae=Table4[,bracket(hb.ae.mid,hb.ae.lo,hb.ae.hi,digits=1,s)]
                      )
fwrite(Table4pretty,file=here('figures/Table4.csv'))
