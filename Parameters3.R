#' ---
#' title: Analysis to create parameters for PTR/B
#' author: Pete Dodd
#' date: October, 2019
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' 
#' #Pre-amble
#' (Last updated: `r Sys.time()`)
#' 
#' This describes analyses of various sources of data to inform decision model parameters
#'
#' N.B. This file should render to html with `rmarkdown::render('Parameters3.R')` or from the command line with ` R -q -e "rmarkdown::render(\"Parameters3.R\")"`
#'
#' ## Load libraries
#'
#' Relevant libraries for loading and manipulating data:
#'

rm(list=ls())
library(here)
library(data.table)
library(readxl)
library(HEdtree)
library(data.tree)
library(ggplot2)
library(ggrepel)
library(stringr)
library(viridis)
knitr::opts_chunk$set(cache=FALSE)

## === utilities
logit <- function(x) log(x/(1-x))
ilogit <- function(x) 1/(exp(-x)+1)
relogit <- function(x,RR) ilogit(RR*logit(x))
odds <- function(x) (x/(1-x))
iodds <- function(x) x/(x+1)
reodds <- function(x,RR) iodds(RR*odds(x))

#'
#' ## Load data
#' 
#' Look at Excel data:
XL <- read_excel("Parameters_v1.xlsx", sheet = "Parameters")
XL <- as.data.table(XL)
XX <- XL[!is.na(`Point estimate`),
         .(Parameter,
           Treatment,
           `Age group`,
           `Point estimate`,
           `95% lower`,
           `95% upper`,`Source`)]
XX

#' Age categories
acat <- c("0-17","18-34","35-64","65+")

#'
#' # Preventive therapy
#'
#' ## Regimen = 6H
#'
#' Change 9H to 6H
#'
XX[Treatment=='9H',Treatment:='6H']

#' 
#' First looing at AEs
#'
#' 

(XA <- XX[grepl('AE',Parameter) &
          grepl('6H',Treatment) &
         !is.na(`95% upper`),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:=paste0("AE.6H.",acat)]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:=c("6H AE probability, ages 0-17",
                   "6H AE probability, ages 18-34",
                   "6H AE probability, ages 35-64",
                   "6H AE probability, ages 65+")]

## get parms
tmp <- getAB(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2)
XA[,DISTRIBUTION:=paste0("B(",tmp$a,",",tmp$b,")")] #beta distributions
## top 1 has to exponential
(aa <- XA[1,`95% upper`])
(bb <- -log(1-.975)/aa)
XA[1,DISTRIBUTION:=paste0("E(",bb,")")] #beta distributions

#' Start data table of parameters
(PDT <- XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])

#'
#' Efficacy:
#' 
#'
(XA <- XX[grepl('Efficacy',Parameter) &
          grepl('6H',Treatment) &
         !is.na(`95% upper`),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:="OR.6H"]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:="6H efficacy TST+ as OR, all ages"]
tmp <- getLNparms(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2,med=FALSE)
curve(dlnorm(x,tmp$mu,tmp$sig))
XA[,DISTRIBUTION:=paste0("LN(",tmp$mu,",",tmp$sig,")")] #LN distributions
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])

#' 
#' ## Regimen = 4R
#' 
#'
(XA <- XX[grepl('AE',Parameter) &
          grepl('4R',Treatment) &
         !is.na(`95% upper`),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:=paste0("AE.4R.",acat)]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:=c("4R AE probability, ages 0-17",
                   "4R AE probability, ages 18-34",
                   "4R AE probability, ages 35-64",
                   "4R AE probability, ages 65+")]
## get parms
tmp <- getAB(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2)
XA[,DISTRIBUTION:=paste0("B(",tmp$a,",",tmp$b,")")] #beta distributions
## top 1 has to exponential
(aa <- XA[1,`95% upper`])
(bb <- -log(1-.95)/aa)
XA[1,DISTRIBUTION:=paste0("E(",bb,")")] #beta distributions

#' Add to data table of parameters
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])

#'
#' Efficacy:
#' 
(XA <- XX[grepl('Efficacy',Parameter) & grepl('4R',Treatment),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:="OR.4R"]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:="4R efficacy TST+ as OR, all ages"]
tmp <- getLNparms(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2,med=FALSE)
curve(dlnorm(x,tmp$mu,tmp$sig))
XA[,DISTRIBUTION:=paste0("LN(",tmp$mu,",",tmp$sig,")")] #LN distributions

#' Add to data table of parameters
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])



#' 
#' ## Regimen = 3HP
#' 
#'
(XA <- XX[grepl('AE',Parameter) &
          grepl('3HP',Treatment) &
         !is.na(`95% upper`),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:=paste0("AE.3HP.",acat)]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:=c("4R AE probability, ages 0-17",
                   "4R AE probability, ages 18-34",
                   "4R AE probability, ages 35-64",
                   "4R AE probability, ages 65+")]
## get parms
tmp <- getAB(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2)
XA[,DISTRIBUTION:=paste0("B(",tmp$a,",",tmp$b,")")] #beta distributions

#' Add to data table of parameters
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])

#'
#' Efficacy:
#' 
(XA <- XX[grepl('Efficacy',Parameter) & grepl('3HP',Treatment),
         .(`Treatment`, `Age group`, `Point estimate`,
           `95% lower`, `95% upper`,`Source`)])

## add relevant
XA[,NAME:="OR.3HP"]
XA[,SOURCE:=Source]
XA[,NOTE:=""]
XA[,DESCRIPTION:="3HP efficacy TST+ as OR, all ages"]
tmp <- getLNparms(XA[,`Point estimate`],(XA[,`95% upper`]-XA[,`95% lower`])^2/3.92^2,med=FALSE)
curve(dlnorm(x,tmp$mu,tmp$sig))
XA[,DISTRIBUTION:=paste0("LN(",tmp$mu,",",tmp$sig,")")] #LN distributions

#' Add to data table of parameters
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])


#'
#' # Efficacy of Preventive therapy efficacy: TST+ vs -
#'
#' NOTE current decision to consider this as 1.


#'
#' # Progression to TB disease
#'
#' ## TST+
#'
#'
#' Using Trauer's adjusted estimates, however including
#' 1 and 3 month cut-offs at the start time to exclude co-prevalent cases,
#' and 1 and 3 year end-points. Primary analysis decided as 3y/3mo.
#' 
## knitr::include_graphics('Trauer2.jpg')
#' (approximate digitization)
trauer <- data.table(age=rep(c('0-4','5-14','15+'),4),
                     p.hi=c(c(0.21,0.12,.006),
                       c(0.53,0.22,.016),
                       c(.71,.35,.04),
                       c(.71,.36,.05)),
                     p=c(c(0.18,.07,.001),
                         c(0.39,0.15,.006),
                         c(.56,.28,0.02),
                         c(.56,0.28,.03)),
                     p.lo=c(c(0.13,.025,0),
                            c(0.26,.09,0),
                            c(.41,.17,.01),
                            c(.41,0.18,.009)),
                     time=c(rep('1 month',3), #30 d
                            rep('3 month',3), #90 d
                            rep('1 year',3),  #365 d
                            rep('3 year',3))) #1095 d

fwrite(trauer[,.(time,age,p,p.lo,p.hi)],file=here('figures/S_TrauerData.csv'))

#' Add beta parameters and inspect
trauer[,c('a','b'):=getAB(p,(p.hi-p.lo)^2/3.92^2)]
knitr::kable(trauer)

#' Visualize distributions
#' 
## make data
nn <- 1e3
xx <- seq(from=0,to=1,len=nn)
TD <- as.data.table(expand.grid(x=xx,
                                age=trauer[,unique(age)],
                                time=trauer[,unique(time)]))
TD <- merge(TD,trauer[,.(age,time,a,b)],by=c('age','time'),all.x=TRUE)
TD[,PDF:=dbeta(x,shape1 = a,shape2 = b)]
TD$time <- factor(TD$time,levels=c('1 month','3 month','1 year','3 year'),ordered=TRUE)
TD$age <- factor(TD$age,levels=c('0-4','5-14','15+'),ordered=TRUE)

ggplot(TD,aes(x,PDF,col=time,group=paste(age,time))) +
  geom_line() +
  facet_wrap(~age,ncol=1,scales='free')  + scale_x_sqrt()

#' NOTE only 3 time curves are visible in some age groups due to completion of events meaning that 1 year risks are identical to 3 year risks
#'
#' How to approach the difference of these variables?
#'
#' Check reasonable to use mean/variance of differences:
## data
one <- rbeta(1e4,trauer[age=='0-4' & time=='1 month',a],
             trauer[age=='0-4' & time=='1 month',b])
two <- rbeta(1e4,trauer[age=='0-4' & time=='1 year',a],
             trauer[age=='0-4' & time=='1 year',b])

## mean & var
mn <- trauer[age=='0-4' & time=='1 year',a/(a+b)] - trauer[age=='0-4' & time=='1 month',a/(a+b)]
vn <- trauer[age=='0-4' & time=='1 year',a*b/(a+b)^2/(a+b+1)] + trauer[age=='0-4' & time=='1 month',a*b/(a+b)^2/(a+b+1)]
AB <- getAB(mn,vn)

## plot
hist(two-one,freq=FALSE)
curve(dbeta(x,shape1=AB$a,shape2=AB$b),n=1e3,col=2,add=TRUE)

#' Make new data for difference combinations
trauer2 <- copy(trauer[,.(age,time,a,b)])
trauer2[,c('mn','vn'):=.(a/(a+b),a*b/(a+b)^2/(a+b+1))]
T2 <- dcast(trauer2[,.(age,time,mn,vn)],age ~ time,value.var = c('mn','vn'))

## make bits of new data
tz <- c('1 to 12 months','1 to 36 months','3 to 12 months','3 to 36 months')
AB11 <- T2[,getAB(`mn_1 year`-`mn_1 month`,`vn_1 year`+`vn_1 month`)]
AB13 <- T2[,getAB(`mn_1 year`-`mn_3 month`,`vn_1 year`+`vn_3 month`)]
AB31 <- T2[,getAB(`mn_3 year`-`mn_1 month`,`vn_3 year`+`vn_1 month`)]
AB33 <- T2[,getAB(`mn_3 year`-`mn_3 month`,`vn_3 year`+`vn_3 month`)]
AB11[,time:='1 to 12 months']; AB13[,time:='3 to 12 months'];
AB31[,time:='1 to 36 months']; AB33[,time:='3 to 36 months']
ABBA <- rbindlist(list(AB11,AB31,AB13,AB33))
ABBA[,age:=rep(T2$age,4)]
ABBA[,c('mn','vn'):=.(a/(a+b),a*b/(a+b)^2/(a+b+1))]
ABBA

#' Now need weighted average for under 18s
ABBA2 <- data.table(expand.grid(time=ABBA[,unique(time)],age=acat))
ABBA2 <- merge(ABBA2,ABBA[age=='15+',.(time,a,b)],by='time') #adults

## weight over children
AY <- ABBA[age=='0-4']; AO <- ABBA[age=='5-14']
AB <- merge(AY[,.(time,mny=mn,vny=vn)],AO[,.(time,mno=mn,vno=vn)],by='time')
AB[,c('mn','vn'):=.((mny+2*mno)/3,(vny+2*vno)/3)]
AB[,c('a2','b2'):=getAB(mn,vn)]
AB[,.(time,mny,mno,mn,a2,b2)]

## merge and clean
ABBA2 <- merge(ABBA2,AB[,.(time,a2,b2)],by='time',all.x = TRUE)
ABBA2[age=='0-17',c('a','b'):=.(a2,b2)]
ABBA2[,c('a2','b2'):=NULL]
ABBA2

#' Summary of mean progression probabilities for different variants
knitr::kable(dcast(ABBA2[,.(time,age,meanprob=a/(a+b))],age ~ time))

#' Introduce a labelling system for the different variants
vrnt <- c('11','13','31','33')          #label for different variants


#' Use these for progression risks
XA <- data.table(NAME=rep("",length(acat)*length(vrnt)),
                 SOURCE="Trauer et al",
                 NOTE="See DESCRIPTION for variant naming convention",
                 DESCRIPTION="Progression risk: age ",
                 DISTRIBUTION=""
                 )
XA[,NAME:=paste0(rep(paste0("P.",acat),each=4),'.',vrnt)]
XA[,DESCRIPTION:=paste0(DESCRIPTION,rep(acat,each=4))]
XA[,DESCRIPTION:=paste0(DESCRIPTION,'; ',substr(vrnt,2,2),' years - ',substr(vrnt,1,1),' months')]
XA[,y:=rep(substr(vrnt,2,2),4)]; XA[,m:=rep(substr(vrnt,1,1),4)]
XA <- XA[order(y,m)]                    #reorder to align with ABBA2
XA[,DISTRIBUTION:=paste0("B(",ABBA2$a,",",ABBA2$b,")")]
XA[,.(NAME,DISTRIBUTION)]


#' Note also there are mis-matches in age categories
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])


#' 
#' ## TST-
#'
#' Abubakar for 10mm
XT <- data.table(TST=c('-ve','+ve'),N=c(5293,2540),n=c(23,69))
XT
sdlnRR <- sqrt(1/XT[TST=='+ve',n] + 1/XT[TST=='-ve',n] - 1/XT[TST=='+ve',N] - 1/1/XT[TST=='-ve',N]) #sqrt(1/a + 1/c - 1/(a+b) - 1/(c+d))
RRtst <- c(mid=1,lo=exp(-1.96*sdlnRR),hi=exp(+1.96*sdlnRR)) * XT[TST=='-ve',n/N]/XT[TST=='+ve',n/N]
RRtst                                   #much lower than above
tmp2 <- getLNparms(RRtst['mid'],(RRtst['hi']-RRtst['lo'])^2/3.92^2,med=FALSE)

#' Choose which to use:
#' 
tmp <- tmp2
XA <- XA[1]
XA[,NAME:="RR.TST"]
XA[,SOURCE:="Abubakar"]
XA[,NOTE:="Under 10mm vs over; Abubakar Figure 1"]
XA[,DESCRIPTION:="Risk ratio for TB given TST-ve vs TST+ve, all ages"]
## curve(dlnorm(x,tmp$mu,tmp$sig))
XA[,DISTRIBUTION:=paste0("LN(",tmp$mu,",",tmp$sig,")")] #LN distributions
PDT <- rbind(PDT,XA[,.(NAME,DISTRIBUTION,DESCRIPTION,SOURCE,NOTE)])


#'
#' # Recording data
#' 
#' Write out clean input data:
fwrite(PDT,'RBAparms.csv')

#'
#' Read in the data (again), make parameter object, generate test outputs:
P <- parse.parmtable(data = read.csv("RBAparms.csv"),testdir = "testplots")
## P <- parse.parmtable(data = read.csv("RBAparms.csv")) # without plotting


#' Scenario table data
CT <- fread('Table2.csv')
CT <- CT[,lapply(.SD,function(x)1e-2*as.numeric(gsub('%','',x))),by=acat,.SDcols=2:5]
CT


#'Choose progression variant to use:
vt <- '33'                              #3 years - 3 month
touse <- grep(vt,names(P),value=TRUE)
newnm <- gsub(paste0('\\.',vt),'',touse)
for(i in 1:length(touse)) P[[newnm[i]]] <- P[[touse[i]]]
names(P)

#' Now write out a table for inclusion in text - ie join in via parm name
#' against the range data
PR <- fread('testplots/out_parmtable.csv')

PRO <- merge(PDT[,.(NAME,Parameter=DESCRIPTION,Source=SOURCE)],
             PR[,.(NAME=parm,`Estimate (75% uncertainty)` = mqrng)],
             by='NAME')

#' Drop parameters that are not relevant to the variant `vt`:
nvt <- vrnt[vrnt!=vt]
rgx <- paste0(nvt,collapse='|')
PRO <- PRO[!grepl(rgx,NAME)]


#' Write out
fwrite(PRO,file='Table1.csv')

#'
#' # Run model
#'
#'
#' This is setup in a different file, but prints the structure
source("RBAmodel3.R")

#'
#' # Inspect results
#'
#' Summary of means (calculated in other script, now per 1000):
knitr::kable(DFT)

#' Results for text around TB
#' 
DFS[LTBI==1][,.(acat,regimen,tb=round(TB),tb.lo=round(TB.lo),tb.hi=round(TB.hi))]

#' Results in text around AE
#' 
DFS[LTBI==1][,.(acat,regimen,ae=round(AE),ae.lo=round(AE.lo),ae.hi=round(AE.hi))]

#'
#' Plotting the above (generated in other script):
#+ out.width = '100%',fig.align="center"
Figure2

## Figure2 <- Figure2 + theme_bw() #+ ggpubr::grids(linetype = "dashed")
ggsave(Figure2,file='figures/Figure2.pdf',height = 5,width=10)

#'
#' In terms of risk differences wrt no PT, this looks like:
#+ out.width = '100%',fig.align="center"
GP2


#'
#'
#'
#' Table of NNH/NNT
#' 
knitr::kable(DFWM2[,.(LTBI,acat,regimen,
                      NNT=ceiling(NNT),
                      NNH=floor(NNH))])

#' Ways to examine the usefulness of the TST in various cases. First as risk differences
#' 
#+ out.width = '100%',fig.align="center"
Figure3

## Figure3 <- Figure3 + theme_classic() + ggpubr::grids(linetype = "dashed")
ggsave(Figure3,file=here('figures/Figure3.pdf'),w=7,h=7)
ggsave(Figure3v2,file=here('figures/Figure3v2.pdf'),w=7,h=7)


#' Figure for text
#'
#' 0-17 & 25% LTBI
DLM4a[acat=='0-17' & LTBI=='25 %'][,.(regimen,qty,
                                      mid=round(1e3*mid),
                                      lo=round(1e3*lo),
                                      hi=round(1e3*hi))]

#' >65 & 25% LTBI
DLM4a[acat=='65+' & LTBI=='25 %'][,.(regimen,qty,
                                    mid=round(1e3*mid),
                                    lo=round(1e3*lo),
                                    hi=round(1e3*hi))]

#'
#' One against other
#+ out.width = '100%',fig.align="center"
GP4b

ggsave(GP4b,file=here('figures/S_GP4b.pdf'))

#'
#' # Analyses for hypothetical populations
#'
#' Suggested table 3 (slightly modified)
knitr::kable(Table3pretty)

#' Suggested table 4
#' 


#' Output table
knitr::kable(Table4pretty)


