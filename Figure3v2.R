## improvements to figure 3 under review
library(here)
library(ggplot2)
library(ggrepel)
library(viridis)

load(here('figures/DLM4a.Rdata'))

rotx <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

pd <- position_dodge(width=0.5)

GP4ab <- ggplot(DLM4a,
                aes(x=acat,y=mid*1e3,ymin=lo*1e3,ymax=hi*1e3,
                    col=LTBI,shape=regimen,group=paste(LTBI,regimen)))+
  geom_point(size=2,position=pd) +
  facet_grid(regimen~QTY) +
  geom_line(lty=1,alpha=.5,position=pd) +
  geom_errorbar(alpha=.99,width=.2,position=pd) +
  ylab('Change in outcomes per 1,000 contacts by not using TST') +
  xlab('Age group') +
  scale_color_viridis_d(name='TST positivity')+
  geom_hline(yintercept=0,col='grey',lty=1) +
  theme_bw() + rotx# + ggpubr::grids()
GP4ab

ggsave(GP4ab,file=here('figures/Figure3v3.pdf'))

ggsave(GP4ab,file=here('figures/Figure3v3.eps'),device = cairo_ps)
