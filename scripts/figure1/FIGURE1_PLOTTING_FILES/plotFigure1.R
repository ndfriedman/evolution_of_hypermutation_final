#written by Noah Friedman
#the code to plot all the figures for figure 1
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)

#TO BEAUTIFY OUR PLOTS WE ADD "EMPTY THEME"
emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
#
###
######
###########
#######
###
#

plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/'

#PLOT FIGURE 1A

#MADE IN ILLUSTRATOR?

#
###
######
###########
#######
###
#

#PLOT FIGURE 1B

plot_n_cases_figure <- function(df){
  bottomScaleFactor <- 10
  p <- ggplot(df, aes(x=reorder(cancerType, -nTotal)))+
    geom_bar(aes(y=nTotal, fill='Total Cases'), stat='identity')+
    #optional include the n high mut burden as well
    geom_bar(aes(y=-1*bottomScaleFactor*nHypermutated - 1*bottomScaleFactor*nHighMutBurden, fill='High mutation burden'), stat='identity')+
    
    geom_bar(aes(y=-1*bottomScaleFactor*nHypermutated, fill='Hypermutated'), stat='identity')+
    
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    scale_fill_manual(values=c('#858585', 'black', 'light gray'))+
    #THIS IS MANUALLY SET YOU NEED TO CHANGE THE SCALE FACTOR IF YOU CHANGE IT
    scale_y_continuous(breaks= c(-4000, -2000, 0, 2000, 4000, 6000), labels=c('400', '200', '0', '2000', '4000', '6000'))+
    emptyTheme+
    guides(fill=guide_legend(title="Tumor Classification"))+
    xlab('Cancer Type')+
    ylab('Total Cases')
  return(p)
}

figure1bDataFrame <- read.table(paste(plottingFilePath, 'figure1bCancerTypeSummary.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_n_cases_figure(dfCancerTypes)
ggsave(paste(plottingFilePath, 'figure1b.pdf', sep=''), plot=p,  width = 6, height = 6)

#
###
######
###########
#######
###
#

#PLOT FIGURE 1C

plot_signatures_figure <- function(df){
  p <- ggplot(df, aes(x=reorder(cancerType, -orderingVal), y=nCases))+
    
    geom_bar(aes(fill=
                   factor(signature, levels=c('APOBEC', 'MMR',
                                              'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV')))
             , stat='identity')+
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    xlab('Cancer Type')+
    ylab('N Cases')+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
    guides(fill=guide_legend(title="Signature"))+
    emptyTheme
  return(p)
}

figure1cDataFrame <- read.table(paste(plottingFilePath, 'figure1cSignatureSummary.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_signatures_figure(dfSigs)
ggsave(paste(plottingFilePath, 'figure1c.pdf', sep=''),
       plot=p,  width = 5, height = 5)

#PIE CHART INSET
sigsInset <- ggplot(dfSigs, aes(x=1, y=nCases, fill=
                                factor(signature, levels=c('APOBEC', 'MMR',
                                                           'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV'))))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+
  guides(fill=guide_legend(title="Signature"))+
  scale_fill_manual(values=c("#FF0000","#267574",
                             'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
  emptyTheme+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  xlab('')+
  ylab('')
  
ggsave(paste(plottingFilePath, 'figure1c_pieChart.pdf', sep=''),
       plot=sigsInset,  width = 5, height = 5)


#
###
######
###########
#######
###
#

#PLOT FIGURE 1D

plot_data <- function(df){
  #plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nHotspots))+
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    #geom_boxplot(fatten = NULL, outlier.shape=NA)+
    geom_violin(bw = 1, aes(fill=factor(cohort,
                                  levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
                                  'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))))+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    #scale_colour_manual(values =  c('black', "#267574", 'gray', '#ADFF2F', "#9acd32", '#2A52BE'), name="Dominant\nSignature")+
    ylab('N Driver Mutations')+
    #ylab('N hotspot mutations')+
    xlab('Cancer Type')+
    emptyTheme+
    coord_cartesian(ylim=c(0,50))+
    #geom_jitter(aes(colour=factor(cohort,
    #                              levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
    #                                         'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))),
    #            shape=16, position=position_jitter(0.1), alpha=0.75)+
    
    scale_fill_manual(values =c('orange', '#b36200', 'lavender', '#301934', '#add8e6', 'blue', 'gray', '#333333'), name='Cohort')
  
  return(plt)
}

df <- read.table(paste(plottingFilePath, 'figure1d_nOncMutByCohort.tsv', sep=''), sep='\t', header=TRUE)
plt <- plot_data(df)
ggsave(paste(plottingFilePath, 'figure1d.pdf', sep=''),
       plot=plt,  width = 6, height = 4, units = c("in"), limitsize = FALSE)



#
####
#########
##############
#########
####
#

#plot figure 1e
make_obs_exp_plot <- function(includeLegend = TRUE){
  plt <- ggplot()+
    
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=obsHotspot, colour='Observed_MMR'), bins=10)+
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=expectedHotspot, colour='Expected_MMR'), bins=10)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=obsHotspot, colour='Observed_MMR'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=expectedHotspot, colour='Expected_MMR'), method='loess', se=FALSE, span =25)+
    
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=obsHotspot, colour='Observed_POLE'), bins=5)+
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=expectedHotspot, colour='Expected_POLE'), bins=5)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=obsHotspot, colour='Observed_POLE'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=expectedHotspot, colour='Expected_POLE'), method='loess', se=FALSE, span =25)+
    
    stat_summary_bin(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=obsHotspot, colour='Observed_TMZ'), bins=5)+
    stat_summary_bin(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=expectedHotspot, colour='Expected_TMZ'), bins=5)+
    geom_smooth(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=obsHotspot, colour='Observed_TMZ'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=expectedHotspot, colour='Expected_TMZ'), method='loess', se=FALSE, span =25)+
    
    xlab('N Nonsynonymous Mutations\nin IMPACT 341 Genes')+
    scale_x_continuous(breaks=c(0,100,200,400))+
    ylab('N hotspots per case')+
    #scale_x_log10()+
    #ggtitle('Observed and expected hotspot\nburden in hypermutated tumors')+
    emptyTheme+
    #scale_color_manual(values=c('gray', 'black'))+
    scale_color_manual(values=c("#9CA89C",
                                "gray", '#6699cc',"#267574","#ADFF2F", '#2A52BE', "#FFF600"))+
    labs(colour = "Number of Hotspots:")
  if(includeLegend == FALSE){
    plt <- plt + theme(legend.position = 'none')
  }
  return(plt)
}

dfObsExp <- read.table(paste(plottingFilePath, 'figure1e_observedVsExpected.tsv', sep=''), sep = '\t', header=TRUE)

plt <- make_obs_exp_plot(includeLegend = FALSE)
legend <- get_legend(make_obs_exp_plot(includeLegend = TRUE))

nCasesHistogram <- ggplot(dfObsExp, aes(x=nmut))+
  geom_histogram(bins=100)+
  scale_y_log10()+
  emptyTheme+
  xlab('')+
  ylab('N cases')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

alignedPlot <- plot_grid(nCasesHistogram, plt, nrow=2, rel_heights = c(.3, 1))
alignedPlotWithLegend <- plot_grid(alignedPlot, legend, ncol=2, rel_widths = c(1, .5))

ggsave(paste(plottingFilePath, 'figure1e.pdf', sep=''),
       plot=alignedPlotWithLegend,  width = 6, height = 5, units = c("in"))

#
####
#########
##############
#########
####
#

#plot figure 1f
df <- read.table(paste(plottingFilePath, 'figure1f_truncatingObsExp.tsv', sep=''),sep = '\t', header=TRUE)

s = 1
p <- ggplot()+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), method = 'lm', se=FALSE)+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), method = 'lm', se=FALSE)+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), method = 'lm', se=FALSE)+
  
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'),  method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'),  method = 'lm', se=FALSE)+
  
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'))+
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'))+
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'))+
  
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), bins=5)+
  
  
  xlab('TMB')+
  ylab('Difference between observed and\nexpected n truncating muts')+
  scale_color_manual(values=c('#013220', '#037D50', '#FF6347',
                             '#FF8C00', '#82CFFD', '#0D4F8B'))+
  labs(colour = 'Signature and\ngene type')+
  emptyTheme

ggsave(paste(plottingFilePath, 'figure1f.pdf', sep=''),
       plot=p,  width = 4, height = 5, units = c("in"))
