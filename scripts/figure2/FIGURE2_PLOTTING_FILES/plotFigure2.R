#written by Noah Friedman

plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/'
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)
library(ggrepel)

make_dnds_weak_drivers_plot <- function(dndsData, title){
  plotThresh <- .1
  minCoord = .9
  
  dndsData <- dndsData[(dndsData$qglobal_cv.Normal > .01) & (dndsData$qglobal_cv.Normal > .01),]
  ggplot()+
    geom_text_repel(data = dndsData[(dndsData$qglobal_cv.Normal < plotThresh) | (dndsData$qglobal_cv.Hypermutated < plotThresh),],
                    aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated), label=gene_name, colour=cancerType))+
    geom_point(data = dndsData[(dndsData$qglobal_cv.Normal >= plotThresh) & (dndsData$qglobal_cv.Hypermutated >= plotThresh),],
               aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated), colour=cancerType))+
    
    geom_point(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
                   y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), colour='orange', size=3)+
    #geom_text_repel(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
    #                    y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), label='mean q_val', colour='orange')+
    
    xlim(minCoord,3)+
    ylim(minCoord,3)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= 3), colour='black', linetype=2)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= 3), colour='black')+
    ylab('1 minus q value in hypermutators')+
    xlab('1 minus q value in non-hypermutators')+
    ggtitle(title)
}


#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2A

make_dnds_plot <- function(dndsData, title){
  
  
  emptyTheme <- theme(axis.line = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  
  capThresh <- 1e-9
  signifThresh = .0001
  plotThreshHyper <- .001
  plotThreshNormal <- .000001
  minCoord = .9
  
  maxCoord = 2-log10(capThresh)
  
  ggplot()+
    geom_text_repel(data = dndsData[((dndsData$qglobal_cv.Normal > signifThresh) & (dndsData$qglobal_cv.Hyper < plotThreshHyper)),],
                    aes(x=1- log10(qglobal_cv.Normal + capThresh), y=1-log10(qglobal_cv.Hypermutated+ capThresh),
                        label=gene_name, colour=cancerType))+
    geom_point(data = dndsData,
               aes(x=1- log10(qglobal_cv.Normal + capThresh), y=1-log10(qglobal_cv.Hypermutated + capThresh),
                   colour=cancerType))+
    
    scale_x_continuous(breaks=c(1,2,3,5,1-log10(capThresh)), labels=c('1', '2', '3', '5', '>10'))+
    scale_y_continuous(breaks=c(1,2,3,5,1-log10(capThresh)), labels=c('1', '2', '3', '5', '>10'))+
    
    #BELLS and whistles for the plot
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= maxCoord), colour='black', linetype=2)+
    
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= maxCoord), colour='black')+
    ylab('1 minus log(q) value in hypermutators')+
    xlab('1 minus log(q) value in non-hypermutators')+
    emptyTheme+
    labs(caption='runDnDsCv.R\nmake_dnds_figure_2c.ipynb')+
    ggtitle(title)
}

#current file lives @ '/Users/friedman/Desktop/WORK/dataForLocalPlotting/dndscvSummary.tsv'
figure2aDataFrame <- read.table(paste(plottingFilePath, 'figure2aDNDSSummary.tsv', sep='') , sep = '\t', header=TRUE)
p <- make_dnds_plot(figure2aDataFrame, 'DNDS-CV Comparing Hypermutated and\n Non-Hypermutated Cancers')
ggsave(paste(plottingFilePath, 'figure2a.pdf', sep=''), plot=p,  width = 6, height = 6)

#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2B


#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2C


#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2D



#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2E
