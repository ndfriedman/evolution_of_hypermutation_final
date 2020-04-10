#written by Noah Friedman
#scripts for plotting supplementary figures for figure 1

#
###
######
#############
##################
##########################
#################
############
#######
###
#


#SUPPLEMENTARY FIGURES

#TEMP

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/observedVsExpectedByType.tsv', sep='\t', header=TRUE)
s=1
p <- ggplot(df[(df$nmut > 5),], aes(x=nmut))+
  geom_smooth(aes(y = nmutNormal/dif, colour='Non-hypermutated driver genes'), method='loess', span=s)+
  geom_smooth(aes(y = nmutNormalAndHyperStrong/dif, colour='Non-hypermutated driver genes\nand strong hypermutant drivers'), method='loess', span=s)+
  geom_smooth(aes(y = nmutNormalAndHyperAll/dif, colour='Non-hypermutated driver genes\nand strong hypermutant drivers\nand weak hypermutant drivers'), method='loess', span=s)+
  scale_colour_manual(values=c('black', 'red', 'orange'))+
  ylim(0,1)+
  emptyTheme+
  xlab('Nmut in IM-341')+
  ylab('Fraction of divergence between\nobserved and expected\nexplained by:')+
  labs(caption='plotFigure1.R\nmake_obs_vs_expected_by_gene_type.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))



df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/giniData.tsv', sep='\t', header=TRUE)
ggplot(df, aes(x=iter))+
  geom_path(aes(y=oncogeneGini), colour='red')+
  geom_path(aes(y=tsgGini), colour='blue')+
  ylim(0,1)








#WORK on Monday Dec 2

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/essentialGenesMuts.tsv', ,sep = '\t', header=TRUE)

ggplot(df, aes(x= nmut))+
  
  #geom_smooth(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'), span=25)+
  #geom_smooth(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'), span=25)+
  #geom_smooth(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'), span=25)
  
  geom_point(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'))+
  geom_point(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'))+
  geom_point(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'))























p <- ggplot()+
  geom_smooth(data = df[(df$signature == 'MMR') & (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), method = 'lm', se=FALSE)+
  
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), alpha=0.35)+
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), alpha=0.35)+
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), alpha=0.35)+
  
  geom_smooth(data = df[(df$signature == 'POLE') & (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), method = 'lm', se=FALSE)+
  
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), alpha=0.35, bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), alpha=0.35, bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), alpha=0.35, bins=5)+
  
  xlab('TMB')+
  ylab('Difference between observed and\nexpected n truncating muts')+
  scale_color_manual(values=c('#013220', '#037D50', '#FF6347',
                              '#FF8C00', '#82CFFD', '#0D4F8B'))+
  labs(colour = 'Signature and\ngene type')+
  emptyTheme+
  scale_x_log10()

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))
