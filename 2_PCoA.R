rm(list=ls())
library(vegan)
library(ggplot2)
library(ggrepel)

# load data
df <- read.delim('PCOA_bacteria.txt', row.names = 1, sep = '\t', head = TRUE, check.names = FALSE)
group <- read.delim('group_pcoa.txt', row.names = 1, sep = '\t', head = TRUE, check.names = FALSE)

df1 <- t(df)

distance <- vegdist(df1, method='bray')
pcoa <- cmdscale(distance, k=(nrow(df1)-1), eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]

# rename 
names(plot_data)[1:2] <- c('PCoA1','PCoA2') 
eig <- pcoa$eig
group1 <- group['group']
data <- plot_data[match(rownames(group), rownames(plot_data)),]
data <- data.frame(group, plot_data)

# plot
ggplot(data,
	aes(x=PCoA1,y=PCoA2,shape=group,color=group)) + \
    geom_point(alpha=1,size=8)+stat_ellipse(level=0.95,size=3) + \
    labs(
    	x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),
    	y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep="")) + \
    geom_vline(aes(xintercept=0),linetype="dotted") + \
    geom_hline(aes(yintercept=0),linetype="dotted") + \
    theme(
    	panel.background = element_rect(fill='white',colour = 'black'),
    	axis.title.x=element_text(colour = 'black',size=20),
    	axis.title.y = element_text(colour = 'black',size=20),
    	legend.text = element_text(size = 15)
    	)
