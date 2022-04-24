####Read the networks and Mixed-SCORE results###############
library(R.matlab)
Data <- readMat('output/CiteeDynamicFinal.mat')
NetworkSeries <- Data$CiteeAdjAggregate
authorNames <- as.matrix(read.table('data/author_name.txt', header=F, sep="\n"))
authorNames <- authorNames[Data$keepNodeID]
rm(list=c('Data'))
load('data/Mixed-SCORE-results.RData')


#### Obtain the trajectories of a few authors ###############
selectNames <- c("Peter Bickel", "Raymond Carroll", "Jianqing Fan",
                 "Peter Hall", "Xihong Lin", "Danyu Y. Lin", "Jun 1 Liu",
                 "Xiao-li Meng", "Donald B. Rubin", "Robert Tibshirani",
                 "Larry Wasserman", "Zhiliang Ying", "Bin Yu")
selectID <- match(selectNames, authorNames)
selectNames[7] <- "Jun Liu"  ##change "Jun 1 Liu" to "Jun Liu"

eigVecs <- MMres.summary$eigVecs
eigVals <- MMres.summary$eigVals
m <- length(selectID)
Tn <- length(NetworkSeries)
trajdata <- matrix(nrow=m*Tn, ncol=3)
namedata <- vector("character", length=m*Tn)
for (t in 1:Tn) {
    ind_fill <- t+ (0:(m-1))*Tn
	At <- as.matrix(NetworkSeries[[t]][[1]])
	Y <- At[selectID,] %*% eigVecs %*% diag(1/eigVals)
	tempPos <- Y[,2:3]/cbind(Y[,1],Y[,1])
	trajdata[ind_fill,] <- cbind(rep(t,m),tempPos)
	namedata[ind_fill] <- selectNames
}
traj_select <- data.frame(name=namedata,window=trajdata[,1],posX=trajdata[,2],posY=trajdata[,3])
save('traj_select', file='trajectories_select.RData')





#### Generate the trajectory plot #####################
library("ggplot2")
library("ggrepel")
library("deldir")
library("ggvoronoi")

### plot the background ##########
Rt <- MMres.summary$R[MMres.summary$label_kmeans!=Inf,]
plotdata.Rt <- data.frame(R1=Rt[,1],R2=Rt[,2])
p <- ggplot(data = plotdata.Rt, aes(x = R1, y = R2)) + theme_bw() + 
  geom_point(shape = 16, colour = "gray89", size = 1,stroke = 1) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-6, 2.4) + ylim(-4,5)
#### create voronoi line segments  #############
center_kmeans <-MMres.summary$center_kmeans
plotdata.centers <- data.frame(R1=center_kmeans[,1], R2=center_kmeans[,2])
voronoi <- deldir(plotdata.centers)
p <- p + geom_point(data = plotdata.centers, shape = 4,colour = "blue", size = 1.2,stroke = 1.6) +
  geom_text(data = plotdata.centers, label = 1:nrow(plotdata.centers),
            hjust = 0, vjust = 0,color = 'blue',cex = 10) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    size = 0.6,
    data = voronoi$dirsgs,
    linetype = 2,
    color= "blue")
p


### plot the trajectories #####
col.list <-  c('gold2','green','blue','cyan3','red','pink2','black',
'purple','darkgreen', 'gold4','magenta2','peru','slateblue')
pch.list <- c(2:(length(selectNames) + 1))

plot.trajdata <- data.frame(R1=traj_select$posX, R2=traj_select$posY)
plot.labels <- as.character(traj_select$window) 
plot.color <- as.vector(matrix(rep(col.list,Tn), nrow=Tn, byrow=T))
plot.pch <- as.vector(matrix(rep(pch.list,Tn), nrow=Tn, byrow=T))

##plot positions#####
p1 <- p + geom_point(data = plot.trajdata, shape = plot.pch,
                     aes(x = R1,y = R2),
                     colour = plot.color, size = 0.4,stroke = 2)

##add labels 1,11,21####
tempind <- which(traj_select$window %in% c(1,11,21))
p1 <- p1 + geom_text_repel(data = plot.trajdata[tempind,],                          
                           aes(label=plot.labels[tempind]), 
                           size = 4, colour=plot.color[tempind])

##plot the paths######
plot.path <- data.frame(R1=traj_select$posX, R2=traj_select$posY, name=traj_select$name)
p2 <- p1 + geom_path(data = plot.path, aes(x=R1,y=R2, group = name, color=name),      
                     linetype = 1, size = 0.8) 
### In the above plot, "color" specifies both the legend and color
### Need to change the actual color
names(col.list) <- selectNames
names(pch.list) <- selectNames
p2 <- p2 + scale_color_manual(name = '', values = col.list) + scale_shape_manual("", values = pch.list) + scale_linetype_manual("", values = 1) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),legend.text=element_text(size=15),
        legend.background = NULL)
p2





##restrict to a rectangle ###
p3 <- p2 + coord_cartesian(xlim = c(-2.5,1),ylim = c(-1.2,2.5))
p3
dev.copy2pdf(file=paste("trajectory_all.pdf",sep = ''), width=9, height = 7.5)





