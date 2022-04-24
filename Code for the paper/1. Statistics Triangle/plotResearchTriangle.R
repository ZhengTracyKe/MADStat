library("ggplot2")
library("ggrepel")
library("deldir")
library("ggvoronoi")
load("output/Mixed-SCORE-results.RData")



############# plot rows of R ##############################
R <- MMres.summary$R
id_keepNodes <- which(MMres.summary$label_kmeans!=Inf)
plotdata.Rt <- data.frame(R1=R[id_keepNodes,1],R2=R[id_keepNodes,2])
p <- ggplot(data = plotdata.Rt, aes(x = R1, y = R2)) + theme_bw() + 
  geom_point(shape = 16, colour = "gray89", size = 1,stroke = 1) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-6, 2.4) + ylim(-4,5)
p

############# create the voronoi line segments#############
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


############# plot the high-degree authors###############
selectAuthors <- MMres.summary$selectAuthors
L <- nrow(selectAuthors)
label_point <- c()
label_index <- c()
for (i in 1:L) {
	label_point[[i]] <- paste(selectAuthors[i,],collapse = "\n")
	label_index[[i]] <- match(selectAuthors[i,], MMres.summary$authorNames)
}
plotdata.select <- data.frame(R1 = R[unlist(label_index),1],R2 = R[unlist(label_index),2])
set.seed(0)
label_size <- c(rep(4.5,L))
p <- p + geom_point(data = plotdata.select,
                     shape = 16,
                     colour = c(rep("blue",nrow(plotdata.select))),
                     size = 1,stroke = 1) + 
  geom_text_repel(data = plotdata.centers,
                  aes(label = label_point),
                  segment.size = 0.8,
                  colour = c(rep("red",nrow(plotdata.centers))),
                  size = label_size, 
                  box.padding = 0.8)

p


############# add cluster names ###############
cluster_names <- MMres.summary$clusterNames
cluster_names <- paste(as.character(1:length(cluster_names)),". ", cluster_names,sep = "")
cluster_names_box <- c()
cluster_names_box <- paste(cluster_names,collapse = "\n")


p <- p + geom_label_repel(data = data.frame(R1 = -6,R2 = 3),
                            aes(label = cluster_names_box),
                            size = 4.5, 
                            color = "blue",
                            box.padding = 0.1) 
                            
p


############## plot the triangle #######################
vertices <- MMres.summary$vertices
vertex_label <- c("Non-parametric", "Bayes", "Biostatistics")
index_adjust <- c(2,3,1)  
vertices <- vertices[index_adjust,]
plotdata.simplex <- data.frame(R1 = vertices[,1], R2 = vertices[,2])
p <- p + geom_label_repel(data = plotdata.simplex,
                            aes(label = vertex_label),
                            arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                            size = 9,
                            segment.size = 1,
                            color = "green4",
                            box.padding = 1,
                            nudge_y = c(-2,0.5,-1.5),
                            nudge_x = c(-0,0.5,-0.5))

plotdata.simplex2 <- rbind(plotdata.simplex, vertices[1,])
p <- p + geom_path(data = plotdata.simplex2,aes(x=R1,y=R2),linetype = 2, color = "green4", size = 0.8)
p



dev.copy2pdf(file=paste("ResearchTriangle.pdf",sep = ''),width = 12,height = 10)





