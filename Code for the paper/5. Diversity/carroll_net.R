##### transform accent letters to english letters, need to be opened with encoding utf-8 to show normally
to_english <- function(s) {
  # 1 character substitutions
  old1 <- "ŚŠŞşśšŽŻŻŻżžźþÅȦàáăặắâãäąåḃČĆÇçćčĐèééėễėêëęȩİİìıíĭîïķľłŁðñņńØÖòóôõöőùůūúûůüÜýÿğĝǧḡeřţt"
  new1 <- "SSSsssZZZZzzzyAAaaaaaaaaaabCCCcccDeeeeeeeeeeIIiiiiiikllLdnnnOOoooooouuuuuuuUyyggggertt"
  s1 <- chartr(old1, new1, s)
  
  # 2 character substitutions
  old2 <- c("œ", "ß", "æ", "ø","a̧","g̃","t́")
  new2 <- c("oe", "ss", "ae", "oe","a","g","t")
  s2 <- s1
  for(i in seq_along(old2)) s2 <- gsub(old2[i], new2[i], s2, fixed = TRUE)
  
  return(s2)
}

##### some functions for network testing
tr=function(A){
  a <- sum(diag(A))
  return(a)
}
test0<-function(id_interest,A_all){
  #select the neighbor of id_interest
  id_neighbor <- which(A_all[,id_interest] > 0)
  
  #build the coauthorship network adjacency matrix
  sz <- length(id_neighbor)
  A <- as.matrix(A_all[id_neighbor,id_neighbor])
  
  diag(A) <- 0
  A[which(A > 0)] <- 1
  
  #compute B/T/Q test statistics
  c <- rep(1,sz)
  V <- as.numeric(sqrt(t(c) %*% A %*% c))
  eta <- as.vector((A %*% c) / V)
  ah <- A - eta %*% t(eta)
  i <- c(1:sz)
  j <- c(1:sz)
  x <- diag(ah)
  diagah <- sparseMatrix(i,j,x=x)
  
  Q <- as.numeric(tr(ah %*% ah %*% ah %*% ah) - 4 * tr(ah * (ah %*% ah %*% ah)) + 8 * tr(ah * ah * (ah %*% ah)) - 
                 6 * tr(ah * ah * ah * ah) - 2 * tr((ah %*% ah) * (ah %*% ah)) + 
                 2 * t(c) %*% (diagah %*% (ah * ah) %*% diagah) %*% c + t(c) %*% (ah * ah * ah * ah) %*% c)
  
  m <- as.numeric(t(eta) %*% eta)
  
  Qtest <- (Q - 2 * ((m - 1)^2)) / sqrt(8 * (m - 1)^4)
  
  #compute p-value
  pqvalue <- as.numeric(1-pnorm(Qtest))
  
  ps <- list()
  ps$pq <- pqvalue
  ps$size <- sz
  
  return(ps)
}


#### load data ####
library('R.matlab')
library(igraph)

rm(list=ls())
setwd("E:/Dropbox/CMU/ADA/Data/data2019/")

community <- readMat(paste("CoauAdjs.mat",sep = ''),
                     sparseMatrixClass=c("Matrix"))

A_cumu <- community$adjCoauCumu
A_all <- A_cumu[[41]][[1]]


######### Raymond Carroll ############

name <- unlist(community$name)
name_interest <- 'Raymond Carroll'
##### authors who cited name_interest
id_neighbor<-which(A_all[,which(name==name_interest)]>0)
A <- as.matrix(A_all[id_neighbor,id_neighbor])
A_orig <- A
name_neighbor <- name[id_neighbor]
"Raymond Carroll" %in% name_neighbor
dim(A)

diag(A) <- 0
A[which(A>0)] <- 1
test0(which(name == name_interest), A_all)

##### first apply community detection to it #####
source('mixSCORE.R')
library(RSpectra)

eig.out <- eigs(A, k = 10)
eig.values <- eig.out$values
eig.values <- sort(abs(eig.values),decreasing = TRUE)
plot(1:10,eig.values)

community_detection <- function(A){
  id_interest <- which(name_neighbor == name_interest)
  id_rest = setdiff(c(1:dim(A)[1]), id_interest)
  
  A = A[id_rest, id_rest]
  g <- graph_from_adjacency_matrix(A, mode = 'undirected')
  cl <- components(g)
  cl$csize
  # cl$membership
  
  gc_id <- which(cl$membership == 1)
  A <- A[gc_id,gc_id]
  A_orig <- A_orig[gc_id,gc_id]
  name_neighbor <- name_neighbor[gc_id]
  ##### first apply community detection to it
  source('mixSCORE.R')
  library(RSpectra)
  help(eigs)
  eig.out <- RSpectra::eigs(A, k = 10)
  eig.values <- eig.out$values
  eig.values <- sort(abs(eig.values),decreasing = TRUE)
  plot(1:10,eig.values[1:10])
  
  set.seed(0)
  result <- SCORE(A, K = 2)
  table(result$labels)
  label = rep(1, length(id_rest) + 1)
  label[id_rest][gc_id][which(result$labels == 2)] = 2
  label[id_interest] = 2
  return(label)
}

label = community_detection(A)
name_neighbor[label == 1]
col.list <- c("blue","black")[label]


##### draw the graph of the network #####
# distribution of node degrees
summary(colSums(A))
summary(rowSums(A))

deg_threshold <- 40
id_deg <- which(colSums(A_orig) >= deg_threshold)
c <- name[id_neighbor[id_deg]]
m <- length(c)

col.list <- col.list[id_deg]
#store the weight
A_1 <- as.matrix(A_all[id_neighbor,id_neighbor][id_deg,id_deg])
diag(A_1) <- 0
l <- rep(0,m)
for (i in 1:m) {
  l[i]=as.numeric(test0(which(name==c[i]),A_all)$pq)
}

c <- iconv(c, from = "", to = "utf-8")
c <- iconv(c, from = "utf-8", to = "")
c <- to_english(c)

draw_graph <- function(seed){
  set.seed(seed)
  permut <- sample(1:m, m, replace = FALSE)
  A_g <- A_1[permut,permut]
  col.list <- col.list[permut]
  deg_A_g <- colSums(A_g)
  l_g <- l[permut]
  name_g <- c[permut]
  str<-graph_from_adjacency_matrix(A_g,mode="undirected",weighted = TRUE)
  
  summary(colSums(A_1))
  summary(rowSums(A_1))
  #select edges by weight
  str_trim=delete.edges(str,which(E(str)$weight<1))
  node_label <- paste(name_g,paste(signif(l_g,digits=2)),sep = "\n")
  node_label[deg_A_g < quantile(deg_A_g,0.4)] <- ""
  a <- colSums(A_g)^0.3
  label_size <- 2 / max(a) * a
  
  g <- str_trim
  g_layout <- layout_on_sphere(g)
  par(mar = c(1,1,1,1))
  label_distance <- rep(0.4,m)
  label_position <- rep(0,m)
  
  p <- plot(g,
            edge.arrow.size = .1, 
            layout = g_layout,
            vertex.color = "red", 
            vertex.size = 3, 
            vertex.frame.color = "gray", 
            vertex.label.color= col.list, 
            vertex.label.cex = label_size,
            vertex.label = node_label,
            vertex.label.font = 2,
            vertex.label.family = "Helvetica",
            vertex.label.dist = label_distance,
            vertex.label.degree = label_position)  
}

draw_graph(160)



#### adjust labels ####
set.seed(160)
permut <- sample(1:m, m, replace = FALSE)
A_g <- A_1[permut,permut]
col.list <- col.list[permut]
deg_A_g <- colSums(A_g)
l_g <- l[permut]
name_g <- c[permut]
str<-graph_from_adjacency_matrix(A_g,mode="undirected",weighted = TRUE)

summary(colSums(A_1))
summary(rowSums(A_1))
#select edges by weight
str_trim=delete.edges(str,which(E(str)$weight<1))
node_label <- paste(name_g,paste(signif(l_g,digits=2)),sep = "\n")
node_label[deg_A_g < quantile(deg_A_g,0.4)] <- ""
a <- colSums(A_g)^0.29
label_size <- 2 / max(a) * a

g <- str_trim
g_layout <- layout_on_sphere(g)
par(mar = c(1,1,1,1))
label_distance <- rep(0.4,m)
label_position <- rep(0,m)

name_need_adjust <- c("Matt P. Wand","Wolfgang Hardle","Laurence Freedman","Nilanjan Chatterjee","Bhramar Mukherjee")
id_need_adjust <- c()
for (i in 1:length(name_need_adjust)){
  id_need_adjust <- append(id_need_adjust,which(name_g == name_need_adjust[i]))
}
adjust_position <- c(pi/2,pi,0,pi,pi * 7/8)
adjust_distance <- c(2,2,2,1.5,2)
label_position[id_need_adjust] <- adjust_position
label_distance[id_need_adjust] <- adjust_distance

p <- plot(g,
          edge.arrow.size = .1, 
          layout = g_layout,
          vertex.color = "red", 
          vertex.size = 3, 
          vertex.frame.color = "gray", 
          vertex.label.color= col.list, 
          vertex.label.cex = label_size,
          vertex.label = node_label,
          vertex.label.font = 2,
          vertex.label.family = "Helvetica",
          vertex.label.dist = label_distance,
          vertex.label.degree = label_position)  

dev.copy2pdf(file=paste("net-Carroll.pdf",sep = ''))