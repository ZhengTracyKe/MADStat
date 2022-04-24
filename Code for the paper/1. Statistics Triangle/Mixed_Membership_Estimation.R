####Read the adjacency matrix and author names################
library(R.matlab)
Data <- readMat('output/CiteeAdjFinal.mat')
A <- as.matrix(Data$A)
authorID <- Data$keepNodeID
rm('Data')
authorNames <- as.matrix(read.table('data/author_name.txt', header=F, sep="\n"))
authorNames <- authorNames[authorID]



####Estimating mixed memberships by Mixed-SCORE################

##SCORE embedding#####
K <- 3
eig.out <- RSpectra::eigs(A, k = K)
eigVecs <- as.matrix(eig.out$vectors[,1:K])
eigVals <- eig.out$values
R <- eigVecs[,2:K] / kronecker(matrix(1,1,K-1),eigVecs[,1])

##Clustering#########
n <- dim(R)[1]
r <- apply(R,MARGIN = 2,FUN = median)
norms <- rowSums(abs(R - matrix(rep(1,dim(R)[1]),ncol = 1) %*% r))
remove_percent <- 0.03   ##remove outliers
id_keepNodes <- which(norms <= quantile(norms,1 - remove_percent)) 
rm(".Random.seed") 
set.seed(100)
L <- 15
kmeansRes <- kmeans(R[id_keepNodes,], L, iter.max=400, nstart=800)
center_kmeans <- kmeansRes$centers
label_kmeans <- vector("numeric", length=n)
label_kmeans[id_keepNodes] <- kmeansRes$cluster
label_kmeans[-id_keepNodes] <- Inf   ##for outliers, cluster labels are Inf
 
##Vertex hunting#####
source('mixSCORE.R')
index.matrix <- combn(L, K)  # K * (L K)  
dist <- sapply(1:ncol(index.matrix), function(i){
    getMaxDist(centers = center_kmeans, 
               vertex.ind = index.matrix[,i])
})
tempind <- index.matrix[, which.min(dist)]
vertices <- center_kmeans[tempind,]

##Membership estimation#####
source('mixSCORE.R')
memb.out <- getMembership(R = R, 
                         vertices = vertices, 
                         K = K,
                         eig.values = eigVals,
                         eig.vectors = eigVecs)
memberships <- memb.out$memberships





####Re-ordering the clusters (for later plotting)#################
degrees <- rowSums(A)
m <- 5   ##select 5 authors per cluster
selectAuthors <- matrix(nrow=L, ncol=m)
for (i in 1:L) {
	tempdegrees <- degrees[label_kmeans==i]
	tempauthornames <- authorNames[label_kmeans==i]
	selectind <- sort(tempdegrees, decreasing=T, index.return=T)$ix[1:m]
	selectAuthors[i,] <- tempauthornames[selectind]
}

ref_authors <- c("Theo Gasser", "Peter Hall", "David Donoho", "Michael Martin",
                 "Tze Leung Lai", "Ross Prentice", "Raymond Carroll", "Peter Green",
                 "Nicholas Polson", "Robert Kass", "Adrian Smith", "Scott Zeger",
                 "Lueping Zhao", "John Oquigley", "Michael Proschan")
permute <- match(ref_authors, selectAuthors[,1])
center_kmeans <- center_kmeans[permute,]
selectAuthors <- selectAuthors[permute,]
oldlabel <- label_kmeans
for (k in 1:L) {
	label_kmeans[which(oldlabel==permute[k])] <- k
}
cluster_names <- c("Kernel Estimation", "Non-parametric I", "Wavelets", 
                   "Decision Theory", "Survival Analysis I", "Survival Analysis II", 
                   "Non-parametric II", "Bayes I(MCMC)", "Bayes II",
                   "Bayes III", "Bayes IV(GLMM)", "Longitudinal I(GEE)",
                   "Longitudinal II","Epidemiology","Clinical Trial")



####Save the results#######################################
MMres.summary <- list(eigVals=eigVals, eigVecs=eigVecs, authorNames=authorNames, R=R, b1=memb.out$b1, vertices=vertices, center_kmeans=center_kmeans, label_kmeans=label_kmeans, clusterNames= cluster_names, selectAuthors=selectAuthors)

save('MMres.summary', file='Mixed-SCORE-results.RData')


