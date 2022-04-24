#' mixedSCORE
#' 
#' @param A n-by-n binary symmtric adjacency matrix.
#' @param K number of communities.
#' @param verbose whether generate message
#'
#' @return A list containing \describe{
#'   \item{R}{n-by-(K-1) ratio matrix.}
#'   \item{L}{Selected tunning parameter used for vertex hunting algorithm.}
#'   \item{thetas}{A vector of the estimated degree heterogeniety parameters}
#'   \item{vertices}{K-by-(K-1) K vertices of the found convex hull}
#'   \item{centers}{L-by-(K-1) L centers by kmeans}
#'   \item{memberships}{n-by-K membership matrix.}
#'   \item{purity}{A vector of maximum membership of each node}
#'   \item{hard.cluster.labels}{A vector of integers indicating hard clutering labels, by assigning the node to the cluster with max membership}
#' }
#' 
#' @import RSpectra, combinat, limSolve
mixedSCORE <- function(A, K, threshold = NULL, verbose = F, remove_outliers = 0.01,L.candidate = (K + 1):(3 * K)){
  
  # check package RSpectra installed, otherwise install it 
  if (!'RSpectra' %in% installed.packages()){
    install.packages('RSpectra', quiet = T)
    if (!'RSpectra' %in% installed.packages()){
      stop("Package: 'RSpectra' not found!")
    }
  }
  
  # check package combniat installed, otherwise install it 
  if (!'combinat' %in% installed.packages()){
    install.packages('combinat', quiet = T)
    if (!'combinat' %in% installed.packages()){
      stop("Package: 'combinat' not found!")
    }
  }
  
  # check package combniat installed, otherwise install it 
  if (!'limSolve' %in% installed.packages()){
    install.packages('limSolve', quiet = T)
    if (!'limSolve' %in% installed.packages()){
      stop("Package: 'limSolve' not found!")
    }
  }
  
  
  if(verbose){
    cat('Get ratios --------\n')
  }
  score.out = SCORE(A = A, K = K, threshold = threshold)
  R = score.out$R
  if(verbose){
    cat('Vertex hunting --------\n')
  }
  vh.out = vertexHunting(R = R, K = K,verbose = verbose,remove_outliers = remove_outliers,L.candidate = L.candidate)
  if(verbose){
    cat('Get the membership ----\n')
  }
  memb.out = getMembership(R = R, 
                           vertices = vh.out$vertices, 
                           K = K,
                           eig.values = score.out$eig.values,
                           eig.vectors = score.out$eig.vectors)
  memberships = memb.out$memberships
  degrees = memb.out$degrees
  b1 <- memb.out$b1
  
  if(verbose){
    cat('Get the purity scores and hard clustering results ----\n')
  }
  purity = apply(X = memberships, MARGIN = 1, FUN = max)
  major.labels = apply(X = memberships, MARGIN = 1, FUN = which.max)
  
  return(list(R = R,
              L = vh.out$L,
              vertices = vh.out$vertices,
              centers = vh.out$centers,
              memberships = memberships,
              degrees = degrees,
              puritys = purity,
              major.labels = major.labels,
              b1 = b1,
              eig.values = score.out$eig.values,
              eig.vectors = score.out$eig.vectors))
}



#' Spectral Clustering On Ratios-of-Eigenvectors (SCORE)
#' 
#' @param A n-by-n binary symmtric adjacency matrix.
#' @param K number of communities.
#' @param threshold (optional) the threshold of ratio matrix. By defalt is \code{log(n)}.
#'
#' @return A list containing \describe{
#'   \item{R}{n-by-(K-1) ratio matrix.}
#'   \item{labels}{A vector of integer indicating the cluster to which each point allocated.}
#' }
#' 
#' @export
SCORE <- function(A, K, threshold = NULL){
  
  # check A is matrix
  if (!is.matrix(A)){
    A = as.matrix(A)
    if (any( A < 0)){
      stop('Entries of adjacency matrix A should be nonegative!')
    }
  }
  
  # check symetry
  if(any(A != t(A))){
    stop("Aadjacency matrix A is not symmetric!")
  }
  
  # check connectivity of the network
  if (!'igraph' %in% installed.packages()){
    install.packages('igraph', quiet = T)
  }
  if (!'igraph' %in% installed.packages()){
    warnings('igraph package not found. Skip connectivity check.')
  } else {
    if(!igraph::is.connected(igraph::graph_from_adjacency_matrix(A, mode = 'undirected'))){
      warnings("Network disconnected!")
    }
  } 
  
  
  eig.out = RSpectra::eigs(A, k = K)
  ev = eig.out$vectors[,1:K]
  eig.values = eig.out$values
  R = ev[,2:K] / matrix(rep(ev[,1], K-1), ncol = K-1, byrow = F)
  R = as.matrix(R, nrow = nrow(ev))
  
  if (is.null(threshold)){
    threshold =  Inf #log(nrow(A)) 
  }
  
  
  if (min(eig.values[1:K]) < 0){
    min.eig = min(eig.values[1:K])
    threshold = sqrt(abs(eig.values[1]) / abs(min.eig))
  }
  
  
  R[R > threshold] = threshold
  R[R < -threshold] = -threshold
  
  labels.hat = kmeans(R, K, nstart = 100, iter.max = 100)$cluster
  
  return(list(R = R, 
              labels = labels.hat,
              eig.values = eig.out$values[1:K],
              eig.vectors = ev))
}



draw_embedding <- function(A,K,name,d_min,d_max,remove_outliers){
  # check A is matrix
  if (!is.matrix(A)){
    A = as.matrix(A)
    if (any( A < 0)){
      stop('Entries of adjacency matrix A should be nonegative!')
    }
  }
  
  # check symetry
  if(any(A != t(A))){
    stop("Aadjacency matrix A is not symmetric!")
  }
  
  # check connectivity of the network
  if (!'igraph' %in% installed.packages()){
    install.packages('igraph', quiet = T)
  }
  if (!'igraph' %in% installed.packages()){
    warnings('igraph package not found. Skip connectivity check.')
  } else {
    if(!igraph::is.connected(igraph::graph_from_adjacency_matrix(A, mode = 'undirected'))){
      warnings("Network disconnected!")
    }
  } 
  
  
  d <- rowSums(A)
  valid_ix <- which(d >= d_min & d <= d_max)
  A <- A[valid_ix,valid_ix]
  name <- name[valid_ix]
  
  g <- igraph::graph_from_adjacency_matrix(A, mode = 'undirected')
  co <- igraph::components(g, mode = 'STRONG')
  gi <- igraph::induced.subgraph(g, which(co$membership == which.max(co$csize)))
  name <- name[which(co$membership == which.max(co$csize))]
  # if you want here you can decide if you want values only
  # in the upper or lower triangle or both
  ad <- igraph::get.adjacency(gi)
  A <- as.matrix(ad)
  
  R <- get_embedding(A,K)$R
  
  r <- apply(R,MARGIN = 2,FUN = median)
  norms <- rowSums(abs(R - matrix(rep(1,dim(R)[1]),ncol = 1) %*% r))
  
  id_outlier <- which(norms >= quantile(norms,1 - remove_outliers))
  R <- R[-id_outlier,]
  A <- A[-id_outlier,-id_outlier]
  name <- name[-id_outlier]
  
  
  if(K == 3){
    plot(R, col='grey', lwd = 2, xlab = 'R[,1]', ylab = 'R[.2]',bty="n")
    # lines(ms.out$vertices[c(1,2,3,1),1], ms.out$vertices[c(1,2,3,1),2], 
    #       lty = 2, lwd = 2, col = 'black')
    # points(ms.out$centers, lwd = 2, col = 'blue')
    
    d <- rowSums(A)
    # sort(d,'decreasing' = TRUE,index.return = TRUE)$ix[1:10]
    
    sort_id <- sort(d,'decreasing' = TRUE,index.return = TRUE)$ix
    name_sorted = name[sort_id]
    m <- min(60,dim(R)[1])
    points(R[sort_id[1:m],1], R[sort_id[1:m],2], lwd = 2, col = 'darkorange')
    m <- min(30,dim(R)[1])
    text(R[sort_id[1:m],1], R[sort_id[1:m],2],  
         name_sorted[1:m], cex=0.5, pos=3,col="red") 
    points(R[sort_id[1:m],1], R[sort_id[1:m],2], lwd = 2, col = 'red')
  }
  if (K == 4){
    rgl::plot3d(R, col='grey', lwd = 2,
        xlab = 'R[,1]', ylab = 'R[,2]',zlab = 'R[,3]')
    # lines(ms.out$vertices[c(1,2,3,1),1], ms.out$vertices[c(1,2,3,1),2], 
    #       lty = 2, lwd = 2, col = 'black')
    # points(ms.out$centers, lwd = 2, col = 'blue')
    
    d <- rowSums(A)
    # sort(d,'decreasing' = TRUE,index.return = TRUE)$ix[1:10]
    
    sort_id <- sort(d,'decreasing' = TRUE,index.return = TRUE)$ix
    name_sorted = name[sort_id]
    
    m <- min(60,dim(R)[1])
    rgl::points3d(R[sort_id[1:m],], lwd = 2, col = 'darkorange',bty="n")
    
    m <- min(30,dim(R)[1])
    rgl::text3d(R[sort_id[1:m],],
           texts = name_sorted[1:m], cex=0.65, pos=3,col="red")
    rgl::points3d(R[sort_id[1:m],], lwd = 2, col = 'red',bty="n")
  }
    
  return(list(R = R,name = name,A = A))
}



get_embedding <- function(A, K, threshold = NULL){
  
  # check A is matrix
  if (!is.matrix(A)){
    A = as.matrix(A)
    if (any( A < 0)){
      stop('Entries of adjacency matrix A should be nonegative!')
    }
  }
  
  # check symetry
  if(any(A != t(A))){
    stop("Aadjacency matrix A is not symmetric!")
  }
  
  # check connectivity of the network
  if (!'igraph' %in% installed.packages()){
    install.packages('igraph', quiet = T)
  }
  if (!'igraph' %in% installed.packages()){
    warnings('igraph package not found. Skip connectivity check.')
  } else {
    if(!igraph::is.connected(igraph::graph_from_adjacency_matrix(A, mode = 'undirected'))){
      warnings("Network disconnected!")
    }
  } 
  
  
  
  
  
  eig.out = RSpectra::eigs(A, k = K)
  ev = eig.out$vectors[,1:K]
  eig.values = eig.out$values
  R = ev[,2:K] / matrix(rep(ev[,1], K-1), ncol = K-1, byrow = F)
  R = as.matrix(R, nrow = nrow(ev))
  
  # if (is.null(threshold)){
  #   threshold =  Inf #log(nrow(A)) 
  # }
  # 
  # 
  # if (min(eig.values[1:K]) < 0){
  #   min.eig = min(eig.values[1:K])
  #   threshold = sqrt(abs(eig.values[1]) / abs(min.eig))
  # }
  
  
  # R[R > threshold] = threshold
  # R[R < -threshold] = -threshold
  
  return(list(R = R, 
              eig.values = eig.out$values[1:K],
              eig.vectors = ev))
}



#' Vertex hunting algorithm to find the cluster centers
#' 
#' @param R n-by-(K-1) ratio matrix
#' @param K number of communities.
#' 
#' @return A list containing \describe{
#'   \item{L}{Selected tunning parameter.}
#'   \item{vertices}{{K-by-(K-1) cluster center matrix shows the K community centers.}
#' }
#' 
#' @export
vertexHunting <- function(R, K, verbose = F, remove_outliers, L.candidate = (K + 1):(3 * K)){
  
  r <- apply(R,MARGIN = 2,FUN = median)
  norms <- rowSums(abs(R - matrix(rep(1,dim(R)[1]),ncol = 1) %*% r))
  
  id_outliers <- which(norms > quantile(norms,1 - remove_outliers))
  if(length(id_outliers) > 0){
    R <- R[-id_outliers,]
  }
  
  # L.candidate = (K + 1):(3 * K) # candidate Ls
  n.L = length(L.candidate)
  
  out.list = lapply(L.candidate, function(L){
    if(verbose){
      cat('L = ',L,'\n')
    }
    centers = kmeans(R, L, iter.max = 100, nstart = 100)$centers # L*(K-1)
    out = vertexSearch(centers, K = K)
    if(verbose){
      cat( 'dist = ', out$dist,'\n')
    }
    return(out)
  })
  
  
  delta.L = sapply(1:n.L, function(i){
    if( i == 1){
      V.Lminus1 =  kmeans(R, K, iter.max = 100, nstart = 100)$centers
    } else {
      V.Lminus1 = out.list[[i-1]]$vertices
    }
    V.L = out.list[[i]]$vertices
    
    perm = combinat::permn(1:K) # list of all the permutations
    delta = min(sapply(1:length(perm), function(p){
      max(rowSums((V.L[perm[[p]],] - V.Lminus1)^2))
    }))
    return(delta/(1+out.list[[i]]$dist))
  })
  
  L.ind = which.min(delta.L)
  L.select = L.candidate[L.ind]
  
  if(verbose){
    cat('Select L =', L.select,'\n')
  }
  
  return(list(L = L.select,
              vertices = out.list[[L.ind]]$vertices,
              centers = out.list[[L.ind]]$centers))
}

#' select the \code{K} vertices from given \code{L} centers
#' 
#' @param centers L-by-(K-1) center matrix
#' @param K number of communities.
#' 
#' @return A list containing \describe{
#'   \item{ind}{a vector of \code{K} integers indicating the index of selected \code{K} vertices out of \code{L} centers.}
#'   \item{dist}{The maximum distance from centers to the convex hull formed by the \code{K} selected vertice}
#' }
#' 
vertexSearch <- function(centers, K){
  L = nrow(centers)
  
  index.matrix = combn(L, K)  # K * (L K) 
  n.c =  ncol(index.matrix )
  
  dist = sapply(1:n.c, function(i){
    getMaxDist(centers = centers, 
               vertex.ind = index.matrix[,i])
  })
  
  ind = index.matrix[, which.min(dist)[1]]
  vertices = matrix(centers[ind,],nrow = K)
  
  return(list(ind = ind, 
              dist = min(dist), 
              vertices = vertices,
              centers = centers))
}



reorder <- function(centers){
  centers[order(apply(centers, 1, function(x) norm(x,'2')),decreasing = T),]
}

#' find the maxinum distance from the convex hull formed by the chosen K vertices
#' 
#' @param centers L-by-(K-1) center matrix
#' @param vertex.ind index of the \code{K} centers that forms the convex hull
#' 
#' @return the maximum distance
#' 
getMaxDist <- function(centers, vertex.ind){
  L = nrow(centers)
  K = length(vertex.ind)
  
  nonvertex = matrix(centers[-vertex.ind,], nrow = L-K)
  vertex = matrix(centers[vertex.ind,], nrow = K)
  
  # A = t(vertex)
  # eig_A <- base::eigen(crossprod(A,A) + diag(rep(1,K) * 1e-08,nrow = K),symmetric = TRUE)$values
  # if (min(eig_A) <= 1e-08){
  #   return(+Inf)
  # }
  
  # else{
  dist = sapply(1:(L-K), function(i){
    node = nonvertex[i,]
    # calculate the distance from the convex hull
    # out = robust_lsei(A = t(vertex), 
    #                      B = node, 
    #                      E = rep(1,K), F = 1,
    #                      G = diag(K), H = rep(0,K),
    #                      type = 2)
    # calculate the distance from the convex hull
    out = limSolve::lsei(A = t(vertex), 
                         B = node, 
                         E = rep(1,K), F = 1,
                         G = diag(K), H = rep(0,K),
                         type = 2)
    return(out$solutionNorm)
  })
  return(max(dist))
  # }
}





robust_lsei <- function (A = NULL, B = NULL, E = NULL, F = NULL, G = NULL, H = NULL, 
          Wx = NULL, Wa = NULL, type = 1, tol = sqrt(.Machine$double.eps), 
          tolrank = NULL, fulloutput = FALSE, verbose = TRUE) {
  if (is.vector(E) & length(F) == 1) 
    E <- t(E)
  else if (!is.matrix(E) & !is.null(E)) 
    E <- as.matrix(E)
  if (is.vector(A) & length(B) == 1) 
    A <- t(A)
  else if (!is.matrix(A) & !is.null(A)) 
    A <- as.matrix(A)
  if (is.vector(G) & length(H) == 1) 
    G <- t(G)
  else if (!is.matrix(G) & !is.null(G)) 
    G <- as.matrix(G)
  if (!is.matrix(F) & !is.null(F)) 
    F <- as.matrix(F)
  if (!is.matrix(B) & !is.null(B)) 
    B <- as.matrix(B)
  if (!is.matrix(H) & !is.null(H)) 
    H <- as.matrix(H)
  if (is.null(A) && is.null(E)) {
    if (is.null(G)) 
      stop("cannot solve least squares problem - A, E AND G are NULL")
    A <- matrix(data = 0, nrow = 1, ncol = ncol(G))
    B <- 0
  }
  else if (is.null(A)) {
    A <- matrix(data = E[1, ], nrow = 1)
    B <- F[1]
  }
  Neq <- nrow(E)
  Napp <- nrow(A)
  Nx <- ncol(A)
  Nin <- nrow(G)
  if (is.null(Nx)) 
    Nx <- ncol(E)
  if (is.null(Nx)) 
    Nx <- ncol(G)
  if (is.null(Nin)) 
    Nin <- 1
  if (is.null(Neq)) {
    Neq <- 0
    if (verbose & type == 1) 
      warning("No equalities - setting type = 2")
    type = 2
    F <- NULL
  }
  else {
    if (ncol(E) != Nx) 
      stop("cannot solve least squares problem - A and E not compatible")
    if (length(F) != Neq) 
      stop("cannot solve least squares problem - E and F not compatible")
  }
  if (is.null(G)) 
    G <- matrix(data = 0, nrow = 1, ncol = Nx)
  if (is.null(H)) 
    H <- 0
  if (ncol(G) != Nx) 
    stop("cannot solve least squares problem - A and G not compatible")
  if (length(B) != Napp) 
    stop("cannot solve least squares problem - A and B not compatible")
  if (length(H) != Nin) 
    stop("cannot solve least squares problem - G and H not compatible")
  if (!is.null(Wa)) {
    if (length(Wa) != Napp) 
      stop("cannot solve least squares problem - Wa should have length = number of rows of A")
    A <- A * Wa
    B <- B * Wa
  }
  Tol <- tol
  if (is.null(Tol)) 
    Tol <- sqrt(.Machine$double.eps)
  IsError <- FALSE
  if (type == 1) {
    ineq <- Nin + Nx
    mIP <- ineq + 2 * Nx + 2
    lpr <- 1
    if (fulloutput) 
      lpr <- lpr + 3
    if (!is.null(tolrank)) 
      lpr <- lpr + 6
    if (!is.null(Wx)) {
      lw <- length(Wx)
      lpr <- lpr + 2 + lw
    }
    ProgOpt <- rep(1, lpr)
    if (lpr > 1) {
      ipr <- 1
      if (fulloutput) {
        ProgOpt[ipr:(ipr + 2)] <- c(ipr + 3, 1, 1)
        ipr <- ipr + 3
      }
      if (!is.null(tolrank)) {
        if (length(tolrank) == 1) 
          tolrank <- rep(tolrank, len = 2)
        ProgOpt[ipr:(ipr + 5)] <- c(ipr + 6, 4, tolrank[1], 
                                    ipr + 6, 5, tolrank[2])
        ipr <- ipr + 6
      }
      if (!is.null(Wx)) {
        lw <- length(Wx)
        if (lw == 1) {
          ProgOpt[ipr:(ipr + 2)] <- c(ipr + 3, 2, 1)
        }
        else {
          if (lw != Nx) 
            stop("cannot solve least squares problem - number of weighs should be =1 or =number of unknowns")
          lw <- lw + ipr + 1
          ProgOpt[ipr:lw] <- c(lw + 1, 3, Wx)
        }
      }
    }
    mdW <- Neq + Napp + ineq
    if (fulloutput) 
      mdW <- max(mdW, Nx)
    mWS <- 2 * (Neq + Nx) + max(Napp + ineq, Nx) + (ineq + 
                                                      2) * (Nx + 7)
    storage.mode(A) <- storage.mode(B) <- "double"
    storage.mode(E) <- storage.mode(F) <- "double"
    storage.mode(G) <- storage.mode(H) <- "double"
    sol <- .Fortran("lsei", NUnknowns = Nx, NEquations = Neq, 
                    NConstraints = Nin, NApproximate = Napp, A = A, B = B, 
                    E = E, F = F, G = G, H = H, X = as.vector(rep(0, 
                                                                  Nx)), mIP = as.integer(mIP), mdW = as.integer(mdW), 
                    mWS = as.integer(mWS), IP = as.integer(rep(0, mIP)), 
                    W = as.double(matrix(data = 0, nrow = mdW, ncol = Nx + 
                                           1)), WS = as.double(rep(0, mWS)), lpr = as.integer(lpr), 
                    ProgOpt = as.double(ProgOpt), verbose = as.logical(verbose), 
                    IsError = as.logical(IsError))
    if (any(is.infinite(sol$nX))) 
      sol$IsError <- TRUE
    if (fulloutput) {
      covar <- matrix(data = sol$W, nrow = mdW, ncol = Nx + 
                        1)[1:Nx, 1:Nx]
      RankEq <- sol$IP[1]
      RankApp <- sol$IP[2]
    }
  }
  else if (type == 2) {
    if (!is.null(Wx)) 
      stop("cannot solve least squares problem - weights not implemented for type 2")
    if (!is.null(Wa)) 
      stop("cannot solve least squares problem - weights not implemented for type 2")
    dvec <- crossprod(A, B)
    Dmat <- crossprod(A, A)
    diag(Dmat) <- diag(Dmat) + 1e-08
    Amat <- t(rbind(E, G))
    bvec <- c(F, H)
    
    
    sol <- tryCatch(solve.QP(Dmat, dvec, Amat, bvec, meq = Neq),
             warning = function(w){list(solution = NULL)},
             error = function(e){list(solution = NULL)})

    sol$IsError <- FALSE
    sol$X <- sol$solution
  }
  else stop("cannot solve least squares problem - type unknown")
  if(is.null(sol$solution)){
    res <- list(X = Inf, residualNorm = 0, solutionNorm = Inf, 
                IsError = TRUE, type = "lsei")
  }
  
  else{
    X <- sol$X
    X[which(abs(X) < Tol)] <- 0
    if (any(is.infinite(X))) {
      residual <- Inf
      solution <- Inf
    }
    else {
      residual <- 0
      if (Nin > 0) {
        ineq <- G %*% X - H
        residual <- residual - sum(ineq[ineq < 0])
      }
      if (Neq > 0) 
        residual <- residual + sum(abs(E %*% X - F))
      if (residual > Tol) 
        sol$IsError <- TRUE
      solution <- 0
      if (Napp > 0) 
        solution <- sum((A %*% X - B)^2)
    }
    xnames <- colnames(A)
    if (is.null(xnames)) 
      xnames <- colnames(E)
    if (is.null(xnames)) 
      xnames <- colnames(G)
    names(X) <- xnames
    res <- list(X = X, residualNorm = residual, solutionNorm = solution, 
                IsError = sol$IsError, type = "lsei")
    if (fulloutput && type == 1) {
      res$covar <- covar
      res$RankEq <- sol$IP[1]
      res$RankApp <- sol$IP[2]
    }
  }
  return(res)
}










#' calculated the membership of each node given ratio matrix and community centers
#' 
#' @param R n-by-(K-1) ratio matrix
#' @param vertices K-by-(K-1) community centers
#' @param K number of communities.
#' 
#' @return n-by-K membership matrix
#' 
#' @export
getMembership <- function(R, vertices, K, eig.values, eig.vectors){
  n = nrow(R)
  memberships = matrix(0, nrow = n, ncol = K)
  if (K == 2) {
    i1 <- which.min(vertices)
    i2 <- which.max(vertices)
    v1 <- min(vertices)
    v2 <- max(vertices)
    R_t <- (R - v1) / (v2 - v1)
    R_t[R_t < 0] <- 0
    R_t[R_t > 1] <- 1
    memberships[,i2] <- R_t
    memberships[,i1] <- 1 - R_t
    print('v2:')
    print(v2)
    print(eig.values[1:2])
    b1 <- sqrt(eig.values[1] + eig.values[2] * v1^2)
    b2 <- sqrt(eig.values[1] + eig.values[2] * v2^2)
    
    memberships[,i2] <- memberships[,i2] / b2
    memberships[,i1] <- memberships[,i1] / b1
    
    degrees <- abs(eig.vectors[,1] * (rowSums(memberships)))
  }
  else {
    for(i in 1:n){
      out = limSolve::lsei(A = t(vertices), 
                           B = R[i,], 
                           E = rep(1,K), F = 1,
                           G = diag(K), H = rep(0,K),
                           type = 2)
      memberships[i, ] = out$X
    }
    
    # truncate and normalize
    memberships[memberships > 1] = 1
    memberships[memberships < 0] = 0
    for(i in 1:n){
      memberships[i,] = memberships[i,]/sum(memberships[i,])
    }
    

    
    # recover the original memberships
    tildeV = cbind(rep(1,K), vertices)
    
    
    b1.inv = sqrt(diag(tildeV%*%diag(eig.values[1:K])%*%t(tildeV))) 
    b1 <- 1 / b1.inv
    
    # get degrees
    degrees = abs(eig.vectors[,1] * rowSums(memberships %*% diag(b1.inv)) )
    
    # get tilted memberships
    memberships = memberships%*%diag( b1.inv) 
    
    
    # normalize again
    for(i in 1:n){
      memberships[i,] = memberships[i,]/sum(memberships[i,])
    }
  }
  
  return(list(memberships = memberships, degrees = degrees, b1 = b1))
}





get_pi <- function(R,vertices,b1.inv){
  # R : n-by-(K-1) 
  # vertices: K-by-(K-1)
  n <- dim(R)[1]
  K <- dim(R)[2] + 1
  memberships = matrix(0, nrow = n, ncol = K)
  for(i in 1:n){
    out = limSolve::lsei(A = t(vertices), 
                         B = R[i,], 
                         E = rep(1,K), F = 1,
                         G = diag(K), H = rep(0,K),
                         type = 2)
    memberships[i, ] = out$X
  }
  
  # truncate and normalize
  memberships[memberships > 1] = 1
  memberships[memberships < 0] = 0
  for(i in 1:n){
    memberships[i,] = memberships[i,]/sum(memberships[i,])
  }
  
  
  
  # recover the original memberships
 
  # get tilted memberships
  memberships = memberships%*%diag( b1.inv) 
  
  
  # normalize again
  for(i in 1:n){
    memberships[i,] = memberships[i,]/sum(memberships[i,])
  }
  return (memberships)
}




####### get the author-level topic weights #######
get_author_tw <- function(){
  av_tw_all <- rowSums(W_hat.normed) / dim(W_hat.normed)[2]
  av_tw_all <- av_tw_all / sum(av_tw_all)
  
  topics <- c("Exp.Design","Inference","Mach.Learn.","Var.Select.","Time Series","Clinic.","Bayes","Math.Stats.",
              "Regression","Hypo.Test","Bio./Med.")
  
  
  topic_id <- 11
  topic_weight <- rep(0,dim(A)[1])
  
  id_select <- id_all
  author_select <- name
  
  author_interest <- id_select
  
  n_interest <- length(author_interest)
  paper_interest <- c()
  for (i in author_interest){
    paper_interest <- append(paper_interest,author_pap_list[[i]])
  }
  
  paper_interest <- unique(paper_interest)
  
  
  time_window <- c(1975:2015)
  paper_interest <- paper[paper$mr %in% paper_interest,]
  paper_interest <- paper_interest[paper_interest$year %in% time_window,]
  paper_author_int <- paper_author[names(paper_author) %in% paper_interest$mr]
  
  #### constuct a dataframe ####
  index_dc <- paper_time$mr[valid_ix] %in% paper_interest$mr
  
  paper_dc <- paper_time[valid_ix,][index_dc,]
  W_dc <- W_hat.normed[,index_dc]
  
  tmp <- data.frame(t(W_dc))
  colnames(tmp) <- topics
  paper_dc <- paper_dc <- paper_time[valid_ix,][index_dc,]
  cd_text_df <- data.frame(paper_dc[,c(1,5)])
  cd_text_df <- cbind(cd_text_df,tmp)
  
  
  paper_edge <- c()
  author_edge <- c()
  for(i in 1:length(paper_interest$mr)){
    aut_list <- paper_author_int[[i]]$id
    author_edge <- append(author_edge,aut_list)
    paper_edge <- append(paper_edge,rep(paper_interest$mr[i],length(aut_list)))
  }
  
  paper_author_edge <- data.frame(mr = paper_edge,author = author_edge)
  
  
  cd_text_df <- merge(x = cd_text_df,y = paper_author_edge,
                      by = "mr",all.y = TRUE)
  cd_text_df[is.na(cd_text_df)] <- 0
  
  
  ##### topic weights of all authors #####
  author_av_tw <- aggregate(.~author,cd_text_df[,colnames(cd_text_df) %in% c(topics,"author")],mean)
  a <- (author_av_tw[,colnames(author_av_tw) %in% c(topics)])
  a.normed <- a
  for (i in 1:dim(a)[1]){
    if (sum(a[i,]) != 0){
      tmp <- a[i,]
      a.normed[i,] <- tmp / sum(tmp)
    }
    else{
      a.normed[i,] <- 0
    }
  }
  author_av_tw[,colnames(author_av_tw) %in% c(topics)] <- a.normed
  
  return(author_av_tw)
}