#' kpodclustr demonstration
#' 
#' This file includes:
#'   - comments for all functions in kpodclustr;
#'   - example output and explanation for each function;
#'   - code to assist making of the vignette for the kpodclustr package.
#'   
#' This file serves the purpose of:
#'   - increase understanding of the purpose of kpodcluster;
#'   - familiarize with the inner-workings of each function/each step of kpod;
#'   - assist making of the vignette for the kpodclustr package.


#' Generate Test Data: complete and noisy
#' @param p number of features for each observation
#' @param n number of observation
#' @param k number of clusters/centroids
#' @param sigma variance
#' @param missing percentage of missing data out of all p*k
#' @param seed memory lot
#' B: centroid features matrix k*p
#' A: centroid assignment
#' 

genTestData <- function(p,n,k,sigma,missingpct,seed=1991){ # sigma is standard deviation, maybe add a mean of generated data?
  
  # Generate complete data
  set.seed(seed) # set memory location
  B <- matrix(rnorm(k*p,0,1),k,p) # Gen centroid features matrix
  A <- sample(1:k,n,replace=TRUE) # Gen centroid assignment for obs', w/ replacement
  X <- B[A,] + sigma*matrix(rnorm(n*p,0,1),n,p) # each observation: assigned centroid + deviation factor(sigma) * random value from standard normal distribution
  
  # Generate noisy data
  
  X_missing <- X # copy the complete data
  missing_index <- sample(1:(n*p),(n*p*missingpct),replace=TRUE) # Gen index for missing observations
  X_missing[missing_index] <- NA # assign NA to data points with missing index
  
  return(list(No_Missing = X, With_Missing = X_missing, Centroids = B, Assignment = A))  
  
  #' Example output when there are k = 3 centroids assigned to n = 10 observations, each with p = 2 features, with standard deviation sigma = 0.25, and 4 missing entries, missingpct = 0.2 -- seed = 1991:
  #' $No_Missing
  # [,1]         [,2]
  # [1,] -1.305041  1.231794130
  # [2,]  1.175285  0.592935726
  # [3,] -0.900364  0.676791746
  # [4,] -2.286674 -0.529542588
  # [5,]  1.046599  0.088004503
  # [6,] -1.982720 -0.006902344
  # [7,] -2.537605 -0.408056101
  # [8,] -2.408308 -0.166181343
  # [9,] -2.223123 -0.610156411
  # [10,] -1.220365  0.911779269
  # 
  # $With_Missing
  # [,1]         [,2]
  # [1,] -1.305041           NA
  # [2,]  1.175285  0.592935726
  # [3,] -0.900364  0.676791746
  # [4,] -2.286674 -0.529542588
  # [5,]  1.046599  0.088004503
  # [6,]        NA -0.006902344
  # [7,]        NA -0.408056101
  # [8,] -2.408308           NA
  # [9,] -2.223123 -0.610156411
  # [10,] -1.220365  0.911779269
  # 
  # $Centroids
  # [,1]       [,2]
  # [1,] -1.0337647  0.9483830
  # [2,] -2.2120213 -0.3159648
  # [3,]  0.8166651  0.2340349
  # 
  # $Assignment
  # [1] 1 3 1 2 3 2 2 2 2 1
}


findMissing <- function(X){
  # find index of all missing data in data matrix X, and put them in a list
  missing_all <- which(is.na(X))
  return(missing_all)
  
  # Example output from same parameters as generation step (seed = 1991)
  # findMissing(X_missing)
  # [1]  6  7 11 18
  # > X_missing
  # [,1]         [,2]
  # [1,] -1.305041           NA
  # [2,]  1.175285  0.592935726
  # [3,] -0.900364  0.676791746
  # [4,] -2.286674 -0.529542588
  # [5,]  1.046599  0.088004503
  # [6,]        NA -0.006902344
  # [7,]        NA -0.408056101
  # [8,] -2.408308           NA
  # [9,] -2.223123 -0.610156411
  # [10,] -1.220365  0.911779269
}

initialImpute <- function(X){
  # Calculate mean value for each feature based on all observation, ignoring the missing data.
  avg <- mean(X,na.rm=TRUE)  
  
  # fill missing with mean value for that feature.
  X[which(is.na(X))] <- avg
  # return filled data matrix
  return(X)
  
  #'Example output when seed = 1991:
  # initialImpute(X_missing)
  # [,1]         [,2]
  # [1,] -1.3050410 -0.462946062
  # [2,]  1.1752849  0.592935726
  # [3,] -0.9003640  0.676791746
  # [4,] -2.2866741 -0.529542588
  # [5,]  1.0465990  0.088004503
  # [6,] -0.4629461 -0.006902344
  # [7,] -0.4629461 -0.408056101
  # [8,] -2.4083082 -0.462946062
  # [9,] -2.2231227 -0.610156411
  # [10,] -1.2203647  0.911779269
  # 
  # Reminder: X_missing
  # [,1]            [,2]
  # [1,] -1.305041           NA
  # [2,]  1.175285  0.592935726
  # [3,] -0.900364  0.676791746
  # [4,] -2.286674 -0.529542588
  # [5,]  1.046599  0.088004503
  # [6,]        NA -0.006902344
  # [7,]        NA -0.408056101
  # [8,] -2.408308           NA
  # [9,] -2.223123 -0.610156411
  # [10,] -1.220365  0.911779269
}

#' Choose initial centers with complete dataset without missing
kmpp <- function(X, k) {
  # retrieve basic statistics about data matrix
  n <- nrow(X)
  p <- ncol(X)
  
  ## Initial Centroids are selected from observations
  # C stores Indices of observations that are chosen to be centroids
  C <- integer(k) # there are k centroids, so, k indices. 
  C[1] <- sample(1:n, 1) # first one is random out of n obs
  
  for (i in 2:k) {
    S <- matrix(NA,n,i-1)
    for (j in 1:(i-1)) {
      S[,j] <- apply(X -
                       matrix(X[C[j],],n,p,byrow=TRUE),1,FUN=function(x)
                       {norm(as.matrix(x),'f')**2})
    }
    D <- apply(S,1,min)
    pr <- D/sum(D)
    C[i] <- sample(1:n, 1, prob = pr)
  }
  return(X[C,])
  #' The result of this kmeans ++ procedure can be different everytime it runs, but none the less it should look similar to this -- 
  #' Example output for seed = 1991:
  #' kmpp(X1991,k)
  # [,1]       [,2]
  # [1,]  1.175285  0.5929357
  # [2,] -2.537605 -0.4080561
  # [3,] -0.900364  0.6767917
  #
  # Reminder: X1991
  # [,1]         [,2]
  # [1,] -1.305041  1.231794130
  # [2,]  1.175285  0.592935726
  # [3,] -0.900364  0.676791746
  # [4,] -2.286674 -0.529542588
  # [5,]  1.046599  0.088004503
  # [6,] -1.982720 -0.006902344
  # [7,] -2.537605 -0.408056101
  # [8,] -2.408308 -0.166181343
  # [9,] -2.223123 -0.610156411
  # [10,] -1.220365  0.911779269
 
}

#'Side note for assigning clusters ++
#'In k-means:
#'Total sum of square = distances between all obs and global center
#'Within-cluster sum of square = sum of distances between obs of one cluster to that cluster center
#'Between-cluster sum of square = distances between cluster centers
#'totss = withinss + betweenss
#'
#'fit 1-(sum(res$withinss)/res$totss) ; The greater the value(closer to 1), the better the fit

assign_clustpp <- function(X,init_centers,kmpp_flag=TRUE,max_iter=20){
  # Run kmeans on filled matrix with given initial centers
  res <- kmeans(X, init_centers)
  
  clusts <- res$cluster # assignment vector, each element is a integer 1:k
  # For seed 1991: res1991$cluster
  #[1] 3 2 3 1 2 1 1 1 1 3
  
  obj <- res$totss # total sum of square
  # For seed 1991: res1991$totss
  # [1] 20.53454
  
  fit <- 1-(sum(res$withinss)/res$totss) # The greater the value(closer to 1), the better the fit
  # for seed 1991: 1-(sum(res1991$withinss)/res1991$totss)
  # [1] 0.9605108
  
  centers <- res$centers # centroids matrix
  # For seed 1991: res1991$centers
  # [,1]       [,2]
  # 1 -2.287686 -0.3441678
  # 2  1.110942  0.3404701
  # 3 -1.141923  0.9401217
  
  
  if (kmpp_flag == TRUE) {
    ## Try to find one better assignment solution, then jump out of interations upon finding ONE.
    for (iter in 1:max_iter) {
      centers_kmpp <- kmpp(X,length(res$size))
      sol <- kmeans(X, centers_kmpp)
      if (sol$totss < obj) {
        obj <- sol$totss
        clusts <- sol$cluster
        fit <- 1-(sum(sol$withinss)/sol$toss)
        centers <- sol$centers
        break
      }
    }
  }
  return(list(clusts=clusts,obj=obj,centers=centers,fit=fit))
}

# The main course served last.
kpod <- function(X,k,kmpp_flag=TRUE,maxiter=10) {
  # retrieve basic statistics about data matrix
  n <- nrow(X) 
  p <- ncol(X)
  
  # memorize the locations of missing indices so we can update them every iteration
  missing <- findMissing(X) 
  
  # initialize tracking variables
  cluster_vals <- vector(mode="list",length=maxiter)
  obj_vals <- double(maxiter)
  fit <- double(maxiter)
  
  # Run the initial imputation to fill missings with mean
  X_copy <- initialImpute(X)
  
  print(list(X_filled = X_copy))
  
  ## Following is basically one iteration of assign_clustpp()
  # Run KMPP first time to select initial centers
  init_centers <- kmpp(X_copy, k)
  
  print(list(initialcenter = init_centers))
  
  # Run first kmeans to initial assignment
  temp <- kmeans(X_copy,init_centers)
  clusts <- temp$cluster
  centers <- temp$centers
  fit[1] <- 1-(sum(temp$withinss)/temp$totss)
  
  print(list(clusters=clusts, centers=centers, fit = fit[1]))
  
  clustMat <- centers[clusts,]  # cluster matrix
  print(list(clustermatrix = clustMat))
  
  # Update originally missing data with first interation k-means centers
  X_copy[missing] <- clustMat[missing]
  
  ## Done with the first interation
  
  #' What these variables may look like
  # $clusters
  # [1] 3 2 3 1 2 3 3 1 1 3
  # 
  # $centers
  # [,1]       [,2]
  # 1 -2.3060350 -0.5342150
  # 2  1.1109420  0.3404701
  # 3 -0.8703324  0.1421333
  # 
  # $fit
  # [1] 0.8647777 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
  # [8] 0.0000000 0.0000000 0.0000000
  # 
  # $clustermatrix
  # [,1]       [,2]
  # 3 -0.8703324  0.1421333
  # 2  1.1109420  0.3404701
  # 3 -0.8703324  0.1421333
  # 1 -2.3060350 -0.5342150
  # 2  1.1109420  0.3404701
  # 3 -0.8703324  0.1421333
  # 3 -0.8703324  0.1421333
  # 1 -2.3060350 -0.5342150
  # 1 -2.3060350 -0.5342150
  # 3 -0.8703324  0.1421333
  
  obj_vals[1] <- sum((X[-missing]-clustMat[-missing])^2)
  cluster_vals[[1]] <- clusts
  
  # Run remaining iterations
  for (i in 2:maxiter){
    temp <- assign_clustpp(X_copy,centers,kmpp_flag)
    clusts <- temp$clusts
    centers <- temp$centers
    fit[i] <- temp$fit
    
    # Impute clusters
    clustMat <- centers[clusts,]
    X_copy[missing] <- clustMat[missing]
    
    obj_vals[i] <- sum((X[-missing]-clustMat[-missing])**2)
    cluster_vals[[i]] <- clusts
    
    if (all(cluster_vals[[i]] == cluster_vals[[i-1]])){
      noquote('Clusters have converged.')
      return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i],centers=clustMat))
      break
    }
  }
  return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i],centers=clustMat))

}

# EX1 small data set
## Retrieve Data 1991
data1991 <- genTestData (p=2,n=10,k=3,sigma=0.25,missingpct=0.90,seed=1991)
X1991 <- data1991[[1]]
Xm1991 <- data1991[[2]]
B1991 <- data1991[[3]]
A1991 <- data1991[[4]]
## run k-mean on 100% data 1991.
init_centers1991 <- kpodclustr::kmpp(X1991,3)
km1991 <- kmeans(X1991,init_centers1991)
fit_km1991 <- 1-(sum(km1991$withinss)/km1991$totss)
## run k-pod on 80% data1991.
kp1991 <- kpodclustr::kpod(Xm1991,3)
fit_kpod1991 <- kp1991$fit
## compare clustering using rand.index
rand1991 <- adj.rand.index(km1991$cluster,kp1991$cluster)

# Ex2: Larger data set 
## Retrieve data 1992 (n = 100, k = 5, p = 2)
data1992 <- genTestData (p=2,n=100,k=5,sigma=0.25,missingpct=0.10,seed=1992)
X1992 <- data1992[[1]]
Xm1992 <- data1992[[2]]
B1992 <- data1992[[3]]
A1992 <- data1992[[4]]
## run k-mean on 100% data 1992.
init_centers1992 <- kpodclustr::kmpp(X1992,5)
km1992 <- kmeans(X1992,init_centers1992)
fit_km1992 <- 1-(sum(km1992$withinss)/km1992$totss)
## run k-pod on 90% data1992.
kp1992 <- kpodclustr::kpod(Xm1992,5)
fit_kpod1992 <- kp1992$fit

# Scatter plot
## EX1 small data set
### original data scatter
ggplot(data = as.data.frame(X1991), mapping = aes(x=X1991[,1],y=X1991[,2]))+ geom_point()
### colored by k means
ggplot(data = as.data.frame(X1991), mapping = aes(x=X1991[,1],y=X1991[,2], color=factor(km1991$cluster))) + geom_point() + theme(legend.position = "none") + geom_point(data = as.data.frame(km1991$centers), mapping = aes(x=km1991$centers[,1],y=km1991$centers[,2], color = "red",size = 3))
### colored by k pod
ggplot(data = as.data.frame(X1991), mapping = aes(x=X1991[,1],y=X1991[,2], color=as.character(kp1991$cluster))) + geom_point() + theme(legend.position = "none") 


# make code k-mean on complete k-pod on missing, and compare result
# missing fraction 10% - 90% performance comparison between k-mean k-pod
# RandIndex 
# ARI adjusted RandIndex assess the similarity between clustrings
# illustrate how to measure similarity between clustrings
# illustrate it gets harder as missing fraction increase