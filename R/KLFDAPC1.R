### Function is built based on DA

requireNamespace("SNPRelate")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!requireNamespace("SNPRelate", quietly=TRUE))
  BiocManager::install("SNPRelate")

KLFDAPC=function(infile,y, n.pc,sample.id=NULL, snp.id=NULL,
                 autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                 algorithm=c("exact", "randomized"),
                 eigen.cnt=ifelse(identical(algorithm, "randomized"), 16L, 32L),
                 num.thread=1L, bayesian=FALSE, need.genmat=FALSE,
                 genmat.only=FALSE, eigen.method=c("DSPEVX", "DSPEV"),
                 aux.dim=eigen.cnt*2L, iter.num=10L, verbose=TRUE,
                 kernel=kernlab::polydot(degree = 1, scale = 1, offset = 1), r=3, tol=1e-30,prior=NULL, CV=FALSE,usekernel = TRUE, fL = 0.5,metric = c('weighted', 'orthonormalized', 'plain'),
                 knn = 6, reg = 0.001,...){


require(SNPRelate)
genofile=snpgdsOpen(infile)
pcadata <- snpgdsPCA(genofile, sample.id=sample.id, snp.id=snp.id,
                     autosome.only=autosome.only, remove.monosnp=remove.monosnp, maf=maf, missing.rate=missing.rate,
                     algorithm=algorithm,
                     eigen.cnt=eigen.cnt,
                     num.thread=num.thread, bayesian=bayesian, need.genmat=need.genmat,
                     genmat.only=genmat.only, eigen.method=eigen.method,
                     aux.dim=aux.dim, iter.num=iter.num, verbose=verbose,...)
snpgdsClose(genofile)

### This function provides the unified solution for mutiple Kernel choice. the kernel choice here will allow users to choose any types of kernel they want, not only constrains on gaussian kernel

showfile.gds(closeall=TRUE)

require(lfda)
klfda_1=function (x, y,kernel=kernlab::polydot(degree = 1, scale = 1, offset = 1), r=3, tol=1e-10,prior=NULL, CV=FALSE,usekernel = TRUE, fL = 0.5,metric = c('weighted', 'orthonormalized', 'plain'),
                  knn = 6, reg = 0.001,...) {## reg regularization parameter (default: 0.001)

  #

  tol=tol
  obj.trainData = x
  obj.trainClass = as.factor(y)
  obj.classes = sort(unique(obj.trainClass));
  obj.nClasses = length(obj.classes)
  requireNamespace("kernlab")
  k=kernlab::kernelMatrix(kernel, obj.trainData,obj.trainData)
  #k=multinomial_kernel(obj.trainData,obj.trainData,order=2)

  obj.nObservations=dim(obj.trainData)[1]
  obj.nFeatures = dim(obj.trainData)[2]
  #k=k/obj.nObservations

  cl <- match.call()
  metric <- match.arg(metric) # the type of the transforming matrix (metric)
  x=as.matrix(k)
  p <- ncol(k)
  n <- nrow(k) # number of samples


  if(n != length(y))
    stop("nrow(x) and length(y) are different")
  g <- as.factor(y)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))

  if(is.null(prior)==FALSE) {
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]

  }


  if(any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(ngettext(length(empty),
                             "group %s is empty",
                             "groups %s are empty"),
                    paste(empty, collapse = " ")), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }

  proportions <- counts/n
  ng <- length(proportions)


  group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)

  ### new section  from local da

  y <- t(as.matrix(y)) # transpose of original class labels
  d=nrow(k)
  if(is.null(r)) r <- n # if no dimension reduction requested, set r to n

  tSb <- mat.or.vec(n, n) # initialize between-class scatter matrix (to be maximized)
  tSw <- mat.or.vec(n, n) # initialize within-class scatter matrix (to be minimized)
  require(lfda)
  # compute the optimal scatter matrices in a classwise manner
  for (i in unique(as.vector(t(y)))) {

    Kcc <- k[y == i, y == i] # data for this class
    Kc <- k[, y == i]
    nc <- nrow(Kcc)

    # Define classwise affinity matrix
    Kccdiag <- diag(Kcc) # diagonals of the class-specific data
    distance2 <- repmat(Kccdiag, 1, nc) + repmat(t(Kccdiag), nc, 1) - 2 * Kcc

    # Get affinity matrix
    A <- lfda::getAffinityMatrix(distance2, knn, nc)

    Kc1 <- as.matrix(rowSums(Kc))
    Z <- Kc %*% (repmat(as.matrix(colSums(A)), 1, n) * t(Kc)) - Kc %*% A %*% t(Kc)
    tSb <- tSb + (Z/n) + Kc %*% t(Kc) * (1 - nc/n) + Kc1 %*% (t(Kc1)/n)
    tSw <- tSw + Z/nc
  }

  K1 <- as.matrix(rowSums(k))
  tSb <- tSb - K1 %*% t(K1)/n - tSw

  tSb <- (tSb + t(tSb))/2 # final between-class cluster matrix
  tSw <- (tSw + t(tSw))/2 # final within-class cluster matrix

  # find generalized eigenvalues and normalized eigenvectors of the problem
  eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw + reg * diag(1, nrow(tSw), ncol(tSw)),tol=tol,bounds = list(a=c(x> 0))) %*% tSb,
                                           k = r,which ='LM')) # r largest magnitude eigenvalues
  eigVec <- Re(eigTmp$vectors) # the raw transforming matrix
  eigVal <- as.matrix(Re(eigTmp$values))

  # options to require a particular type of returned transform matrix
  # transforming matrix (do not change the "=" in the switch statement)
  Tr <- getMetricOfType(metric, eigVec, eigVal, n)


  Z <- t(t(Tr) %*% k) # transformed data




  classVecsTrain = matrix(nrow=obj.nObservations, ncol=obj.nClasses);
  obj.nObsPerClas = matrix(nrow=1, ncol=obj.nClasses);
  for (i in 1:obj.nClasses) {
    clas = obj.classes[i]
    classVecsTrain[, i] = match(obj.trainClass, clas,nomatch = 0)

    obj.nObsPerClas[i] = sum(classVecsTrain[,i])
  }



  redZ=list()

  obj.means1 = list()
  obj.covariances1=list()
  for (i in 1 : obj.nClasses){

    redZ[[i]]=list()
    #redProjData[[i]] = projData[as.logical(classVecsTrain[,i]), (1 : (obj.nClasses - 1))]### here obj,nClasses -1 can be subsetuted by d, the reduced number of dim
    redZ[[i]] = Z[as.logical(classVecsTrain[,i]), (1 : r)]

    obj.means1[[i]]= colMeans(redZ[[i]])
    obj.covariances1[[i]]=list()
    obj.covariances1[[i]] = cov(redZ[[i]])
  }
  ## Generate priors and likelihood functions
  # Priors

  if (is.null(prior)==TRUE){
    obj.priors = matrix(data=0,nrow=1, ncol=obj.nClasses);
    for (i in  1 : obj.nClasses){
      obj.priors[,i] = obj.nObsPerClas[,i] / obj.nObservations;
    }
  }
  if(is.null(prior)==FALSE) {
    priorsTmp = prior;
    obj.priors = matrix(data=0,nrow=1, ncol=obj.nClasses)
    for (i in 1 : obj.nClasses){
      clas = obj.classes[i]
      obj.priors[,i]= priorsTmp[i]
    }
  }
  prior=obj.priors
  names(prior) <- names(counts) <- lev1
  # Likelihood
  obj.likelihoodFuns = list()
  # warning('error', 'nearlySingularMatrix')


  ### the methods of Z
  obj.covInv1=list()
  obj.covDet1=list()
  factor1=list()
  #obj.covInv=solve(obj.covariances[[i]])
  for (i in 1 : obj.nClasses){

    obj.covInv1[[i]]=list()
    obj.covDet1[[i]]=list()
    obj.covInv1[[i]] = solve(obj.covariances1[[i]],tol=tol,bounds = list(a=c(obj.covariances1[[i]]> 0)));# ginv(obj.co)
    obj.covDet1[[i]] = det(obj.covariances1[[i]]);

    #pvaCov = vpa(obj.covariances[[i]]);
    #obj.covInv[[i]] = double((pvaCov));
    #obj.covDet.(clas) = double(det(pvaCov));
    factor1[[i]]=list()
    factor1[[i]] = (2 * pi) ^ (-(obj.nClasses - 1) / 2) * (obj.covDet1[[i]] ^ -0.5)
  }


  likelihoods1 = matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses);

  posteriors1 = matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses)

  dist1=matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses)
  ## here from the  scripts  of Perr, each row of new data has nclass of
  # (rpdTest[j,] - obj.means[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[[i]])
  for (j in 1 : obj.nObservations){

    for (i in 1 : obj.nClasses){
      # clas = obj.classes[i]


      #mean((rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))

      dist1[j,i]=mean((Z[j,] - obj.means1[[i]]) %*% obj.covInv1[[i]] * t(Z[j,] -obj.means1[[i]]))
      # likelihoods[[j]][[i]] =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means[i]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[i]))
      #likelihoods1 =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))

      likelihoods1[j,i] =factor1[[i]] * exp(-0.5 * dist1[j,i])


    }

    posteriors1[j,] = likelihoods1[j,] * obj.priors / sum(likelihoods1[j,] * obj.priors)
  }
  #}
  ### posteriors.class is the normal prodiction using projData posteriors.classZ using Z

  ## Predicting the class of each data point


  posteriors.classZ=factor(obj.classes[max.col(posteriors1)], levels = obj.classes)




  # if(CV) {
  #   x <- t(t(Tr) %*% k) # transformed data Z
  #   dm <- group.means %*% Tr
  ### K is deternmined from n.pc, ng is number of gropus
  #   K = ng
  #   dist <- matrix(0, n, ng) # n row and ng col

  ## dev matrix
  #   for(i in 1L:ng) {
  #     dev <- x - matrix(dm[i,  ], n, r, byrow = TRUE)
  #      dist[, i] <- rowSums(dev^2)### dis or devation of each class
  #    }
  #   ind <- cbind(1L:n, g) ### give the order number (sequence) and group(class )
  #   nc <- counts[g] ## number of class
  #   cc <- nc/((nc-1)*(n-K)) ### proportation of eac
  #    dist2 <- dist
  #   for(i in 1L:ng) {
  #     dev <- x - matrix(dm[i,  ], n, r, byrow = TRUE)
  #     dev2 <- x - dm[g, ]
  #    tmp <- rowSums(dev*dev2)
  #     dist[, i] <- (n-1L-K)/(n-K) * (dist2[, i] +  cc*tmp^2/(1 - cc*dist2[ind]))
  #  }
  ### dist should be discriminat function
  #   dist[ind] <- dist2[ind] * (n-1L-K)/(n-K) * (nc/(nc-1))^2 /(1 - cc*dist2[ind])
  #    dist <- 0.5 * dist - matrix(log(prior), n, ng, byrow = TRUE) # proboloty of discrim density function

  #   dist <- exp(-(dist - min(dist, na.rm = TRUE))) #### density function proporbility of pi this is distribution Normal distribution function
  #   cl <- factor(lev1[max.col(dist)], levels = lev)
  ## calculate distance is
  #    deltaTrain[,i]=deltaTrain[,i]-0.5*t(Mu_k[,i])%*%sigmaMinus1%*%Mu_k[,i]+log(Pi_k[i])
  ##  convert to posterior probabilities

  #  posterior <- dist/drop(dist %*% rep(1, length(prior)))
  #    dimnames(posterior) <- list(rownames(x), lev1)
  #   return(list( eigTmp=eigTmp,eifVec=eigVec,eigVal=eigVal,bayes=bayes,bayes_judgement=bayes_judgement,bayes_assigment=bayes_assigment, post_class=cl,posterior = posterior,Z=x,T=Tr))
  #  }
  # xbar <- colSums(prior %*% group.means)
  #  fac <-1/(ng - 1)
  #  X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% Tr
  # X.s <- svd(X, nu = 0L)
  # rank <- sum(X.s$d > tol * X.s$d[1L])
  # rank=ng
  #  T_lda <- Tr %*% X.s$v[, 1L:rank]
  #  if(is.null(dimnames(x)))
  #   dimnames(T_lda) <- list(NULL, paste("LD", 1L:rank, sep = ""))
  #  else {
  #   dimnames(T_lda) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
  #    dimnames(group.means)[[2L]] <- colnames(x)
  # }
  #  svd = X.s$d[1L:rank]
  #  PLD1=svd[1]^2/sum(svd^2)
  # PLD2=svd[2]^2/sum(svd^2)
  # PLD=as.data.frame(cbind(PLD1,PLD2))
  require(WMDB) ## FOR dbayes function the  function
  require(DA)
  if (!requireNamespace("DA", quietly=TRUE))
  devtools::install_github("xinghuq/DA")
  bayes_judgement=DA::Mabayes(Z,obj.trainClass,var.equal = FALSE,tol=tol)
  require(klaR)
  bayes=klaR::NaiveBayes(as.data.frame(Z), obj.trainClass, prior=obj.priors, usekernel, fL = fL,...)
  bayes_assigment=predict(bayes,threshold = 0.001,dkernel=kernel)
  # t_lda=lda(Z,obj.trainClass,...)
  cl <- match.call()
  cl[[1L]] <- as.name("klfda")
  structure(list(kernel=kernel, tol=tol,usekernel=usekernel,fL=fL,obj.classes=obj.classes,obj.nObservations=obj.nObservations,obj.means1=obj.means1,obj.covInv1=obj.covInv1,factor1=factor1, eigTmp=eigTmp,eigVec=eigVec,eigVal=eigVal, prior = prior, bayes=bayes,posteriors.classZ=posteriors.classZ,posteriors1=posteriors1, bayes_judgement=bayes_judgement,bayes_assigment=bayes_assigment,means=group.means,
                 obj.priors=obj.priors,T=Tr,obj.nClasses=obj.nClasses,obj.trainData=obj.trainData, obj.trainClass=obj.trainClass, Z=Z,lev = lev, N = n, call = cl),class = "klfda")
}



#PCs=data.frame(sample.id = pcadata$sample.id,pcadata$eigenvect[,1:n.pc],stringsAsFactors = FALSE)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#PCs_norm=apply(pcadata$eigenvect[,1:n.pc],2,normalize)

#### Note, it's better not to normzlize the PCs before do kernel computations or you will need to choose different parameters (e.g., sigma) in order to adapt to the transformed matrix.
## It is easy to find a sigma if you do not normalize the PCs. However, if you have normalized the PCs, you probably get a warning:
#A[flag] <- exp(-distance2[flag]/localscale[flag]) :
#NAs are not allowed in subscripted assignments In addition: Warning message:
# In sqrt(kNNdist2) : NaNs produced
print("Doing KLFDAPC")

klfdapcs=klfda_1(pcadata$eigenvect[,1:n.pc],y,kernel,r,tol,...)

out=list(KLFDAPC=klfdapcs,PCA=pcadata)
 class(out) <- "KLFDAPC"
return(out)


}









##### predict, input PCs


predict.klfdapc=function(object,newdata,dimen,...){
  tol=object$KLFDAPC$tol
  nObsTest=dim(newdata)[1]
  nFeaturesTest = dim(newdata)[2]
  obj.nClasses=object$KLFDAPC$obj.nClasses
  obj.nFeatures=object$KLFDAPC$obj.nFeatures
  obj.nObservations=object$KLFDAPC$obj.nObservations
  require(kernlab)
  kernel=object$KLFDAPC$kernel
  knewdata=kernelMatrix(kernel, object$KLFDAPC$obj.trainData,newdata)
  #knewdata=multinomial_kernel(object$obj.trainData,newdata,order=2)
  #newdata=as.matrix(newdata)
  knewdata=knewdata/obj.nObservations
  Z=object$KLFDAPC$Z
  Y=object$KLFDAPC$obj.trainClass
  Trans=as.matrix(object$KLFDAPC$T)
  #Z2=as.matrix(knewdata) %*% Trans
  Z2=t(as.matrix(knewdata)) %*% Trans
  # if (is.null(prior)==TRUE){
  prior=object$KLFDAPC$obj.priors
  # }
  usekernel=object$KLFDAPC$usekernel
  fL=object$KLFDAPC$fL
  # t_lda=object$t_lda
  require(klaR) ## FOR Nativebaye function
  bayes_jud_pred=klaR::Mabayes(Z,Y,TstX = Z2,var.equal = FALSE,tol=tol)
  # bayes=NaiveBayes(Z, Y, prior, usekernel, fL,kernel,bw = "nrd0", adjust = 1,weights = NULL, window = kernel, give.Rkern = FALSE,...)
  Nbayes_assig_pred=predict(object$KLFDAPC$bayes,as.data.frame(Z2),threshold = 0.001,dkernel=kernel,...)
  #  p_lda=predict(t_lda,Z2,...)

  ## Retrieving the likelihood of each test point

  likelihoods1 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);

  posteriors1 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)

  dist1=matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)
  ## here from the  scripts  of Perr, each row of new data has nclass of
  # (rpdTest[j,] - obj.means[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[[i]])

  obj.means1=object$KLFDAPC$obj.means1
  obj.covInv1=object$KLFDAPC$obj.covInv1

  factor1=object$KLFDAPC$factor1
  obj.classes=object$KLFDAPC$obj.classes
  obj.priors=object$KLFDAPC$obj.priors

  for (j in 1 : nObsTest){

    for (i in 1 : obj.nClasses){
      # clas = obj.classes[i]


      dist1[j,i]=mean((Z2[j,] - obj.means1[[i]]) %*% obj.covInv1[[i]] * t(Z2[j,] -obj.means1[[i]]))
      # likelihoods[[j]][[i]] =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means[i]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[i]))
      #likelihoods1 =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))

      likelihoods1[j,i] =factor1[[i]] * exp(-0.5 * dist1[j,i])


    }

    posteriors1[j,] = likelihoods1[j,] * obj.priors / sum(likelihoods1[j,] * obj.priors)## Z
  }
  #}


  ## Predicting the class of each data point


  posteriors.class1=factor(obj.classes[max.col(posteriors1)], levels = obj.classes)

  lev=object$KLFDAPC$lev
  ng=length(lev)

  means <- colSums(prior %*% object$KLFDAPC$means)
  scaling <- object$KLFDAPC$scaling
  x <-Z2
  dm <- scale(object$KLFDAPC$means, center = means, scale = FALSE) %*% Trans

  dimen <- if(missing(dimen)) length(object$KLFDAPC$svd) else min(dimen, length(object$KLFDAPC$svd))
  N <- object$KLFDAPC$N

  dm <- dm[, 1L:dimen, drop = FALSE]
  dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x),
                 length(prior), byrow = TRUE) - x[, 1L:dimen, drop=FALSE] %*% t(dm)
  dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))

  posterior <- dist / drop(dist %*% rep(1, length(prior)))
  nm <- names(object$KLFDAPC$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$KLFDAPC$lev)
  dimnames(posterior) <- list(rownames(x), nm)
  list(class = cl, posterior = posterior,posteriors1=posteriors1, posteriors.class1=posteriors.class1,x = x[, 1L:dimen, drop = FALSE],bayes_jud_pred=bayes_jud_pred,bayes_assig_pred=Nbayes_assig_pred)
}



