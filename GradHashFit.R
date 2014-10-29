## Last Update: 10.22.14

##### Miscellaneous Support Functions ##### 

# For binomial family calculations
expit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(p) {log(p/(1-p))}

# Soft thresholding
soft <- function(x, lambda){
  x.soft <- sign(x) * max(abs(x) - lambda, 0)
  return(x.soft)
}

##### Data Preparation ##### 

# take a vector of filenames or strings and number of available memory locations 
# and hash the text to those locations by word or tuple
hashText <- function(x, type=c("textfile","textstring"), num.mem,tuple){
  library(hashFunction)
  if (type == "textfile") {
    str <- paste(readLines(x))
  } else if (type == "textstring") {
    str <- x
  } else {
    print("Please input a file or string and specify type.") 
  }
  str <- gsub("[[:punct:]]","",str)
  str <- tolower(str)
  str <- strsplit(str, ' ')
  words = c()
  for (element in str[[1]]) {
    if (element != "") {
      words <- c(words, element)
    }
  }
  if (tuple == TRUE) {
    nwords <- length(words)
    tups <- c()
    for (i in c(1:(nwords-1))){
      this.tup <- paste(words[i],words[i+1])
      tups <- c(tups, this.tup)
    }
    words <- c(words, tups)  
  }
  words <- sapply(words, function(x){ murmur3.32(x) %% num.mem })
  return(words)
}

# Hash categories of high-dimensional categorical variable
hashCategories <- function(x, num.mem){
  library(hashFunction)
  cats <- sapply(x, function(y){ murmur3.32(y) %% num.mem })
  return(cats)
}

# take a data.frame and a vector of types, and export a list of lists 
## data.frame of explanatory variables (n rows, p columns) 
## types can be: "category", "smooth", "linear", "textfile", "textstring" 
## family = c("binomial","linear") 
## n.mem.vec is a vector of the number of memory locations (p-length; can have NAs in non-cat/text vars)

data.prep <- function(data, types, n.mem.vec, tuples=FALSE){
  p <- ncol(data) 
  out <- vector("list", p)
  for (i in c(1:p)){
    thistype <- types[i]
    if (thistype == "category"){
      out[[i]] <- lapply(data[,i], FUN = function(x){hashCategories(as.character(x), n.mem.vec[i])})
    } else if (thistype == "smooth" | thistype == "linear"){
      out[[i]] <- data[,i]
    } else if (thistype == "textfile" | thistype == "textstring"){
      out[[i]] <- lapply(data[,i], FUN = function(x){hashText(as.character(x), thistype, n.mem.vec[i], tuples)})
    } else { (out[[i]] <- "Not a recognized type")}
  }
  return(out)
}

##### Gradient Calculation ##### 

# In each case, x is a column vector 
# Resids calculation is the only difference between quantitative & binomial outcomes 
# gradient for smooth is just residuals 

calcGradient.category <- function(hashed.cats, mem.avail, resids){   
  gradient.cat <- rep(0, mem.avail)
  n <- length(resids) 
  for (i in 1:n) {    # For each person 
    this.hash <- hashed.cats[[i]]
    gradient.cat[this.hash] <- gradient.cat[this.hash] + resids[i]
  }
  return(gradient.cat)
}

# text variables (original input as either file or string)
calcGradient.text <- function(hashed.text, mem.avail, resids){
  gradient.text <- rep(0, mem.avail)
  n <- length(resids) 
  for (i in 1:n){
    for (j in 1:length(hashed.text[[i]])){
      this.hash <- hashed.text[[i]][j] 
      gradient.text[this.hash] <- gradient.text[this.hash] + resids[i]
    }
  }
  return(gradient.text)
}

# linear variables 
calcGradient.linear <- function(x, resids){
  gradient.linear <- t(x) %*% resids  
  return(gradient.linear) 
}

# Overall 
calcGradient <- function(datalist, mem.avail.vec, resids, types){
  gradient <- vector("list",length(types))
  for(i in 1:length(types)){
    thistype <- types[i]
    if (thistype=="category"){
      gradient[[i]] <- calcGradient.category(datalist[[i]], mem.avail.vec[i], resids)
    } else if (thistype=="linear"){
      gradient[[i]] <- calcGradient.linear(datalist[[i]], resids)
    } else if (thistype=="smooth"){
      gradient[[i]] <- resids 
    } else if (thistype=="textfile" | thistype=="textstring"){
      gradient[[i]] <- calcGradient.text(datalist[[i]], mem.avail.vec[i], resids)
    }
  }
  return(gradient)
}

##### Step & fitted values ##### 

takeStep <- function(params, gradient, datalist, types, resids, sthlambda, smoothlambda, step.size){
  p <- length(types)
  n <- length(resids)
  params.new <- vector("list",p)
  for(i in 1:p){
    thistype = types[i]
    if (thistype=="category"){
      params.new[[i]] <- soft(params[[i]] + step.size * gradient[[i]], sthlambda * step.size / 2)
    } else if (thistype=="linear"){
      params.new[[i]] <- params[[i]] + step.size * gradient[[i]]
    } else if (thistype=="smooth"){
      perm <- order(datalist[[i]])
      invperm = order(perm)
      ord.dat <- datalist[[i]][perm]
      resids.ord <- resids[perm]
      out.step <- rep(0,n)
      # smoothing spline call
      if(!is.loaded("callSS_Fortran")) { dyn.load("Fit.so") }
      junk <- .C("callSS_Fortran", y = as.double(resids.ord), x = as.double(ord.dat), 
                 sy = as.double(out.step), lambda = as.double(smoothlambda), n_point = as.integer(n))
      params.new[[i]] <- params[[i]] + step.size * junk$sy[invperm]
    } else if (thistype=="textfile" | thistype=="textstring"){
      params.new[[i]] <- soft(params[[i]] + step.size * gradient[[i]], sthlambda * step.size / 2)
    }
  }
  return(params.new)
}

initialize.params <- function(types,n.mem.vec,p,n){
  params <- vector("list",p)
  for(i in 1:p){
    thistype = types[i]
    if (thistype=="category"){
      params[[i]] <- rep(0, n.mem.vec[i])
    } else if (thistype=="linear"){
      params[[i]] <- 0
    } else if (thistype=="smooth"){
      params[[i]] <- rep(0,n)
    } else if (thistype=="textfile" | thistype=="textstring"){
      params[[i]] <- rep(0, n.mem.vec[i])
    }
  }
  return(params)
}

fitFromParams <- function(params, types, datalist, p, n){
  yhat <- rep(0, n)
  for(i in 1:p){
    thistype <- types[i]
    if (thistype=="category"){
      thisX <- as.numeric(datalist[[i]])
      thisbeta <- as.numeric(params[[i]])
      for (j in 1:n){ 
        yhat[j] <- yhat[j] + thisbeta[thisX[j]+1]
      }
    } else if (thistype=="linear"){
      thisX <- as.numeric(datalist[[i]])
      thisbeta <- as.numeric(params[[i]])
      yhat <- yhat + thisX * thisbeta 
    } else if (thistype=="smooth"){
      thisX <- as.numeric(datalist[[i]])
      thisbeta <- as.numeric(params[[i]])
      yhat <- yhat + thisbeta 
    } else if (thistype=="textfile" | thistype=="textstring"){
      thisX <- datalist[[i]]
      thisbeta <- as.numeric(params[[i]])
      for (j in 1:n) {
        nwords <- length(thisX[[j]])
        for (k in c(1:nwords)){
          thishash <- thisX[[j]][k]
          yhat[j] <- yhat[j] + thisbeta[thishash+1]
        }
      }
    }
  }
  return(yhat)
}

calcDiff <- function(beta.old.list, beta.new.list){
  all.params.old <- c() 
  all.params.new <- c() 
  for (i in 1:length(beta.old.list)){
    all.params.old <- c(all.params.old, as.numeric(beta.old.list[[i]]))
    all.params.new <- c(all.params.new, as.numeric(beta.new.list[[i]]))
  }
  diff <- sum(all.params.new - all.params.old)
  return(diff)
}

##### Model Fitting ##### 

hashFit <- function(datalist, types, y, n.mem.vec, family=c("binomial","linear"), 
                    sthlambda, smoothlambda, thresh=1e-3, maxit=1e3, step.size=0.001){
  p <- length(datalist)
  n <- length(y)
  yhat <- rep(0, n)
  params <- initialize.params(types, n.mem.vec, p, n)
  
  count <- 0 
  coef.change <- thresh + 1 
  while((coef.change > thresh) & (count < maxit)){
    if (family == "linear") { resids <- y - yhat }
    if (family == "binomial") { resids <- y - expit(yhat) }
    
    gradient <- calcGradient(datalist, n.mem.vec, resids, types)
    params.new <- takeStep(params, gradient, datalist, types, resids, 
                           sthlambda, smoothlambda, step.size)
    yhat <- fitFromParams(params.new, types, datalist, p, n)
    coef.change <- calcDiff(params, params.new)
    count <- count + 1
    params <- params.new 
  }
  return(yhat)
}


##### Cross Validation ##### 
# Credit: Noah Simon

# Get folds for cross validation. 
selectFolds <- function(n, nfolds){
  perFold <- n/nfolds
  reorder <- sample(1:n, replace=FALSE)
  folds <- list()
  last <- 1
  for (i in 1:nfolds){
    folds[[i]] <- reorder[last:(last + perFold - 1)]
    last <- last + perFold
  }
  return(folds)
}

lambdas <- function(X, Z, y, family=c("binomial","linear"), n.sth.lam=11, n.sm.lam=11, n, I, M, N){
  # Soft thresholding lambdas
  xty <- rep(0,M)
  if (family == "binomial") {
    for (j in 1:n) {
      for (i in 1:I) {
        m = sum(N[1:i] + X[j,i])
        xty[m] <- xty[m] + (y[m] - 1/2)
      }
    }
    zty <- t(Z) %*% (y - 1/2) 
  }
  if (family == "linear") {
    for (j in 1:n) {
      for (i in 1:I) {
        m = sum(N[1:i]) + X[j,i]
        xty[m] <- xty[m] + y[m]
      }
    }
    zty <- t(Z) %*% y 
  }
  
  max.sth.lambda <- max(abs(c(xty, zty)), na.rm=TRUE)
  sth.lambdas <- rev(exp(seq(log(0.1), log(max.sth.lambda), length.out=n.sth.lam)))
  
  # Smoothing splines lambdas
  sm.lambdas <- exp(seq(-5,5,length.out=n.sm.lam))
  
  # Matrix of lambdas
  sm.col <- as.numeric(sapply(sm.lambdas, function(x) rep(x, length(sth.lambdas))))
  sth.col <- rep(sth.lambdas, length(sm.lambdas))
  alllambdas <- cbind("smooth"=sm.col, "sth"=sth.col)
  return(alllambdas)
}

testOneFold <- function(X, Z, y, thisfold, all.lambdas, family=c("binomial","linear"), 
                        thresh=1e-3, maxit=1e3, step.size=0.001) {
  # training and test sets 
  train.X <- as.matrix(X[-thisfold,])
  test.X <- as.matrix(X[thisfold,])
  train.Z <- as.matrix(Z[-thisfold,])
  test.Z <- as.matrix(Z[thisfold,])
  train.y <- y[thisfold]
  train.n <- length(train.y)
  test.n <- n - train.n
  orig.beta <- beta
  
  # fit model for each lambda pair on this fold 
  n.lam <- length(all.lambdas[,1])
  res <- as.list(rep(NA, n.lam))
  old.sm.lam <- 0
  test.eta <- matrix(nrow = n, ncol=n.lam)
  for (l in 1:n.lam) tryCatch( {
    sm.lam <- all.lambdas[l,1]
    sth.lam <- all.lambdas[l,2]
    if (sm.lam != old.sm.lam) {
      if (!is.null(orig.beta)) {
        beta <- orig.beta
      } else {
        beta <- rep(0,M)
      }
    }  # Use the updated beta within the same smoothing lambda
    res[[l]] <- modelfit(train.X, train.Z, train.y, sth.lam, sm.lam, family)
    
    beta <- res[[l]]$beta
    old.sm.lam <- sm.lam
    print(paste("Done with lambdas sth", sth.lam, "sm", sm.lam))
    
    # calculate fitted values on test data 
    test.cat <- rep(0, test.n)
    for (j in 1:test.n) { for (i in 1:I) {
      m = sum(N[1:i]) + test.X[j,i]
      test.cat[j] = test.cat[j] + res[[l]]$beta[m]
    } }
    test.spline <- matrix(nrow=test.n, ncol=G)
    if (!is.matrix(train.Z)) {train.Z <- as.matrix(train.Z)}
    for (g in 1:G) { # Linear interpolation 
      test.spline[,g] <- approx( x = train.Z[,g], y = res[[l]]$theta[,g], xout=test.Z[,g], rule=2)$y
    }
    test.eta[,l] <- test.cat + apply(test.spline, 1, function(x) sum(x) )
  }, error=function(e) {print(paste("Lambdas sth", sth.lam, "sm", sm.lam, "failed."))})
  return(test.eta)
}

crossval <- function(X=NULL, Z=NULL, y, family=c("binomial", "linear"), nfolds=10, thresh=1e-3, maxit=1e3, beta=NULL, intercept=TRUE, adaptive=TRUE, step=0.001) {
  # Variable set-up
  vars <- varprep(X, Z, y, intercept)
  X <- vars$X
  Z <- vars$Z
  y <- vars$y
  n <- vars$n
  I <- vars$I
  N <- vars$N
  M <- vars$M
  G <- vars$G
  
  # Set up folds and lambda values to try  
  folds <- selectFolds(n, nfolds)
  all.lambdas <- lambdas(X, Z, y, family, n=n, N=N, I=I, M=M) 
  
  # Try each pair of lambdas on each fold; store fitted values of test set
  all.test.etas <- matrix(0, nrow=1, ncol=length(all.lambdas[,1]))
  for (i in 1:nfolds){
    thisfold.eta <- testOneFold(X, Z, y, folds[[i]], all.lambdas, family, thresh, maxit, beta, intercept, adaptive, step, n, N, I, M, G)
    all.test.etas <- rbind(all.test.etas, thisfold.eta)
  }
  losses <- apply(all.test.etas, 2, FUN = function(x) { loss(x, y, family) })
  best.loss <- which.min(losses) 
  best.lams <- all.lambdas[best.loss,]
  return(best.lams)
}

make.plot <- function(fitted,true,name){
  pdf(paste(name,".pdf",sep=""),w=6,h=6)
  plot(true,fitted,xlab="true",ylab="fitted",main=name)
  dev.off()
}




