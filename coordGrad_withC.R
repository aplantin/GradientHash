### Install and load required packages 
library(stringi)
library(hashFunction)
library(Rcpp)
library(inline)
library(data.table)

##### Miscellaneous Support Functions ##### 
# For binomial family calculations
expit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(p) {log(p/(1-p))}

# Soft thresholding
soft <- function(xvec, lambda){
  x.soft <- sign(xvec) * pmax(abs(xvec) - lambda, 0)
  return(x.soft)
}

##### Data Preparation ##### 

# **input: single filename or string, # available memory locations 
# hash the text to those locations by word or tuple
hashText <- function(x, type=c("textfile","textstring"), num.mem, tuples){  
  if (type == "textfile") {
    str <- paste(readLines(x))
  } else { str <- x }
  
  str <- stri_replace_all_regex(stri_trans_tolower(str), '[:punct:]', ' ')
  words <- stri_split_regex(str, ' ', omit_empty=TRUE)[[1]]
  if (tuples == TRUE) {
    nwords <- length(words)
    tups <- c()
    for (i in c(1:(nwords-1))){
      this.tup <- paste(words[i],words[i+1])
      tups <- c(tups, this.tup)
    }
    words <- c(words, tups)  
  }
  words <- sapply(words, function(x){ murmur3.32(x) %% num.mem}) ## hash is indexed from 0 
  return(words)
}

hashCategories <- function(x, num.mem){
  cats <- sapply(x, function(y){ murmur3.32(y) %% num.mem })
  thisidx <- sort(unique(cats)) + 1 
  counts <- rep(0, num.mem)
  counts[thisidx] <- table(cats)
  return(list(cats, counts))
}

scaleLinear <- function(x){
  this.var <- var(x)
  rescale <- x/sqrt(this.var)
  return(rescale)
}

countText <- function(x, numMem){
  counts <- rep(0, numMem)
  for(i in 1:length(x)){
    thistab <- as.data.frame(table(x[[i]]))
    thisidx <- sort(unique(x[[i]]))  
    counts[thisidx + 1] <- counts[thisidx + 1] + thistab$Freq
  }
  return(counts)
}

# take a data.frame (n x p) and a vector of types; output a list of lists 
## types can be: "category", "smooth", "linear", "textfile", "textstring" 
## family = c("binomial","gaussian") 
## n.mem.vec = p-vector of number of memory locations (NAs ok for linear/smooth vars)

data.prep <- function(dataframe, types, n.mem.vec, tuples=FALSE, intercept=TRUE){
  p <- ncol(dataframe) 
  n <- nrow(dataframe)
  if(intercept==TRUE){
    out <- vector("list", (p+1))
    counts <- vector("list", (p+1))
  } else {
    out <- vector("list", (p))
    counts <- vector("list", (p))
  }
 
  if(intercept==TRUE){out[[1]] <- rep(1, n)}
  for (i in c(1:p)){
    thistype <- types[i]
    if (thistype == "category"){
      this.out <- hashCategories(as.character(dataframe[,i]), n.mem.vec[i])
      out[[i+intercept]] <- this.out[[1]]
      counts[[i+intercept]] <- this.out[[2]]
    } else if (thistype == "linear" | thistype == "smooth"){
      out[[i+intercept]] <- scaleLinear(dataframe[,i])
    } else if (thistype == "textfile" | thistype == "textstring"){
      out[[i+intercept]] <- sapply(dataframe[,i], FUN = function(x){hashText(as.character(x), thistype, n.mem.vec[i], tuples)})
      counts[[i+intercept]] <- countText(out[[i+intercept]], n.mem.vec[i])
    } else { (out[[i+intercept]] <- "Not a recognized type")}
  }
  if(intercept==TRUE){
    type <- c("intercept",types)
    nMem <- c(NA, n.mem.vec)
  } else{ 
    type <- types
    nMem <- n.mem.vec
  }
  return(list(out, counts, type, nMem))
}

##### Gradient Calculation ##### 
# Categorical variables 
cppFunction('NumericVector gradCat(IntegerVector hash, int mem, IntegerVector count, NumericVector res) {
  int n = res.size(); 
  NumericVector out(mem); 
  for(int i = 0; i < n; ++i) {
    out[hash[i]] += res[i]/sqrt(count[hash[i]]); 
  }
  return out; 
  }
  ')

# text variables (original input as either file or string) 
gradText = cxxfunction(signature(x='List', Mem='integer', Count='integer', Res='numeric'), plugin='Rcpp', body = '
  Rcpp::List xlist(x);
  int mem = Rcpp::as<int>(Mem); 
  NumericVector out(mem); 
  for(int i = 0; i < mem; i++){
    out[i] = 0; 
  }
  NumericVector res(Res); 
  IntegerVector count(Count); 
  int n = res.size(); 

  for(int i = 0; i < n; i++){
    NumericVector ll = xlist[i];  
    int m = ll.size(); 
    for(int j=0; j < m; j++){
      out[ll[j]] += res[i]/sqrt(count[ll[j]]); 
    }
  }
  return(Rcpp::wrap(out)); 
  ')

# Overall 
calcGradient <- function(data, mem, count, resids, thistype){
  if (thistype=="category"){
    gradient <- gradCat(data, mem, count, resids)
  } else if (thistype=="linear" | thistype=="intercept"){
    gradient <- sum(data*resids) 
  } else if (thistype=="smooth"){
    gradient <- resids
  } else if (thistype=="textfile" | thistype=="textstring"){
    gradient <- gradText(data, mem, count, resids)
  } 
  return(gradient)
}

##### Step & fitted values ##### 

initialize.params <- function(types,n.mem.vec,p,n){
  params <- vector("list",p)
  for(i in 1:p){
    thistype = types[i]
    if (thistype=="category" | thistype=="textfile" | thistype=="textstring"){
      params[[i]] <- rep(0, n.mem.vec[i])
    } else if (thistype=="linear"){
      params[[i]] <- 0
    } else if (thistype=="intercept"){
      params[[i]] <- 0
    } else if (thistype=="smooth"){
      params[[i]] <- rep(0,n)
    } 
  }
  return(params)
}

cppFunction('NumericVector fitCat(IntegerVector hash, NumericVector beta, IntegerVector count, NumericVector y) {
  int n = y.size(); 
  for(int i = 0; i < n; ++i) {
    y[i] += beta[hash[i]]/sqrt(count[hash[i]]);
  }
  return y; 
  }')

fitText = cxxfunction(signature(x='List', Beta='numeric', Count='int', Y='numeric'), plugin='Rcpp', body = '
  Rcpp::List xlist(x);
  NumericVector y(Y); 
  IntegerVector count(Count); 
  NumericVector beta(Beta); 
  int n = y.size(); 
  for(int i = 0; i < n; i++){
    NumericVector ll = xlist[i];  
    int m = ll.size(); 
    for(int j=0; j < m; j++){
      y[i] += beta[ll[j]]/sqrt(count[ll[j]]); 
    }
  }
  return(Rcpp::wrap(y)); 
  ')


fitFromParams <- function(thisbeta, thistype, thisX, count){
  zero <- rep(0, length(thisX))
  if (thistype=="category"){ fitted <- fitCat(thisX, thisbeta, count, zero)
  } else if (thistype=="linear" | thistype=="intercept"){ fitted <- thisX * thisbeta 
  } else if (thistype=="smooth"){ fitted <- thisbeta 
  } else if (thistype=="textfile" | thistype=="textstring"){
    fitted <- fitText(thisX, thisbeta, count, zero)
  }
  return(fitted)
}

takeStep <- function(params, gradient, data, thistype, thislambda, step.size){
  n <- length(data)
  # "beta_{j+1}" 
  if (thistype=="category" | thistype=="textfile" | thistype=="textstring"){
    params.new <- soft(params + step.size * gradient, thislambda*step.size)
  } else if (thistype=="linear" | thistype=="intercept"){
    params.new <- soft(params + step.size * gradient, thislambda*step.size)
  } else if (thistype=="smooth"){
    perm <- order(data)
    invperm <- order(perm)
    ord.dat <- data[perm]
    resids.ord <- gradient[perm]
    out.step <- rep(0,n)
    # smoothing spline call
    if(!is.loaded("callSS_Fortran")) { dyn.load("Fit.so") }
    junk <- .C("callSS_Fortran", y = as.double(resids.ord*step.size + params[perm]), x = as.double(ord.dat), 
               sy = as.double(out.step), lambda = as.double(step.size*thislambda), n_point = as.integer(n))
    params.new <- junk$sy[invperm]
  } 
  return(params.new)
}

##### Model Fitting ##### 
## alphas in order (linear, smooth, categorical, text)
hashFit <- function(dataprep.obj, params=NULL, y, family=c("binomial","gaussian"), 
                    lambda, alphas, thresh=1e-5, maxit=1e6, fixed=FALSE, step.size=1){
  datalist <- dataprep.obj[[1]] 
  datacount <- dataprep.obj[[2]]
  types <- dataprep.obj[[3]]
  n.mem.vec <- dataprep.obj[[4]]
  
  lambda.dt <- data.table(type=c("intercept","linear","smooth","category","textfile", "textstring"), lam=lambda*c(0, alphas, alphas[4]))
  setkey(lambda.dt, type)
  
  p <- length(types)
  n <- length(y)
  
  eta <- matrix(0, nrow=n, ncol=p)
  
  if(is.null(params)){
    params <- initialize.params(types, n.mem.vec, p, n)
  } else{
    for(i in 1:p){
      eta[,i] <- fitFromParams(params[[i]], types[i], datalist[[i]], datacount[[i]])
    }
  }
  eta.new <- eta 
  params.new <- params 
  if(length(step.size) != p){ step <- rep(step.size,p) } else{ step <- step.size }  
  step[1] <- 1/n  ### change for logistic!! 
  ### note that if 1 < k < p step sizes are input, this will do some sort of weird vector-filling thing 
  
  count <- 0 
  coef.change <- thresh + 1 
  yhat <- apply(eta, 1, sum)
  
  while(all(coef.change > thresh, count < maxit)){
    coef.change = 0 
    for(j in 1:p){
      resids <- y - yhat
      this.type <- types[j]
      this.grad <- calcGradient(datalist[[j]], n.mem.vec[j], datacount[[j]], resids, this.type)    
      params.new[[j]] <- takeStep(params[[j]], this.grad, datalist[[j]], this.type, lambda.dt[this.type]$lam, step[j])   ## Just update this variable's params
      eta.new[,j] <- fitFromParams(params.new[[j]], this.type, datalist[[j]], datacount[[j]])
      yhat.new <- yhat - eta[,j] + eta.new[,j]
      if(fixed==FALSE){
        if(this.type!="intercept"){
          while( 0.5*sum((y-yhat.new)^2) > (0.5*sum((y-yhat)^2) - sum((params.new[[j]] - params[[j]])*this.grad) + 
                                              (1/(2*step[j])) * sum((params.new[[j]] - params[[j]])^2)) ){
            step[j] <- step[j] * 0.8 
            params.new[[j]] <- takeStep(params[[j]], this.grad, datalist[[j]], this.type, lambda.dt[this.type]$lam, step[j])   ## Just update this variable's params
            eta.new[,j] <- fitFromParams(params.new[[j]], this.type, datalist[[j]], datacount[[j]])
            yhat.new <- yhat - eta[,j] + eta.new[,j]
          }
        }
      }
      coef.change <- coef.change + sum(abs(params.new[[j]] - params[[j]]))
      yhat <- yhat.new 
      eta <- eta.new 
      params[[j]] <- params.new[[j]] 
    }
    count <- count + 1
  }
  return(list("fitted"=yhat, "betas"=params, "delta"=coef.change, "count"=count, "step"=step))
}


warmstart <- function(dataprep.obj, params=NULL, y, 
                      family=c("binomial","gaussian"), minLam, maxLam, nLam, 
                      alphas, thresh=1e-5, maxit=1e6, fixed=FALSE, step.size=1){
  fail <- c() 
  lambdas <- exp( seq(log(minLam), log(maxLam), length.out=nLam) )
  oldbeta <- hashFit(dataprep.obj, params, y, family, maxLam, alphas, thresh, maxit, fixed, step.size)$betas
  for(i in (2:length(lambdas))) tryCatch( {
    thislam = rev(lambdas)[i]
    new.fit <- hashFit(dataprep.obj, params=oldbeta, y, family, 
                       thislam, alphas, thresh=1e-5, maxit=1e6, fixed, step.size)
    oldbeta <- new.fit$betas
    }, error = function(e) fail <- c(fail, i) )
  if(length(fail)>0){ print("Failed for ", fail) }
  return(new.fit)
}

##### Prediction ##### 
predict <- function(dataprep.obj, fit.params){
  datalist <- dataprep.obj[[1]] 
  datacount <- dataprep.obj[[2]]
  types <- dataprep.obj[[3]]
  n.mem.vec <- dataprep.obj[[4]]
  
  n <- length(datalist[[1]])
  zero <- rep(0, n)
  sumfit <- rep(0, n)
  
  for(i in 1:length(types)){
    thistype <- types[i]
    thisX <- datalist[[i]] 
    thisbeta <- fit.params[[i]] 
    count <- datacount[[i]]
    
    if (thistype=="category"){ fitted <- fitCat(thisX, thisbeta, count, zero)
    } else if (thistype=="linear" | thistype=="intercept"){ fitted <- thisX * thisbeta 
    } else if (thistype=="smooth"){ fitted <- thisbeta 
    } else if (thistype=="textfile" | thistype=="textstring"){
      fitted <- fitText(thisX, thisbeta, count, zero)
    }
    sumfit <- sumfit + fitted 
  }
  return(sumfit)
}

traintest.prep <- function(dat, y, proptest, types, mem, tuples=F, intercept=T){
  ntest <- floor(nrow(dat) * proptest)
  which.test <- sample(nrow(dat), ntest)
  which.train <- subset(c(1:nrow(dat)), !c(1:nrow(dat))%in%which.test)
  
  test.dat <- data.frame(dat[which.test,])
  train.dat <- data.frame(dat[which.train,])
  
  test.y <- y[which.test]
  train.y <- y[which.train] 
  
  test.obj <- data.prep(test.dat, types, mem, tuples=F, intercept=T)
  train.obj <- data.prep(train.dat, types, mem, tuples=F, intercept=T)
  
  return(list(test.obj, test.y, train.obj, train.y))
}

predicterr <- function( test.obj, train.obj, test.y, train.y, lambdas, alphas, thresh=1e-5, maxit=1e5, fixed=F ){
  test.dat <- test.obj[[1]] 
  test.count <- test.obj[[2]]
  test.types <- test.obj[[3]]
  test.mem <- test.obj[[4]]

  lambdas <- rev(sort(lambdas)) 
    
  hf <- hashFit(train.obj, params=NULL, train.y, family="gaussian", lambda=lambdas[1], alphas=alphas, thresh=thresh, maxit=maxit, fixed=fixed)
  test.err <- predict(test.obj, hf$betas) - test.y
  
  if(length(lambdas)>1){
    for(i in 2:length(lambdas)){
      lam = lambdas[i]
      hf <- hashFit(train.obj, params=hf$betas, train.y, family="gaussian", lambda=lam, alphas=alphas, thresh=thresh, maxit=maxit, fixed=fixed)
      test.err <- cbind(test.err, predict(test.obj, hf$betas) - test.y )
    }
  }
  return(test.err)
}



