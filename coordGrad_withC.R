## Last Update: 12.14.14

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

# take a vector of filenames or strings and number of available memory locations 
# and hash the text to those locations by word or tuple
hashText <- function(x, type=c("textfile","textstring"), num.mem, tuples){
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
  if (tuples == TRUE) {
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

# take a data.frame (n x p) and a vector of types; output a list of lists 
## types can be: "category", "smooth", "linear", "textfile", "textstring" 
## family = c("binomial","gaussian") 
## n.mem.vec = p-vector of number of memory locations (NAs ok for linear/smooth vars)

data.prep <- function(dataframe, types, n.mem.vec, tuples=FALSE){
  p <- ncol(dataframe) 
  out <- vector("list", p)
  for (i in c(1:p)){
    thistype <- types[i]
    if (thistype == "category"){
      out[[i]] <- sapply(dataframe[,i], FUN = function(x){hashCategories(as.character(x), n.mem.vec[i])})
    } else if (thistype == "smooth" | thistype == "linear"){
      out[[i]] <- dataframe[,i]
    } else if (thistype == "textfile" | thistype == "textstring"){
      out[[i]] <- sapply(dataframe[,i], FUN = function(x){hashText(as.character(x), thistype, n.mem.vec[i], tuples)})
    } else { (out[[i]] <- "Not a recognized type")}
  }
  return(out)
}

##### Gradient Calculation ##### 
# Categorical variables 
library(Rcpp)
cppFunction('NumericVector gradCat(IntegerVector hash, int mem, NumericVector res) {
  int n = res.size(); 
  NumericVector out(mem); 
  for(int i = 0; i < n; ++i) {
    out[hash[i]] += res[i]; 
  }
  return out; 
  }')

# text variables (original input as either file or string) 
library(inline)
gradText = cxxfunction(signature(x='List', Mem='integer', Res='numeric'), plugin='Rcpp', body = '
  Rcpp::List xlist(x);
  int mem = Rcpp::as<int>(Mem); 
  NumericVector out(mem); 
  for(int i = 0; i < mem; i++){
    out[i] = 0; 
  }
  NumericVector res(Res); 
  int n = res.size(); 

  for(int i = 0; i < n; i++){
    NumericVector ll = xlist[i];  
    int m = ll.size(); 
    for(int j=0; j < m; j++){
      out[ll[j]] += res[i]; 
    }
  }
  return(Rcpp::wrap(out)); 
  ')

# linear variables 
gradLin <- function(x, resids){ 
  x <- as.numeric(x)
  gradient.linear <- t(x) %*% resids  
  return(gradient.linear) 
}

# Overall 
calcGradient <- function(data, mem, resids, thistype){
  if (thistype=="category"){
    gradient <- gradCat(data, mem, resids)
  } else if (thistype=="linear"){
    gradient <- gradLin(data, resids)
  } else if (thistype=="smooth"){
    gradient <- resids 
  } else if (thistype=="textfile" | thistype=="textstring"){
    gradient <- gradText(data, mem, resids)
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
    } else if (thistype=="smooth"){
      params[[i]] <- rep(0,n)
    } 
  }
  return(params)
}

cppFunction('NumericVector fitCat(IntegerVector hash, NumericVector beta, NumericVector y) {
  int n = y.size(); 
  for(int i = 0; i < n; ++i) {
    y[i] += beta[hash[i]];
  }
  return y; 
  }')

fitText = cxxfunction(signature(x='List', B='numeric', Y='numeric'), plugin='Rcpp', body = '
  Rcpp::List xlist(x);
  NumericVector y(Y); 
  NumericVector b(B); 
  int n = y.size(); 
  for(int i = 0; i < n; i++){
    NumericVector ll = xlist[i];  
    int m = ll.size(); 
    for(int j=0; j < m; j++){
      y[i] += b[ll[j]]; 
    }
  }
  return(Rcpp::wrap(y)); 
  ')

fitFromParams <- function(thisbeta, thistype, thisX){
  zero <- rep(0, length(thisX))
  if (thistype=="category"){ fitted <- fitCat(thisX, thisbeta, zero)
  } else if (thistype=="linear"){ fitted <- c(thisX) * thisbeta 
  } else if (thistype=="smooth"){ fitted <- thisbeta 
  } else if (thistype=="textfile" | thistype=="textstring"){
    fitted <- fitText(thisX, thisbeta, zero)
  }
  return(fitted)
}
takeStep <- function(params, gradient, data, thistype, resids, sthlambda, smoothlambda, step.size){
  n <- length(resids)
  # "beta_{j+1}" 
  if (thistype=="category" | thistype=="textfile" | thistype=="textstring"){
    params.new <- soft(params + step.size * gradient, sthlambda*step.size/2)
  } else if (thistype=="linear"){
    params.new <- params + step.size * gradient
  } else if (thistype=="smooth"){
    perm <- order(data)
    invperm <- order(perm)
    ord.dat <- data[perm]
    resids.ord <- resids[perm]
    out.step <- rep(0,n)
    # smoothing spline call
    if(!is.loaded("callSS_Fortran")) { dyn.load("Fit.so") }
    junk <- .C("callSS_Fortran", y = as.double(resids.ord), x = as.double(ord.dat), 
               sy = as.double(out.step), lambda = as.double(smoothlambda), n_point = as.integer(n))
    params.new <- params + step.size * junk$sy[invperm]
  } 
  return(params.new)
}

##### Model Fitting ##### 
hashFit <- function(datalist, types, y, n.mem.vec, family=c("binomial","gaussian"), sthlambda, 
                    smoothlambda, thresh=1e-4, maxit=1e3){
  p <- length(types)
  n <- length(y)
  yhat <- rep(0, n)
  eta <- matrix(0, nrow=n, ncol=p)
  eta.new <- matrix(0, nrow=n, ncol=p)
  
  params <- initialize.params(types, n.mem.vec, p, n)
  params.new <- params 
  step <- rep(1, p)
  
  count <- 0 
  coef.change <- thresh + 1 
  yhat <- apply(eta, 1, sum)
  
  while((coef.change > thresh) & (count < maxit)){
    coef.change <- 0 
    for(j in 1:p){
      if (family == "gaussian") { resids <- y - yhat }
      if (family == "binomial") { resids <- y - expit(yhat) }
      this.type <- types[j]
      this.grad <- calcGradient(datalist[[j]], n.mem.vec[j], resids, this.type)    
      params.new[[j]] <- takeStep(params[[j]], this.grad, datalist[[j]], this.type, resids, 
                                  sthlambda, smoothlambda, step[j])   ## Just update this variable's params
      eta.new[,j] <- fitFromParams(params.new[[j]], this.type, datalist[[j]])
      yhat.new <- yhat - eta[,j] + eta.new[,j]
      while( sum((y-yhat.new)^2) > (sum((y-yhat)^2) - sum((params.new[[j]] - params[[j]])*this.grad) - 
                                      1/(2*step.size) * sum((params.new[[j]] - params[[j]])^2)) ){
        step[j] <- step[j] * 0.8 
        params.new[[j]] <- takeStep(params[[j]], this.grad, datalist[[j]], this.type, resids, 
                                    sthlambda, smoothlambda, step[j])   ## Just update this variable's params
        eta.new[,j] <- fitFromParams(params.new[[j]], this.type, datalist[[j]])
        yhat.new <- yhat - eta[,j] + eta.new[,j]
      }
      coef.change <- coef.change + sum(abs(params.new[[j]] - params[[j]]))
      yhat <- yhat.new 
      eta <- eta.new 
      params[[j]] <- params.new[[j]] 
    }
    count <- count + 1
  }
  return(list(yhat, params, coef.change, count))
}
