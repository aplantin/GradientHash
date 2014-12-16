grad.hash <- function(X, y, thresh=1e-5, maxit=1e3, beta=NULL, intercept=FALSE) {
  if (class(X) != "matrix") {
    if (is.null(ncol(X))) {X <- as.matrix(as.numeric(X))}
    else if (ncol(X) > 0) {X <- as.matrix(X)}
  }
  if(intercept == TRUE) { X <- cbind(1, X) }
  
  step.size = 1
  
  I = ncol(X)
  n = length(y) 
  
  if(length(beta) == 0) {
    N = rep(0, I + 1)
    for (i in 1:I) {N[i+1] = length(unique(X[,i]))}
    M = sum(N)
    beta <- rep(0, M)
  }
  count = 0
  diff = thresh + 1
  eta <- rep(0, n)
  
  # Use gradient descent algorithm and after convergence, return betas
  while((count < maxit) & (diff > thresh)) { 
    
    resids <- y - eta 
    # Gradients should be calculated based on all residuals with the appropriate X in that category 
    grad <- rep(0, M)
    for(j in 1:n) { 
      for(i in 1:I){
        m = sum(N[1:i]) + X[j,i]
        grad[m] <- grad[m] + resids[j]
      }
    }
    
    # Calculate updated betas; iterate until appropriate step size is reached 
    beta.new <- rep(0, M)
    
    for (m in 1:M) {
      beta.new[m] <- beta[m] + step.size * grad[m] 
    }
    
    eta.new <- rep(0, n)
    for (j in 1:n) {
      for (i in 1:I) {
        m = sum(N[1:i]) + X[j,i]  # mth category 
        eta.new[j] = eta.new[j] + beta.new[m]  # Sum of betas corresponding to categories of X for this observation
      }
    }
    
    diffs <- abs(beta.new - beta)
    while( sum((y-eta.new)^2) > sum((y-eta)^2) - sum(beta * diffs) - 1/(2*step.size) * sum(diffs^2) ){
      step.size <- step.size * 0.9 
      beta.new <- rep(0, M)
      for (m in 1:M) {
        beta.new[m] <- beta[m] + step.size * grad[m] 
      }
      
      eta.new <- rep(0, n)
      for (j in 1:n) {
        for (i in 1:I) {
          m = sum(N[1:i]) + X[j,i]  # mth category 
          eta.new[j] = eta.new[j] + beta.new[m]  # Sum of betas corresponding to categories of X for this observation
        }
      }
      diffs <- beta.new - beta
    }
    
    diff <- sum(diffs)
    beta <- beta.new
    eta <- eta.new
    count <- count + 1 
  }
  return(list(beta=beta, count=count)) 
}


##### 

hashFit <- function(datalist, types, y, n.mem.vec, family=c("binomial","gaussian"), 
                    sthlambda, smoothlambda, thresh=1e-4, maxit=1e3, adaptive=FALSE, step.size=0.00001){
  p <- length(types)
  n <- length(y)
  yhat <- rep(0, n)
  eta <- matrix(0, nrow=n, ncol=p)
  eta.new <- matrix(0, nrow=n, ncol=p)
  
  params <- initialize.params(types, n.mem.vec, p, n)
  if (adaptive == TRUE) { step <- rep(1, p) } else { step <- rep(step.size, p) }
  
  count <- 0 
  coef.change <- thresh + 1 
  
  while((coef.change > thresh) & (count < maxit)){
    coef.change <- 0 
    yhat <- apply(eta, 1, sum)
    for(j in 1:p){
      if (family == "gaussian") { resids <- y - yhat }
      if (family == "binomial") { resids <- y - expit(yhat) }
      this.type <- types[j]
      this.step <- step[j]
      this.grad <- calcGradient(datalist[[j]], n.mem.vec[j], resids, this.type)    
      params.new[[j]] <- takeStep(params[[j]], this.grad, datalist[[j]], this.type, resids, 
                                  sthlambda, smoothlambda, this.step)   ## Just update this variable's params
      eta.new[,j] <- fitFromParams(params.new[[j]], this.type, datalist[[j]])
      yhat <- yhat - eta[,j] + eta.new[,j]
      coef.change <- coef.change + sum(abs(params.new[[j]] - params[[j]]))
    }
    count <- count + 1
    params <- params.new 
    eta <- eta.new 
  }
  return(list(yhat, params, coef.change, count))
}

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




