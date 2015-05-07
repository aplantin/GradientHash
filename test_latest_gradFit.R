## Test gradFit code 
source("coordgrad_withC.R")

#### only text ####
## all 0 when lambda = 1200, alpha = 1 
test.text <- c() 
for (i in formatC(c(1:100), width=3, format="d", flag="0")){
  test.text <- c(test.text, paste("txt_data/text",i,".txt", sep=""))
}
my.df <- data.frame(test.text)
numMem <- c(2000)
type <- c("textfile")
hash.data <- data.prep(my.df, type, numMem)
y.text <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[1]][[1]][[i]]
  y.text[i] <- y.text[i] + length(this.txt[this.txt < 200])  
}
y.text <- y.text / sd(y.text)
all.fit <- hashFit(hash.data, type, params=NULL, y.text, numMem, family="gaussian", lambda=10, 
                   alphas=rep(1,4), thresh=1e-6, maxit=10000, fixed=TRUE, step.size=0.001) 
plot(all.fit[[1]] ~ y.text)

## with majorization 
all.fit <- hashFit(hash.data, type, params=NULL, y.text, numMem, family="gaussian", lambda=2, 
                   alphas=rep(1,4), thresh=1e-6, maxit=1000) 
plot(all.fit[[1]] ~ y.text)

## warm start 
warm.fit <- warmstart(hash.data, type, params=NULL, y.text, numMem, 
                      family="gaussian", minLam=2, maxLam=200, nLam=75, alphas=rep(1,4), 
                      thresh=1e-5, maxit=1e6)
plot(warm.fit[[1]] ~ y.text)

#### only categorical ##### 
## lambda=10 ish ==> all betas = 0 
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse","parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
my.df <- data.frame(test.cat)
numMem <- c(500) 
type <- c("category")
hash.data <- data.prep(my.df, type, numMem)
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="iguana")] <- 2
all.fit <- hashFit(hash.data, type, params=NULL, y.cat, numMem, family="gaussian", lambda=0.5, 
                   alphas=rep(1,4), thresh=1e-6, maxit=1000, fixed=TRUE, step.size=0.001) 
plot(all.fit[[1]] ~ y.cat) 

## with majorization 
all.fit <- hashFit(hash.data, type, params=NULL, y.cat, numMem, family="gaussian", lambda=0.5, 
                   alphas=rep(1,4), thresh=1e-6, maxit=1000) 
plot(all.fit[[1]] ~ y.cat) 

## with warm start 
warm.fit <- warmstart(hash.data, type, params=NULL, y.cat, numMem, 
                      family="gaussian", minLam=0.1, maxLam=10, nLam=50, alphas=rep(1,4), 
                      thresh=1e-5, maxit=1e6)
plot(warm.fit[[1]] ~ y.cat)

#### only linear #####
## lambda = 8 ish ==> all betas = 0 
test.linear <- rnorm(100,0,2)
my.df <- data.frame(test.linear)
numMem <- c(NA) 
type <- c("linear")
hash.data <- data.prep(my.df, type, numMem)
y.linear <- test.linear * 2
y.linear <- y.linear / sd(y.linear)
all.fit <- hashFit(hash.data, type, params=NULL, y.linear, numMem, family="gaussian", 
                   lambda=1, alphas=rep(1,4), thresh=1e-5, maxit=10000, fixed=TRUE, step.size=0.01) 
plot(all.fit[[1]] ~ y.linear) 

## with majorization --- all betas = 0 at about lambda = 8  
all.fit <- hashFit(hash.data, type, params=NULL, y.linear, numMem, family="gaussian", 
                   lambda=1, alphas=rep(1,4), thresh=1e-5, maxit=10000) 
plot(all.fit[[1]] ~ y.linear) 

## with warm start 
warm.fit <- warmstart(hash.data, type, params=NULL, y.linear, numMem, 
                      family="gaussian", minLam=0.1, maxLam=10, nLam=50, alphas=rep(1,4), 
                      thresh=1e-5, maxit=1e6)
plot(warm.fit[[1]] ~ y.linear)

#### only smooth ####
## lambda = 150, alpha = 0.25 => completely flat.
test.smooth <- runif(100,0,6)
my.df <- data.frame(test.smooth)
numMem <- c(NA)
type <- c("smooth")
hash.data <- data.prep(my.df, type, numMem)
y.smooth <- sin(test.smooth)
all.fit <- hashFit(hash.data, type, params=NULL, y.smooth, numMem, family="gaussian", 
                   lambda=3, alphas=rep(1,4), thresh=1e-6, maxit=1e7, fixed=TRUE, step.size=0.01) 
plot(all.fit[[1]] ~ test.smooth)

## with majorization 
all.fit <- hashFit(hash.data, type, params=NULL, y.smooth, numMem, family="gaussian", 
                   lambda=3, alphas=rep(1,4), thresh=1e-5, maxit=10000) 
plot(all.fit[[1]] ~ test.smooth) 

## with warm start 
warm.fit <- warmstart(hash.data, type, params=NULL, y.smooth, numMem, 
                      family="gaussian", minLam=3, maxLam=50, nLam=100, alphas=rep(1,4), 
                      thresh=1e-5, maxit=1e6)
plot(warm.fit[[1]] ~ test.smooth)

#### linear and text #### 
my.df <- data.frame(test.linear, test.text)
numMem <- c(NA, 1000) 
type <- c("linear", "textfile")
hash.data <- data.prep(my.df, type, numMem)
y <- y.linear + y.text
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=2, alphas=rep(1,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### linear and categorical #### 
my.df <- data.frame(test.linear, test.cat)
numMem <- c(NA, 1000) 
type <- c("linear", "category")
hash.data <- data.prep(my.df, type, numMem)
y <- y.linear + y.cat
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", 
                   lambda=2, alphas=c(1,1,1,1), thresh=1e-5, maxit=1e4) 
plot(all.fit[[1]] ~ y)

#### linear and smooth #### 
my.df <- data.frame(test.linear, test.smooth)
numMem <- c(NA, NA) 
type <- c("linear","smooth")
hash.data <- data.prep(my.df, type, numMem)
y <- y.linear + y.smooth
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### text and categorical #### 
my.df <- data.frame(test.text, test.cat)
numMem <- c(1000, 1000) 
type <- c("textfile","category")
hash.data <- data.prep(my.df, type, numMem)
y <- y.text + y.cat 
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=10, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### text and smooth #### 
my.df <- data.frame(test.text, test.smooth)
numMem <- c(1000, NA) 
type <- c("textfile","smooth")
hash.data <- data.prep(my.df, type, numMem)
y <- y.text + y.smooth
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### categorical and smooth #### 
my.df <- data.frame(test.cat, test.smooth)
numMem <- c(1000, NA) 
type <- c("category","smooth")
hash.data <- data.prep(my.df, type, numMem)
y <- y.cat + y.smooth
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=c(1,0.25,0,1), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### everything but smooth #### 
## recall order of alphas is (linear, smooth, categorical, text)
my.df <- data.frame(test.linear, test.text, test.cat)
numMem <- c(NA, 1000, 1000) 
type <- c("linear","textfile","category")
hash.data <- data.prep(my.df, type, numMem)
y <- y.text + y.cat + y.linear
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=c(0.5, 1, 0, 0), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### everything but text #### 

#### everything but categorical #### 

#### everything but linear #### 

#### All four variable types #####   
my.df <- data.frame(test.cat, test.linear, test.text, test.smooth)
numMem <- c(1000, NA, 1000, NA)
type <- c("category", "linear", "textfile", "smooth")
hash.data <- data.prep(my.df, type, numMem)

y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="llama")] <- 2
y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[1]][[3]][[i]]
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt < 200])  
}
y.linear <- test.linear 
y.smooth <- sin(test.smooth)
err <- rnorm(100, 0, 1)

## Without error 
y <- y.cat + y.linear + y.smooth + y.txt
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=c(1,1,1,10), thresh=1e-5, maxit=1e7) 
plot(all.fit[[1]] ~ y)


##### Testing cross-validation ##### 
## For speed, will use one categorical and one smooth 

set.seed(13)
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse","parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
test.smooth <- runif(100,0,6)
my.df <- data.frame(test.cat, test.smooth)
numMem <- c(1000, NA) 
type <- c("category", "smooth")
hash.data <- data.prep(my.df, type, numMem)

y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="llama")] <- 2
y.smooth <- sapply(test.smooth, FUN=sin)
err <- rnorm(100,0,1)
y <- y.cat + y.smooth + err 

system.time(all.fit <- hashFit(hash.data, type, y, numMem, family="gaussian", sthlambda=0, 
                               smoothlambda=1, thresh=1e-6, maxit=1e7) )
plot(all.fit[[1]] ~ y)
