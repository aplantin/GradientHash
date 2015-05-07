## Test gradFit code 
source("coordgrad_withC.R")


#### only text ####
## all 0 when lambda = 1700, alpha = 0.25 
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
all.fit <- hashFit(hash.data, type, params=NULL, y.text, numMem, family="gaussian", lambda=10, 
                   alphas=rep(1,4), thresh=1e-6, maxit=10000, fixed=TRUE, step.size=0.001) 
plot(all.fit[[1]] ~ y.text)

## with majorization 
all.fit <- hashFit(hash.data, type, params=NULL, y.text, numMem, family="gaussian", lambda=10, 
                   alphas=rep(0.25,4), thresh=1e-6, maxit=1000) 
plot(all.fit[[1]] ~ y.text)

## warm start 
warm.fit <- warmstart(hash.data, type, params=NULL, y.text, numMem, 
                      family="gaussian", minLam=4, maxLam=450, nLam=75, alphas=rep(1,4), 
                      thresh=1e-5, maxit=1e6)
plot(warm.fit[[1]] ~ y.text)





#### only categorical ##### 
## lambda=45, alpha = 0.25 ==> all betas = 0 
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
## lambda = 4 (alpha = 1) ==> all betas = 0 
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
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=2, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### linear and categorical #### 
my.df <- data.frame(test.linear, test.cat)
numMem <- c(NA, 1000) 
type <- c("linear", "category")
hash.data <- data.prep(my.df, type, numMem)
y <- y.linear + y.cat
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", 
                   lambda=2, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7, fixed=TRUE, step.size=0.001) 
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
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=rep(0.25,4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)

#### everything but smooth #### 
## recall order of alphas is (linear, smooth, categorical, text)
my.df <- data.frame(test.linear, test.text, test.cat)
numMem <- c(NA, 1000, 1000) 
type <- c("linear","textfile","category")
hash.data <- data.prep(my.df, type, numMem)
y <- y.text + y.cat + y.linear
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=rep(0.25, 4), thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y)


#### All four variable types #####   
my.df <- data.frame(test.cat, test.linear, test.text, test.smooth)
numMem <- c(1000, NA, 1000, NA)
type <- c("category", "linear", "textfile", "smooth")
hash.data <- data.prep(my.df, type, numMem)

### Y-vectors 
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="llama")] <- 2
y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[1]][[3]][[i]]
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt < 200])  
}
y.linear <- test.linear 
y.smooth <- sin(test.smooth)
err <- rnorm(100, 0, 1)


##### Test functions #####
### Without error 
y <- y.cat + y.linear + y.smooth + y.txt
all.fit <- hashFit(hash.data, type, params=NULL, y, numMem, family="gaussian", lambda=1, alphas=c(0.3, 0, 0.3, 0.4), thresh=1e-5, maxit=1e7) 
plot(all.fit[[1]] ~ y)

##### Categorical time comparisons ##### 

### Number of memory locations -- categorical, constant number of categories matter (here, 2 out of 10 categories)
test.time.cat <- function(cat.df, y.cat, mem){
  numMem <- c(mem) 
  type <- c("category")
  system.time(hash.data <- data.prep(my.df, type, numMem))
  t <- system.time(all.fit <- hashFit(hash.data, type, params=NULL, y.cat, numMem, family="gaussian", sthlambda=0, 
                                 smoothlambda=0, thresh=1e-6, maxit=1e7)) ## 26.45 seconds 
  b <- unique(all.fit$betas[[1]]) ## there was a conflict 
  c <- all.fit$count
  return(c(t, b, c))
}

## Testing number of memory locations -- constant number of categories matter (here, 2 out of 10 categories)
cat.df <- data.frame(test.cat)
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="iguana")] <- 2
test.time.cat(cat.df, y.cat, 100)
test.time.cat(cat.df, y.cat, 200)
test.time.cat(cat.df, y.cat, 500)
test.time.cat(cat.df, y.cat, 1000)
test.time.cat(cat.df, y.cat, 10000)

### Number of memory locations -- categorical, constant number of memory locations (100), 
### constant number of categories matter (here, first 2 of 10), 
test.time.cat.n <- function(k, n, mem){  ## k is number of categorical variables, n is number of observations
  test.categories <- data.frame(replicate(k, sample(c(1:mem),10,replace=TRUE)))
  test.dat <- data.frame(apply(test.categories, 2, FUN=function(x) sample(x, 100, replace=TRUE)))
  sig.cats <- test.categories[1:2,]
  y.mat <- matrix(nrow=nrow(test.dat), ncol=ncol(test.dat))
  for(i in 1:k){
    y.mat[,i] <- apply(test.dat, 1, FUN=function(x) {2*(x[i]==test.dat[1,i]) + 3*(x[i]==test.dat[2,i])})
  }
  y.cat <- apply(y.mat, 1, FUN=function(x) sum(x))
  test.time.cat(test.dat, y.cat, 100)
}

test.time.cat.n(1,100,100)
test.time.cat.n(2,100,100)
test.time.cat.n(4,100,100)
test.time.cat.n(6,100,100)
test.time.cat.n(8,100,100)
test.time.cat.n(10,100,100)

### Doesn't appear (in these tests) to add much time with categorical. I wonder if text is the main problem. 

##### Text-file time comparisons; 10% of hash categories matter ##### 
type=c("textfile")
numMem=10000
test.text <- c() 
for (i in formatC(c(1:100), width=3, format="d", flag="0")){
  test.text <- c(test.text, paste("txt_data/text",i,".txt", sep=""))
}
hash.text <- data.prep(data.frame(test.text), type, numMem)

y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.text[[1]][[i]]
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt < numMem/10])  
}


system.time(all.fit <- hashFit(hash.data, type, y.txt, numMem, family="gaussian", sthlambda=0.5, smoothlambda=0.5, 
                               thresh=1e-5, maxit=1e7))
plot(all.fit[[1]] ~ y.txt)

### With R profiler  
Rprof(filename="testprof.out")
## run code for text, above 
Rprof(NULL) 
summaryRprof("testprof.out")$by.total[1:8,]

##### Smooth time comparisons ##### 
type=c("smooth")
numMem=c(NA)
test.smooth <- runif(100,0,6)
hash.smooth <- data.prep(data.frame(test.smooth), type, numMem)
y.smooth <- sapply(test.smooth, FUN=sin)
system.time(all.fit <- hashFit(hash.smooth, type, y.smooth, numMem, family="gaussian", sthlambda=0, 
                               smoothlambda=1, thresh=1e-6, maxit=1e7) )


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


##### How about with the Arabic tweets? ##### 
arabic.txt <- c()
for(i in 1:1000){
  arabic.txt <- c(arabic.txt, readLines(paste("../Twitter/Negative/negative",i,".txt", sep="")))
}
for(i in 1:1000){
  arabic.txt <- c(arabic.txt, readLines(paste("../Twitter/Positive/positive",i,".txt",sep="")))
}
arabic.txt <- arabic.txt[c(1:102,104:115,117:175,177,179,181:183,185,187:188,190,192:2000)] # the omitted files had only ????? 

arabic.nchar <- c() 
for(i in 1:length(arabic.txt)) {
  arabic.nchar <- c(arabic.nchar, nchar(arabic.txt[[i]]))
}

arabic.y <- c(rep(0,991),rep(1,1000)) # 0 is negative, 1 is positive; all missing files were neg. 

# make data frame 
arabic.df <- data.frame(arabic.nchar, arabic.txt)
arabic.hash <- data.prep(arabic.df, c("smooth","textstring"), c(NA, 1000))
arabic.fit <- hashFit(arabic.hash, c("smooth","textstring"), arabic.y, 
                   c(NA, 1000), family="gaussian", sthlambda=2, smoothlambda=0.5, maxit=1e5) 
plot(arabic.fit[[1]] ~ y)
