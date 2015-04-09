## Test gradFit code 
source("coordgrad_withC.R")

##### Test Data: Yelp ##### 

### X-vectors 
set.seed(13)
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse","parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
test.text <- c() 
for (i in formatC(c(1:100), width=3, format="d", flag="0")){
  test.text <- c(test.text, paste("txt_data/text",i,".txt", sep=""))
}
test.linear <- runif(100,0,1)
test.smooth <- runif(100,0,6)

### only categorical
my.df <- data.frame(test.cat)
numMem <- c(1000) 
type <- c("category")
hash.data <- data.prep(my.df, type, numMem)
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="iguana")] <- 2
all.fit <- hashFit(hash.data, type, y.cat, numMem, family="gaussian", sthlambda=0, smoothlambda=0, thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y.cat)

### only linear 
my.df <- data.frame(test.linear)
numMem <- c(NA) 
type <- c("linear")
hash.data <- data.prep(my.df, type, numMem)
y.linear <- test.linear * 2 
all.fit <- hashFit(hash.data, type, y.linear, numMem, family="gaussian", sthlambda=0, smoothlambda=1, thresh=1e-6, maxit=1e7) 
plot(all.fit[[1]] ~ y.linear)

##### All four variable types #####   
my.df <- data.frame(test.cat, test.linear, test.text, test.smooth)
numMem <- c(1000, NA, 1000, NA)
type <- c("category", "linear", "textfile", "smooth")
hash.data <- data.prep(my.df, type, numMem)

### Y-vectors 
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="llama")] <- 2
y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[3]][[i]]
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt < 200])  
}
y.linear <- test.linear 
y.smooth <- sin(test.smooth)
err <- rnorm(100, 0, 1)


##### Test functions #####
### Without error 
y <- y.cat + y.linear + y.smooth + y.txt
all.fit <- hashFit(hash.data, type, y, numMem, family="gaussian", sthlambda=0, smoothlambda=0.5, thresh=1e-5, maxit=1e7) 
plot(all.fit[[1]] ~ y)

### And with error 
y <- y.cat + y.linear + y.smooth + y.txt + err
all.fit <- hashFit(hash.data, type, y, numMem, family="gaussian", sthlambda=1, smoothlambda=0.5, thresh=1e-5, maxit=1e7) 
plot(all.fit[[1]] ~ y)

##### Categorical time comparisons ##### 

### Number of memory locations -- categorical, constant number of categories matter (here, 2 out of 10 categories)
test.time.cat <- function(cat.df, y.cat, mem){
  numMem <- c(mem) 
  type <- c("category")
  system.time(hash.data <- data.prep(my.df, type, numMem))
  t <- system.time(all.fit <- hashFit(hash.data, type, y.cat, numMem, family="gaussian", sthlambda=0, 
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
