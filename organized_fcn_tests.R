## Test gradFit code 
source("coordgrad_withC.R")

#### Set-up X and y #### 

## X-variables 
test.text <- c() 
for (i in formatC(c(1:100), width=3, format="d", flag="0")){
  test.text <- c(test.text, paste("txt_data/text",i,".txt", sep=""))
}
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse","parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
test.linear <- rnorm(100,0,2)
test.smooth <- runif(100,-3,3)


## Y-variables 
y.text <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[1]][[1]][[i]]  ## from only hashing a text variable
  y.text[i] <- y.text[i] + length(this.txt[this.txt < 200])  
}
y.text <- y.text / sd(y.text)

y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="iguana")] <- 2
y.cat <- y.cat/sd(y.cat)

y.linear <- test.linear * 2
y.linear <- y.linear / sd(y.linear)

y.smooth <- sin(test.smooth)
y.smooth <- y.smooth / sd(y.smooth)

#### single variable models ####
### for each test, run the code at the bottom to fit the model and see results

## text 
d <- data.frame(test.text); t <- c("textfile"); m <- c(2000)
y <- y.text
lam <- 1
alph <- rep(1,4)

## categorical
d <- data.frame(test.cat); t <- c("category"); m <- c(500)
y <- y.cat
lam <- 0.5
alph <- rep(1,4)

## linear
d <- data.frame(test.linear); t <- c("linear"); m <- c(NA)
y <- y.linear
lam <- 1
alph <- rep(1,4)

## smooth 
d <- data.frame(test.smooth); t <- c("smooth"); m <- c(NA)
y <- y.smooth
lam <- 5
alph <- rep(1,4)

#### two-variable models #### 
### alphas in order are (linear, smooth, cat, text) 

## Text, Cat 
d <- data.frame(test.text, test.cat); t <- c("textfile","category"); m <- c(2000,500)
y <- y.text+y.cat
lam <- 1
alph <- c(1,1,1,1)

## Text, Linear
d <- data.frame(test.text, test.linear); t <- c("textfile","linear"); m <- c(2000,NA)
y <- y.text+y.linear
lam <- 1
alph <- c(1,1,1,2)

## Text, Smooth
d <- data.frame(test.text, test.smooth); t <- c("textfile","smooth"); m <- c(2000,NA)
y <- y.text+y.smooth
lam <- 1
alph <- c(1,0.1,1,1)

## Category, Linear
d <- data.frame(test.cat, test.linear); t <- c("category","linear"); m <- c(500,NA)
y <- y.cat + y.linear
lam <- 110
alph <- c(1,1,1,1)

## Category, Smooth 
d <- data.frame(test.cat, test.smooth); t <- c("category","smooth"); m <- c(500,NA)
y <- y.cat + y.smooth
lam <- 1
alph <- c(1,0.1,1,1)

## Linear, Smooth
d <- data.frame(test.linear, test.smooth); t <- c("linear","smooth"); m <- c(NA,NA)
y <- y.linear + y.smooth
lam <- 50
alph <- c(1,1,1,1)

#### all four variables #### 
d <- data.frame(test.linear, test.smooth, test.cat, test.text)
t <- c("linear","smooth","category","textfile")
m <- c(NA,NA,500,2000)
y <- y.linear + y.smooth + y.cat + y.text
lam <- 0.5
alph <- c(1,1,0,1)

## text, linear, smooth
d <- data.frame(test.linear, test.smooth, test.text)
t <- c("linear","smooth","textfile")
m <- c(NA,NA,2000)
y <- y.linear + y.smooth + y.text
lam <- 0.5
alph <- c(1,1,1,1)

## text, cat, smooth 
d <- data.frame(test.smooth, test.cat, test.text)
t <- c("smooth","category","textfile")
m <- c(NA,500,2000)
y <- y.smooth + y.cat + y.text
lam <- 0.5
alph <- c(1,1,0,1)

## linear, cat, smooth
d <- data.frame(test.linear, test.smooth, test.cat)
t <- c("linear","smooth","category")
m <- c(NA,NA,500)
y <- y.linear + y.smooth + y.cat
lam <- 0.2
alph <- c(1,1,0,1)

## more than one of some types? 
test.linear2 <- rnorm(100,0,1)
test.cat2 <- sample(c("purple","red","blue","green","magenta","fuschia","ecru","lilac","navy","orchid","pink","brown","gray","black"), 100, replace=TRUE)
test.smooth2 <- rnorm(100,0,3)

y.cat2 <- rep(0,100); y.cat2[which(test.cat2=="purple")] <- 1; y.cat2[which(test.cat2=="magenta")] <- 4; y.cat2[which(test.cat2=="orchid")] <- 2
y.cat2 <- y.cat2/sd(y.cat2)

y.linear2 <- test.linear2 * 5
y.linear2 <- y.linear2 / sd(y.linear2)

y.smooth2 <- test.smooth2^2 
y.smooth2 <- y.smooth2 / sd(y.smooth2)

d <- data.frame(test.linear, test.smooth, test.cat, test.text, test.linear2, test.cat2, test.smooth2)
t <- c("linear","smooth","category","textfile", "linear","category","smooth")
m <- c(NA,NA,500,2000,NA,100,NA)
y <- y.linear + y.smooth + y.cat + y.text + y.linear2 + y.cat2 + y.smooth2
lam <- 2
alph <- c(1,1,1,1)

##### Code for running each test ##### 
hd <- data.prep(d, t, m, tuples=F)
system.time(all.fit <- hashFit(hd, t, params=NULL, y, m, family="gaussian", lambda=lam, alphas=alph))
all.fit$step; all.fit$count
plot(all.fit[[1]] ~ y)




