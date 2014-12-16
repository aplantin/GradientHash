## Test gradFit code 
source("coordgrad_withC.R")

##### Make and Process Test Data ##### 

##### X-vectors ##### 
set.seed(1456)
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse","parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
test.text <- c() 
for (i in formatC(c(1:100), width = 3, format = "d", flag = "0")){
  test.text <- c(test.text, paste("txt_data/text",i,".txt", sep=""))
}
test.linear <- runif(100,0,1)
test.smooth <- runif(100,0,6)

my.df <- data.frame(test.cat, test.linear, test.text, test.smooth)
numMem <- c(100, NA, 1000, NA)
type <- c("category", "linear", "textfile", "smooth")
all.hash <- data.prep(my.df, type, numMem)

##### Y-vectors #####
hash.data <- all.hash
y.cat <- rep(0,100); y.cat[which(test.cat=="alpaca")] <- 3; y.cat[which(test.cat=="llama")] <- 2
y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[3]][[i]]
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt==535])  
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt==220]) 
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt==112]) 
}
y.linear <- test.linear 
y.smooth <- sin(test.smooth)
err <- rnorm(100, 0, 1)
y <- y.cat + y.linear + y.smooth + y.txt + err

##### Test functions ##### 
all.fit <- hashFit(hash.data, type, y, numMem, family="gaussian", sthlambda=2, smoothlambda=0.5, maxit=1e5) 
plot(all.fit[[1]] ~ y)
