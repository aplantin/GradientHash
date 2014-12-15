## Test gradFit code 
source("GradHashFit.R")

##### Make and Process Test Data ##### 
test.cat <- sample(c("cat","dog","fish","turtle","iguana","horse",
                     "parrot","macaw","rabbit","alpaca","llama"),100,replace=TRUE)
hashed.cat <- hashCategories(test.cat, 15) 
test.text <- c() 
for (i in formatC(c(1:100), width = 3, format = "d", flag = "0")){
  test.text <- c(test.text, paste("../txt_data/text",i,".txt", sep=""))
}
test.smooth <- rnorm(100,0,1)
test.df <- data.frame(test.cat, test.text, test.smooth)
hash.data <- data.prep(test.df, c("category","textfile","smooth"), c(15,1000,NA), tuples=FALSE)

y.smooth <- sin(test.smooth)
y.cat <- rep(0,100)
y.cat[which(test.cat=="alpaca")] <- 3
y.cat[which(test.cat=="llama")] <- 2
y.txt <- rep(0,100)
for(i in 1:100){ 
  this.txt <- hash.data[[2]][[i]]
  y.txt[i] <- 2*length(this.txt[this.txt==534]) # happy 
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt==219]) # friendly 
  y.txt[i] <- y.txt[i] + length(this.txt[this.txt==111]) # amazing
}
err <- rnorm(100,0,0.1)
y <- y.smooth + y.cat + y.txt + err

##### Test functions ##### 
calcGradient.category(hash.data[[1]],15,y)
calcGradient.text(hash.data[[2]],1000,y)
calcGradient.linear(rep(1,100),y)
# gradient for smooth is just residuals 

calcGradient(hash.data,c(15,1000,NA),y,c("category","textfile","smooth"))


##### Test Fit with Category & Textfile ##### 
test.df <- test.df <- data.frame(test.cat)
hash.data <- data.prep(test.df, c("category"), c(15))
y <- y.cat
fitted <- hashFit(hash.data, c("category"), y, c(15), family="linear",
        sthlambda=1, smoothlambda=0.1, thresh=1e-3, maxit=1e3, step.size=0.01)

