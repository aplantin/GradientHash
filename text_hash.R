library(stringi)
# read string character-by-character
# for each character: 
## if space: 
  ## paste this.word into a single string (not vector of characters)
  ## append hash(this.word) to words 
  ## if tuple=TRUE, append hash("last.word this.word") to words 
  ## store this.word as last.word
  ## make this.word an empty string 
## if a letter: make lower case and append to this.word 
## if punctuation: move on to next character (do nothing with this character)




