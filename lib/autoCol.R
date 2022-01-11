#returns rainbow colors in response to an input of characters, factors or numbers
autoCol <- function(x){
  if(!is.factor(x)){x <- as.factor(x)}
  outCol <- rainbow(length(unique(x)))[as.numeric(x)]
  return(outCol)
}