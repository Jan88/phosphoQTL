### cScore:                 extracts a combined importance score from a randomForest object.
### input:
#  rf:                      A randomForest-object that contains both permutation importance (PI) and the increase in node purity (RSS).
#                           Random Forests that were generated with importance=TRUE are suited for this.
#  predictorsOfInterest:    Predictors for which the combined score should be computed. As a default all predictors are selected.
#  normalize:               Should the scores be normalized with the trait variance? This makes the importance scores comparable across different traits.
### output:
#                           The output consists of a numeric vector containing a score for each predictor of interest.
#                           The score is the product of RSS and PI. Negative RSS- or PI-scores are set to zero first.
#                           These scores are influenced by correlation between predictors and therefore contain a marker-specific bias.
cScore <- function(rf=NULL,predictorsOfInterest=NULL,normalize=TRUE){
  if(class(rf)!="randomForest"){
    stop("A randomForest object is required as input!")
  }
  if(is.null(predictorsOfInterest)){
    predictorsOfInterest <- 1:nrow(rf$importance)
  }else{
    if((!is.numeric(predictorsOfInterest))|(!is.vector(predictorsOfInterest))){
      stop("Predictors of interest have to be entered as a numeric vector!")
    }
  }
  if(ncol(rf$importance)<2){
    stop("A randomForest object with both permutation importance and increase in node purity (importance=TRUE) is required as input!")
  }
  if(normalize){
    if(is.null(rf$y)|rf$type!="regression"){
      stop("Normalization can only be performed if the randomForest object was created by supervised regression!")
    }
    tVar <- var(rf$y)
    PI <- rf$importance[predictorsOfInterest,1]/tVar
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]/tVar
    RSS[which(RSS<0)] <- 0
  }else{
    PI <- rf$importance[predictorsOfInterest,1]
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]
    RSS[which(RSS<0)] <- 0
  }
  combinedScore <- PI*RSS
  names(combinedScore) <- rownames(rf$importance)[predictorsOfInterest]
  return(combinedScore)
}
