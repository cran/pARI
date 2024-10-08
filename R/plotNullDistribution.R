#' @title Plot permutation \eqn{p}-values distribution
#' @description Create a plot of permutation-based \eqn{p}-values with corresponding specified critical vectors.
#' @usage plotNullDistribution(P=NULL,family="simes",alpha = 0.05, 
#' path = getwd(), name = "plot", delta = 0,
#' copes=NULL,mask=NULL, alternative = "two.sided", rand = FALSE, B = 1000)
#' @param P Matrix of \eqn{p}-values with dimensions \eqn{m \times B} where \eqn{m} is the number of variables 
#' and \eqn{B} the number of permutations used instead of the data matrix \code{X}. Default to \code{NULL}.
#' @param family String character. Name of the family confidence envelope to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"}, \code{"higher.criticism"}, and \code{"power"}.
#' Default to "simes". If more than one critical vector are considered, it must be a vector. 
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @param path Character string. Path to save the plot. The path does not must end with \code{/}. Default to \code{getwd()}.
#' @param name Character string. The name of file that will be used to save the plot. Default to "plot".
#' @param delta Numeric value. \eqn{\delta} value. Please see the reference below. Default to 0. 
#' If more than one critical vector are considered, \code{delta} must be a vector having length equals to the length of the vector specified in \code{family}.
#' @param copes List of NIfTI file. The list of copes, i.e., contrasts maps, one for each subject used to compute the statistical tests.
#' @param mask NIfTI file or character string. 3D array of logical values (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that none of the voxels have to be excluded.
#' @param alternative Character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param rand Boolean value. Default to \code{FALSE}. If \code{rand = TRUE}, the \eqn{p}-values are computed by \code{rowRanks}.
#' @param B Numeric value. Number of permutations, default to 1000. 
#' @author Angela Andreella
#' @return Save a plot in \code{path} with name specified in \code{name} describing the \eqn{p}-values null distribution with critical value curve and observed \eqn{p}-values in red.
#' @export
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @references Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
#' @examples 
#' \dontrun{
#'db <- simulateData(pi0 = 0.8, m = 100, n = 20, rho = 0)
#'out <- signTest(X = db)
#'pv <- cbind(out$pv, out$pv_H0)
#'plotNullDistribution(P = pv)
#' }
#' 
plotNullDistribution <- function(P=NULL,family="simes",alpha = 0.05, 
                                 path = getwd(), name = "plot", delta = 0,
                                 copes=NULL,mask=NULL, alternative = "two.sided", rand = FALSE, B = 1000){
  
  family_set <- c("simes", "aorc", "beta", "higher.criticism", "power")
  fam_match <- function(x) {match.arg(tolower(x), family_set)}
  alternative_set <- c("two.sided", "greater", "lower")
  alternative <- match.arg(tolower(alternative), alternative_set)
  if(!is.null(family)){family <- unlist(lapply(family, fam_match))}
  if(is.null(copes) & is.null(P)){stop('Please insert pvalues matrix or copes images')}
  
  if(!is.null(copes)){
    
    #check copes
    if(!is.list(copes)){stop("Please insert the list of copes as list class object")}
    
    img_dims <- dim(copes[[1]])
    
    #check mask
    if(!is.null(mask)){
      if(!is.character(mask) && !is.array(mask)){stop("mask must be an array or a path")}
      if(is.character(mask)){mask = readNifti(mask)}
      if(!all(dim(mask) == img_dims)){stop("incompatible dimensions of mask and copes")}
    }else{
      mask <- array(1, img_dims)
    }
    
    #create scores matrix
    img <- array(NA, c(img_dims, length(copes)))
    
    for (sid in 1:length(copes)) {  
      img[,,,sid] <- copes[[sid]]
      
    }
    
    scores <- matrix(img,nrow=(img_dims[1]*img_dims[2]*img_dims[3]),ncol=length(copes))
    scores <- scores[which(mask==1),]
    res <- signTest(X=scores, B = B,alternative = alternative, rand = rand) #variables times number of permutation
    
    P <- cbind(res$pv,res$pv_H0)
    rm(res)
    rm(scores)
    rm(copes)
    rm(img)
  
    
  }

  if(is.null(family)){
    if(is.unsorted(P[,1])){pvalues_ord <- colSortC(P)}else{pvalues_ord <- P}
    
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[,1], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:ncol(pvalues_ord)){
      
      lines(pvalues_ord[,i],col='black',type="l")
      
    }
    lines(pvalues_ord[,1], lwd =2, col= 'red')
    dev.off()
  }else{
    if(is.unsorted(P[,1])){pvalues_ord <- colSortC(P)}else{pvalues_ord <- P}
    
    lcv <- function(family,delta=NULL, cols = "blue"){
      lambdaO <- lambdaOpt(pvalues = P,family=family,alpha=alpha, delta = delta)
      cvO<- criticalVector(pvalues = P, family = family, alpha = alpha, lambda = lambdaO, delta = delta)
      lines(cvO, lwd =2, col= cols)
    }
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    lcvV <- Vectorize(lcv,vectorize.args = c("family", "delta", "cols"))
    cols = rainbow(length(family))
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[,1], type = 'l', col = ' red', lty = "dashed", xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:ncol(pvalues_ord)){
      
      lines(pvalues_ord[,i],col='black',type="l")
      
    }
    lines(pvalues_ord[,1], lwd =2, col= 'red', lty = "dashed")
    #lines(cvO, col= 'blue', lwd =2)
    mapply(lcv, family, delta, cols)
    family <- firstup(family)
    family <- ifelse(family == "Aorc", "AORC", family)
    legend('top',legend=c(sapply(c(1:length(family)), 
                                 function(x) as.expression(bquote(~ .(family[x]) ~ delta == .(delta[x]) ))), 
                          " Observed Pvalues"), col= c(cols, "red"),lwd =2, lty =c(rep("solid", length(family)), "dashed"))
    
    dev.off()
  }
  
  
}
