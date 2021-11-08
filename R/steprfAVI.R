#' @title Select predictive variables for random forest by AVI and accuracy
#' in a stepwise algorithm
#'
#' @description This function is to select predictive variables for random forest
#' by their averaged variable importance (AVI) that is calculated for each model
#'  after excluding the least important variable, and returns the corresponding predictive
#'   accuracy. It is also developed for 'steprf' function.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param rpt iteration number of cross validation.
#' @param predacc "VEcv" for vecv for numerical data, or "ccr" (i.e., correct
#'  classification rate) or "kappa" for categorical data.
#' @param importance importance of predictive variables.
#' @param maxk maxk split value. By default, 4 is used.
#' @param nsim iteration number for 'avi'. By default, 100 is used.
#' @param min.n.var minimum number of predictive variables remained in the final
#' predictive model the default is 2. If 1 is used, then warnings: 'invalid mtry:
#'  reset to within valid range' will be issued, which should be ignored.
#' @param corr.threshold correlation threshold and the defaults value is 0.5.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A list with the following components: 1) variable.removed:
#' variable removed based on AVI, 2) predictive.accuracy: averaged predictive
#'  accuracy of the model after excluding the variable.removed, 3)
#'  delta.accuracy: contribution to accuracy by each variable.removed, and 4)
#'   predictive.accuracy2: predictive accuracy matrix of the model after
#'    excluding the variable.removed for each iteration.
#'
#' @references Li, J. 2013. Predicting the spatial distribution of seabed
#' gravel content using random forest, spatial interpolation methods and
#' their hybrid methods. Pages 394-400  The International Congress on
#' Modelling and Simulation (MODSIM) 2013, Adelaide.
#'
#' Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z., Przeslawski, R.,
#' Radke, L., Howard, F. and Nichol, S. (2017). "Application of random forest,
#'  generalised linear model and their hybrid methods with geostatistical
#'  techniques to count data: Predicting sponge species richness."
#'  Environmental Modelling & Software 97: 112-129.
#'
#' Li, J., Siwabessy, J., Huang, Z., Nichol, S. (2019). "Developing an optimal
#'  spatial predictive model for seabed sand content using machine learning,
#'  geostatistics and their hybrid methods." Geosciences 9 (4):180.
#'
#' Li, J., Siwabessy, J., Tran, M., Huang, Z. and Heap, A., 2014. Predicting
#'  Seabed Hardness Using Random Forest in R. Data Mining Applications with R. Y.
#'   Zhao and Y. Cen. Amsterdam, Elsevier: 299-329.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Smith, S.J., Ellis, N., Pitcher, C.R., 2011. Conditional variable
#' importance in R package extendedForest. R vignette
#' [http://gradientforest.r-forge.r-project.org/Conditional-importance.pdf].
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#'
#' data(petrel)
#' set.seed(1234)
#' steprf1 <- steprfAVI(trainx = petrel[, c(1,2, 6:9)], trainy = petrel[, 5],
#'  rpt = 2, predacc = "VEcv", importance = TRUE, nsim = 3, min.n.var = 2)
#' steprf1
#'
#' steprf2 <- steprfAVI(trainx = hard[, -c(1, 17)], trainy = hard[, 17],
#'  rpt = 2, predacc = "ccr", importance = TRUE, nsim = 3, min.n.var = 2)
#' steprf2
#'
#' #plot steprf1 results
#' library(reshape2)
#' pa1 <- as.data.frame(steprf1$predictive.accuracy2)
#' names(pa1) <- steprf1$variable.removed
#' pa2 <- melt(pa1, id = NULL)
#' names(pa2) <- c("Variable","VEcv")
#' library(lattice)
#' par (font.axis=2, font.lab=2)
#' with(pa2, boxplot(VEcv~Variable, ylab="VEcv (%)", xlab="Predictive variable removed"))
#'
#' barplot(steprf1$delta.accuracy, col = (1:length(steprf1$variable.removed)),
#' names.arg = steprf1$variable.removed, main = "Predictive accuracy vs variable removed",
#' font.main = 4, cex.names=1, font=2, ylab="Increase in VEcv (%)")
#' dev.off()
#' }
#'
#' @export
steprfAVI <- function (trainx, trainy, cv.fold = 10, ntree = 500, rpt = 20, predacc = "VEcv", importance = TRUE, maxk = c(4), nsim = 100, min.n.var = 2, corr.threshold = 0.5, ...) {
  n <- nrow(trainx); p <- ncol(trainx)
  impvar <- NULL; rmvar2 <- NULL
  predictive.accuracy <- NULL; delta.accuracy <- NULL
  #predictive.accuracy2 <- vector((p - min.n.var), mode = "list")
  predictive.accuracy2 <- array(NA, dim = c(rpt, (p - min.n.var)))
  rfcv1 <- NULL
  for (j in 1:rpt) {
  rfcv1[j] <- RFcv2(trainx, trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc)
  }
  predictive.accuracy1 <- mean(rfcv1)
  avi1 <- spm::avi(trainx, trainy, ntree = ntree, nsim = nsim, maxk = maxk, importance = importance, corr.threshold = corr.threshold)
  rmvar <- avi1$impvar[p]
  rmvar2[1] <- avi1$impvar2[p]
  trainx1 <- trainx[, -rmvar]
  for (i in 1:(p - min.n.var)) {
  rfcv1 <- NULL
  for (j in 1:rpt) {
  rfcv1[j] <- RFcv2(trainx1, trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc)
  }
  predictive.accuracy [i] <- mean(rfcv1)
  delta.accuracy[i] <- predictive.accuracy1 - predictive.accuracy [i]
  predictive.accuracy1 <- predictive.accuracy [i]
  predictive.accuracy2[ , i] <- rfcv1
  if (i == p - min.n.var)
    { break
    } else
    {
    avi1 <- spm::avi(trainx1, trainy, ntree = ntree, nsim = nsim, importance = importance)
    rmvar <- avi1$impvar[p-i]
    rmvar2[1+i] <- avi1$impvar2[p-i]
    trainx1 <- trainx1[, -rmvar, drop = FALSE]
    }
  }
  #rmvar2 <- rmvar2 [- length(rmvar2)]
  list(variable.removed = rmvar2, predictive.accuracy = predictive.accuracy,
  delta.accuracy = delta.accuracy,  predictive.accuracy2 =  predictive.accuracy2)
}
